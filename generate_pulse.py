import click
from datetime import datetime
import json
import logging
import numpy as np
import os
import requests
import time

import simpulse

# Logging Config
LOGGING_CONFIG = {}
logging_format = "[%(asctime)s] %(process)d-%(levelname)s "
logging_format += "%(module)s::%(funcName)s():l%(lineno)d: "
logging_format += "%(message)s"

# Configure Logging
logging.basicConfig(format=logging_format, level=logging.DEBUG)
log = logging.getLogger()


def gaussian(arr, x_0, f_low, f_hi, fwhm):
    assert len(arr.shape) == 2, "Error: Input array must be 2D!"
    x = np.linspace(f_low, f_hi, arr.shape[0])
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    exp = np.exp(-1 / 2 * ((x - x_0) / sigma) ** 2)
    # TODO: properly normalizing?
    # a = 1/(np.sum(exp))
    a = 1
    # Do some transposing so multiplication works along frequency axis
    return (arr.T * a * exp).T


@click.command()
@click.option(
    "--nt",
    "nt",
    required=True,
    type=int,
    default=1024,
    help="Number of time samples in the generated pulse array",
)
@click.option(
    "--nfreq",
    "nfreq",
    required=True,
    type=int,
    default=16384,
    help="Number of (equally spaced) frequency samples in the generated pulse array",
)
# TODO: should default f_lo/f_hi be 400.1x/800.1x?
@click.option(
    "--freq_lo_MHz",
    "freq_lo_MHz",
    required=True,
    type=float,
    default=400,
    help="Frequency value of the bottom channel in the generated array",
)
@click.option(
    "--freq_hi_MHz",
    "freq_hi_MHz",
    required=True,
    type=float,
    default=800,
    help="Frequency value of the bottom channel in the generated array",
)
@click.option(
    "--ra",
    "ra",
    required=True,
    type=float,
    help="J2000 right ascension in decimal degrees",
)
@click.option(
    "--dec",
    "dec",
    required=True,
    type=float,
    help="J2000 declination in decimal degrees",
)
@click.option(
    "--beam_no", "beam_no", required=True, type=int, help="Beam number to inject into"
)
@click.option(
    "--dm",
    "dm",
    required=True,
    type=float,
    default=0,
    help="DM of the pulse to be generated",
)
@click.option(
    "--sm",
    "sm",
    required=True,
    type=float,
    default=0,
    help="SM (scattering time in ms at 1GHz) of the pulse to be generated",
)
@click.option(
    "--intrinsic_width",
    "width",
    required=True,
    type=float,
    default=0.001,
    help="Width of the pulse to be generated, in seconds",
)
@click.option(
    "--fluence",
    "fluence",
    required=True,
    type=float,
    # TODO: figure out fluence units
    help="Fluence of the pulse to be generated in units * s",
)
@click.option(
    "--spectral_index",
    "spindex",
    required=True,
    type=float,
    default=0,
    help="Spectral index of the pulse to be generated",
)
@click.option(
    "--t_arrival",
    "t_arrival",
    required=True,
    type=int,
    default=0,
    help="Arrival time of the pulse in seconds as freq -> infinity relative to an arbitrary origin",
)
@click.option(
    "--gaussian_central_freq",
    "gaussian_central_freq",
    required=False,
    type=float,
    default=None,
    help="The frequency at which the gaussian peaks (in MHz) for a gaussian spectrum modulation",
)
@click.option(
    "--gaussian_fwhm",
    "gaussian_fwhm",
    required=False,
    type=float,
    default=None,
    help="The FWHM (in MHz) for a gaussian spectrum modulation",
)
@click.option(
    "--frb_master_url",
    "frb_master_url",
    required=True,
    type=str,
    default="http://frb-vsop.chime:8001/v1/code/beam-model/get-sensitivity",
    help="The url used to retrieve the sensitivities from the beam model",
)
@click.option(
    "--distributor_url",
    "distributor_url",
    required=True,
    type=str,
    default="http://frb-vsop.chime:8002/distributor/work/mimic",
    help="The url used to add work to the mimic distributor",
)
@click.option(
    "--unique_id",
    "unique_id",
    required=True,
    type=str,
    help="Unique id used to identify the given injection",
)
def main(
    nt,
    nfreq,
    freq_lo_MHz,
    freq_hi_MHz,
    ra,
    dec,
    beam_no,
    dm,
    sm,
    width,
    fluence,
    spindex,
    t_arrival,
    gaussian_central_freq,
    gaussian_fwhm,
    frb_master_url,
    distributor_url,
    unique_id,
):
    # Generate a pulse using simpulse and the specified params
    pulse = simpulse.single_pulse(
        nt, nfreq, freq_lo_MHz, freq_hi_MHz, dm, sm, width, fluence, spindex, t_arrival
    )

    # Figure out what the array size needs to be to contain the dispersed pulse
    # TODO: should sampling time be 0.0009xx?
    dt_sample = 0.001  # in s
    t_end = pulse.get_endpoints()[1]
    n_chunks = int(t_end / (nt * dt_sample) + 1)
    # Create the data array as a memory map, so that the data can be streamed
    # to the distributor as a binary file in-memory
    log.info("Creating memory map...")
    data = np.memmap("pulse", mode="w+", shape=(nfreq, nt * n_chunks), dtype=np.float32)
    # TODO: should I make n chunks, or just one big intensity array?
    chunk_t0 = 0
    chunk_t1 = n_chunks * nt * dt_sample
    log.info("Adding pulse to timestream...")
    pulse.add_to_timestream(data, chunk_t0, chunk_t1, freq_hi_to_lo=True)
    if gaussian_central_freq and gaussian_fwhm:
        # Modulate the frequency spectrum by a gaussian profile
        log.info("Modulating pulse by a gaussian profile...")
        data = gaussian(
            data, gaussian_central_freq, freq_lo_MHz, freq_hi_MHz, gaussian_fwhm
        )
    # Retrieve the sensitivities from the beam model and modulate the frequency
    # spectrum once again
    # TODO: make date more flexible? What does date even end up changing?
    date = datetime.now().strftime("%Y-%m-%d")
    payload = {"ra": ra, "dec": dec, "date": date, "beam": beam_no}
    resp = requests.post(frb_master_url, json=payload)
    sensitivities = np.array(resp.json()["sensitivities"])
    data = (data.T * sensitivities).T
    # Make a dictionary object to send our pulse array as a list, marked by id
    # (because python lists are JSON encodable, but np arrays are not)
    # and make a request to add work to the mimic distributor
    # (not using chime_frb_api as this container is python 2.7)
    # payload = {"work": [json.dumps(data.tolist())]}
    # log.info("Sending pulse data to the distributor...")
    # log.debug("distributor_url: {}".format(distributor_url))
    # resp = requests.post(distributor_url, json=payload)

    # log.debug(resp.raise_for_status())
    # log.debug(resp.json())


if __name__ == "__main__":
    main()
