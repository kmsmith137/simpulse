import click
from datetime import datetime
import json
import logging
import numpy as np
import os
import requests
import time

import simpulse
from ch_frb_l1.rpc_client import RpcClient

# Logging Config
LOGGING_CONFIG = {}
logging_format = "[%(asctime)s] %(process)d-%(levelname)s "
logging_format += "%(module)s::%(funcName)s():l%(lineno)d: "
logging_format += "%(message)s"

# Configure Logging
logging.basicConfig(format=logging_format, level=logging.DEBUG)
log = logging.getLogger()

# Define the dispersion constant used for pulsar studies (Manchester & Taylor 1972)
# NOTE: should probably this import from frb_common, but want to keep container light
k_DM = 1.0 / 2.41e-4


####################
# HELPER FUNCTIONS #
####################


def gaussian(nfreq, f_c, f_low, f_hi, fwhm):
    """
    Simple function to generate a gaussian array which will
    modulate the frequency component of an injected pulse
    with a gaussian envelope

    INPUT:
    ------
    
    nfreq : int
        The number of frequency channels in the pulse
    f_c : float
        The central frequency of the gaussian
    f_low : float
        The lowest frequency (in MHz) of the pulse (eg: 400 MHz for CHIME)
    f_hi : float
        The highest frequency (in MHz) of the pulse (eg: 800 MHz for CHIME)
    fwhm : float
        The full width at half-maximum for the Gaussian envelope (in MHz)
    
    OUTPUT:
    -------

    gaussian_modulation: numpy array of np.float32
        Numpy array containing a gaussian profile normalized to an amplitude of 1
    """
    freq = np.linspace(f_low, f_hi, nfreq, dtype=np.float32)
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    exp = np.exp(-1 / 2 * ((freq - f_c) / sigma) ** 2)
    # TODO: properly normalizing?
    # a = 1/(np.sum(exp))
    a = 1
    gaussian_modulation = a * exp
    return gaussian_modulation


def get_frame0_time(url="http://carillon.chime:54321/get-frame0-time"):
    """
    Make a request to recieve the GPS time for frame number 0

    INPUT:
    ------

    url : string
        The url from which to request the frame0 ctime

    OUTPUT:
    -------

    frame0_ctime : np array of np.datetime64
        The GPS time at frame 0, given as a us precision datetime in a numpy array
    """
    resp = requests.get(url)
    tinfo = resp.json()
    frame0_time = np.array(tinfo["frame0_ctime"] * 1e6, dtype="int").astype("<M8[us]")
    return frame0_time


def timestamp2fpga(time, frame0_ctime):
    """
    Returns the fpga frame number for a given timestamp and frame0 ctime.

    INPUT:
    ------

    time : datetime
        The timestamp to turn to an FPGA frame number

    frame0_ctime : datetime
        The frame0 ctime with which to correct the timestamp

    OUTPUT:
    -------

    fpga_time : int
        The FPGA number corresponding to `time`
    """
    fpga_time = int((time - frame0_ctime).astype(int) / 2.56)
    return fpga_time


########
# MAIN #
########


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
    default=400.1953125,
    help="Frequency value of the bottom channel in the generated array",
)
@click.option(
    "--freq_hi_MHz",
    "freq_hi_MHz",
    required=True,
    type=float,
    default=800.1953125,
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
    "--t_injection",
    "t_injection",
    required=True,
    type=str,
    help="A date string in the form in ISO format (`%y-%m-%dT%H:%M:%SZ`) representing the injection time referenced to infinite frequency",
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
    "--frb_master_base_url",
    "frb_master_base_url",
    required=True,
    type=str,
    default="http://frb-vsop.chime:8001",
    help="The url used to access frb-master",
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
    t_injection,
    gaussian_central_freq,
    gaussian_fwhm,
    frb_master_base_url,
    unique_id,
):
    log.info("Beginning generate_pulse routine...")

    # Seconds per DM unit of delay (across the frequency band)
    dmdt = k_DM * (1 / freq_lo_MHz ** 2 - 1 / freq_hi_MHz ** 2)

    # Generate a pulse using simpulse and the specified params
    t_inf = 0
    pulse = simpulse.single_pulse(
        nt, nfreq, freq_lo_MHz, freq_hi_MHz, dm, sm, width, fluence, spindex, t_inf
    )

    if gaussian_central_freq and gaussian_fwhm:
        # Modulate the frequency spectrum by a gaussian profile
        log.info("Generating gaussian profile...")
        gaussian_modulation = gaussian(
            nfreq, gaussian_central_freq, freq_lo_MHz, freq_hi_MHz, gaussian_fwhm
        )
    else:
        log.info("No gaussian component specified, continuing...")
        gaussian_modulation = np.ones(nfreq, dtype=np.float32)

    # Retrieve the sensitivities from the beam model and modulate the frequency
    # spectrum once again
    # TODO: make date more flexible? What does date even end up changing?
    log.info("Retrieving the sensitivity from the beam model...")
    frb_master_url = frb_master_base_url + "/v1/code/beam-model/get-sensitivity"
    log.info("frb_master_url: {}".format(frb_master_url))
    payload = {"ra": ra, "dec": dec, "date": t_injection, "beam": beam_no}
    resp = requests.post(frb_master_url, json=payload)
    sensitivities = np.float32(np.array(resp.json()["sensitivities"]))
    log.info("Modulating the gaussian profile by the beam model sensitivity...")
    freq_modulation = gaussian_modulation * sensitivities

    # Get the tcp server address for the given beam_no from frb-master
    log.info("Retrieving rpc_server address for beam {}...".format(beam_no))
    frb_master_url = frb_master_base_url + "/v1/parameters/get-beam-info/{}".format(
        beam_no
    )
    resp = requests.get(frb_master_url)
    rpc_server = resp.json()["rpc_server"]
    if rpc_server.endswith("5555"):
        rpc_server = rpc_server.replace("5555", "5556")
    log.info("rpc_server: {}".format(rpc_server))

    # Create an RpcClient object which will send a pulse injection request to L1
    servers = {"injections-test": rpc_server}
    client = RpcClient(servers=servers)

    # Determine the fpga frame number for the injection
    date_format = "%Y-%m-%dT%H:%M:%S.%fZ"
    t_injection = np.datetime64(datetime.strptime(t_injection, date_format), "us")
    frame0_ctime = get_frame0_time()
    fpga0 = timestamp2fpga(t_injection, frame0_ctime)
    log.info("fpga0: {}".format(fpga0))

    # Make the request to the RpcClient to inject the pulse at fpga0
    # NOTE: offsetting beam_no by 10000 as that is currently how we verify
    # injection beams
    # TODO: modify RpcClient to take in a frequency modulation array when
    # injecting a pulse
    injection_beam = beam_no + 10000
    resp = client.inject_single_pulse(beam_no, pulse, fpga0, wait=True, nfreq=nfreq)
    log.info("Injection results: {}".format(resp))


if __name__ == "__main__":
    main()
