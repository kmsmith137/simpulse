# Set non-interactice plotting backend
import matplotlib
matplotlib.use("Agg")

import click
from datetime import datetime
import json
import logging
import os
import requests
import time
from flask import Flask, request

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

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

# Initialize the Flask application
app = Flask(__name__)

# TODO: testing with a single beam, so forcing the rpc client to the beam that's
# running duplication code now for testing
servers = {"injections-test": "tcp://10.6.205.12:5556"}
client = RpcClient(servers=servers)

# Set the location to save plots
plot_dir = "/frb-archiver/frb-injections"


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
    gaussian_modulation = np.exp(-1 / 2 * ((freq - f_c) / sigma) ** 2)
    return gaussian_modulation


def empirical_spectral_model(nfreq, f_low, f_hi, spindex, running):
    """
    Simple function to generate an array which can modulate the
    spectrum of an injected pulse using the empirical FRB model
    with spectral index and spectral running from the morphologies
    working group

    INPUT:
    ------

    nfreq : int
        The number of frequency channels in the pulse
    f_low : float
        The lowest frequency (in MHz) of the pulse (eg: 400 MHz for CHIME)
    f_hi : float
        The highest frequency (in MHz) of the pulse (eg: 800 MHz for CHIME)
    spindex : float
        The spectral index of the pulse to be used in the empirical model
    running : float
        The spectral running of the pulse to be used in the empirical model

    OUTPUT:
    -------

    empirical_modulation: numpy array of np.float32
        Numpy array containing the spectral modulation of the model
    """
    freq = np.linspace(f_low, f_hi, nfreq, dtype=np.float32)
    empirical_modulation = (freq / f_low) ** (spindex + running * np.log(freq / f_low))
    return empirical_modulation

def make_pulse_plot(pulse, spectral_modulation=None, fn=None):
    """
    A function to create a plot of the pulse to be injected. The plot
    is created by generating a dedispersed version of the pulse, and
    injecting it into a nfreq by nt array of zeroes.

    INPUT:
    ------

    pulse : simpulse object
        The simpulse object to be injected
    spectral_modulation : numpy array, optional
        A numpy array to modulate the spectrum of the pulse. None by default
    fn : str, optional
        The name of the image to be saved, None by default
    """

    # Set up the axes for plotting
    fig = plt.figure(figsize=(10, 7))
    gs = gridspec.GridSpec(
        2, 2, width_ratios=[3, 1], height_ratios=[1, 3], hspace=0.0, wspace=0.0
    )
    ax_im = plt.subplot(gs[2])
    ax_ts = plt.subplot(gs[0])
    ax_spec = plt.subplot(gs[3])
    for ax in [ax_ts, ax_spec]:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("")
        ax.set_ylabel("")
    ax_im.tick_params(axis="both", which="major", labelsize=14)

    # Set a dm 0 pulse with an undispersed arrival time of 0.5 to appear in the middle of a diagnostic plot
    dm_0 = simpulse.single_pulse(
        pulse.pulse_nt,
        pulse.nfreq,
        pulse.freq_lo_MHz,
        pulse.freq_hi_MHz,
        0.0,
        pulse.sm,
        pulse.intrinsic_width,
        pulse.fluence,
        pulse.spectral_index,
        0.5,
    )
    data = np.zeros((pulse.nfreq, pulse.pulse_nt))
    dm_0.add_to_timestream(data, 0, 1)

    # Show a cut of the data depending on how expansive the pulse is in time
    max_nt = int(max(1024.0 * 6.0 / 10.0, np.where(np.sum(data, axis=0) != 0)[0].max()))
    min_nt = int(1024.0 * 4.0 / 10.0)
    # Reference time to the pulse time, and put in ms
    max_t = max_nt / 1024. * 1000. - 500.
    min_t = min_nt / 1024. * 1000. - 500.
    data_cut = data[:, min_nt:max_nt]
    if spectral_modulation is not None:
        old_sum = data_cut.sum()
        data_cut_copy = data_cut.copy()
        new_data = data_cut_copy * spectral_modulation[:, np.newaxis]
        new_sum = new_data.sum()
        # Conserve fluence
        factor = old_sum / new_sum
        data_cut = new_data * factor

    # Set values for the time series/spectrum
    x_times = np.linspace(min_t, max_t, data_cut.shape[1])
    y_freqs = np.linspace(pulse.freq_lo_MHz, pulse.freq_hi_MHz, data_cut.shape[0])
    data_ts = data_cut.sum(axis=0)
    data_spec = data_cut.sum(axis=1)

    # Plot the data
    ax_im.imshow(
        data_cut,
        aspect="auto",
        extent=[min_t, max_t, pulse.freq_lo_MHz, pulse.freq_hi_MHz],
        origin="lower",
    )
    ax_ts.plot(x_times, data_ts)
    ax_spec.plot(data_spec, y_freqs)

    # Set labels/axes
    ax_ts.set_xlim(min_t, max_t)
    ax_ts.set_ylim(0, data_ts.max() * 1.1)
    ax_spec.set_xlim(0, data_spec.max() * 1.1)
    ax_spec.set_ylim(pulse.freq_lo_MHz, pulse.freq_hi_MHz)
    ax_im.set_xlabel("Time (ms)", fontsize=16)
    ax_im.set_ylabel("Frequency (MHz)", fontsize=16)

    # Save the plot
    if fn is not None:
        plt.savefig(fn)


###################
# PULSE INJECTION #
###################


# An endpoint to inject a pulse into a CHIME/FRB L1 node
@app.route("/inject-pulse", methods=["POST"])
def inject_pulse():
    log.info("Beginning generate_pulse routine...")
    # Unpack the data from the injection service
    data = request.get_json()
    nt = data["nt"]
    nfreq = data["nfreq"]
    freq_lo_MHz = data["freq_lo_MHz"]
    freq_hi_MHz = data["freq_hi_MHz"]
    ra = data["ra"]
    dec = data["dec"]
    beam_no = data["beam_no"]
    dm = data["dm"]
    sm = data["tau_1_ghz"]
    width = data["pulse_width_ms"] / 1000
    fluence = data["fluence"]
    spindex = data["spindex"]
    running = data["running"]
    t_injection = data["injection_time"]
    gaussian_central_freq = data["gaussian_central_freq"]
    gaussian_fwhm = data["gaussian_fwhm"]
    frb_master_base_url = data["frb_master_base_url"]
    fpga0 = data["timestamp_fpga_injection"]
    plot = data["plot"]

    # Make some basic assertions
    assert (
        0 < beam_no < 255
        or 1000 < beam_no < 1255
        or 2000 < beam_no < 2255
        or 3000 < beam_no < 3255
    ), "ERROR: invalid beam_no!"
    t_injection_dt = datetime.strptime(t_injection, "%Y-%m-%dT%H:%M:%S.%fZ")
    assert (
        t_injection_dt > datetime.utcnow()
    ), "ERROR: t_injection must not be in the past!"

    # Seconds per DM unit of delay (across the frequency band)
    dmdt = k_DM * (1 / freq_lo_MHz ** 2 - 1 / freq_hi_MHz ** 2)

    # Generate a pulse using simpulse and the specified params
    t_inf = 0
    pulse = simpulse.single_pulse(
        nt, nfreq, freq_lo_MHz, freq_hi_MHz, dm, sm, width, fluence, spindex, t_inf
    )

    if gaussian_central_freq and gaussian_fwhm:
        # Modulate the frequency spectrum by a gaussian profile
        log.info("Generating spectral gaussian profile...")
        gaussian_modulation = gaussian(
            nfreq, gaussian_central_freq, freq_lo_MHz, freq_hi_MHz, gaussian_fwhm
        )
    else:
        log.info("No gaussian component specified, continuing...")
        gaussian_modulation = None

    # Retrieve the sensitivities from the beam model and modulate the frequency
    # spectrum once again
    # TODO: make date more flexible? What does date even end up changing?
    if ra and dec:
        log.info("Retrieving the sensitivity from the beam model...")
        frb_master_url = frb_master_base_url + "/v1/code/beam-model/get-sensitivity"
        log.info("frb_master_url: {}".format(frb_master_url))
        payload = {"ra": ra, "dec": dec, "date": t_injection, "beam": beam_no}
        resp = requests.post(frb_master_url, json=payload)
        sensitivities = np.float32(np.array(resp.json()["sensitivities"]))
        if gaussian_modulation is not None:
            log.info(
                "Setting spectral modulation to gaussian profile multiplied by beam model sensitivity..."
            )
            spectral_modulation = gaussian_modulation * sensitivities
        else:
            log.info(
                "Setting spectral modulation to the beam model sensitivity at given RA/Dec..."
            )
            spectral_modulation = sensitivities
    else:
        log.info("RA and/or DEC == 0.0, not applying beam model, continuing...")
        spectral_modulation = gaussian_modulation

    # Get the tcp server address for the given beam_no from frb-master
    # log.info("Retrieving rpc_server address for beam {}...".format(beam_no))
    # frb_master_url = frb_master_base_url + "/v1/parameters/get-beam-info/{}".format(
    #    beam_no
    # )
    # resp = requests.get(frb_master_url)
    # rpc_server = resp.json()["rpc_server"]
    # Use the "heavy" rpc port for data injection
    # if rpc_server.endswith("5555"):
    #    rpc_server = rpc_server.replace("5555", "5556")
    # log.info("rpc_server: {}".format(rpc_server))

    # Make the request to the RpcClient to inject the pulse at fpga0
    # NOTE: offsetting beam_no by 10000 as that is currently how we verify
    # injection beams
    beam_no_offset = beam_no + 10000
    log.info("Injecting into beam {}...".format(beam_no_offset))
    log.debug("Pulse argument: {}".format(pulse))
    resp = client.inject_single_pulse(
        beam_no_offset,
        pulse,
        fpga0,
        wait=True,
        nfreq=nfreq,
        spectral_modulation=spectral_modulation,
    )
    log.info("Inject single pulse response: {}".format(resp))
    log.info("Injection completed!")
    if plot:
        log.info("Request to plot pulse given, plotting...")
        fn = os.path.join(
            plot_dir, "injection_{}".format(t_injection_dt.strftime("%Y%m%d_%H%M%S"))
        )
        make_pulse_plot(pulse, spectral_modulation=spectral_modulation, fn=fn)
        log.info("Plot generated at {}!".format(fn))
    return {"injection": True}


if __name__ == "__main__":
    SIMPULSE_HOST = os.environ.get("SIMPULSE_HOST", "0.0.0.0")
    SIMPULSE_PORT = os.environ.get("SIMPULSE_PORT", "8221")
    app.run(host=SIMPULSE_HOST, debug=True, port=SIMPULSE_PORT)
