# Set non-interactice plotting backend
import matplotlib

matplotlib.use("Agg")

from datetime import datetime
import logging
import os
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

# Initialize the Flask application
app = Flask(__name__)

# TODO: testing with a single beam, so forcing the rpc client to the beam that's
# running duplication code now for testing
servers = {"injections": "tcp://10.8.210.16:5556"}
client = RpcClient(servers=servers)

# Set the location to save plots
plot_dir = "/intensity-data/processed/frb-injections"


####################
# HELPER FUNCTIONS #
####################


def make_pulse_plot(pulse, spectral_model=None, beam_model=None, fn=None):
    """
    A function to create a plot of the pulse to be injected. The plot
    is created by generating a dedispersed version of the pulse, and
    injecting it into a nfreq by nt array of zeroes.

    INPUT:
    ------

    pulse : simpulse object
        The simpulse object to be injected
    spectral_model : numpy array, optional
        A numpy array to modulate the spectrum of the pulse. None by default
    spectral_model : numpy array, optional
        A numpy array to apply the beam model to the pulse. None by default
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
    max_t = max_nt / 1024.0 * 1000.0 - 500.0
    min_t = min_nt / 1024.0 * 1000.0 - 500.0
    data_cut = data[:, min_nt:max_nt]
    if spectral_model is not None:
        spectral_model = np.array(spectral_model)
        data_cut = data_cut * spectral_model[:, np.newaxis]
    if beam_model is not None:
        beam_model = np.array(beam_model)
        data_cut = data_cut * beam_model[:, np.newaxis]

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
    log.info("Unpacking data from mimic request...")
    data = request.get_json()
    nt_per_sample = data["nt_per_sample"]
    nfreq = data["nfreq"]
    freq_lo_MHz = data["freq_lo_MHz"]
    freq_hi_MHz = data["freq_hi_MHz"]
    beam_no = data["beam_no"]
    dm = data["dm"]
    tau_1_ghz_ms = data["tau_1_ghz_ms"]
    pulse_width_ms = data["pulse_width_ms"]
    # Simpulse expects width in s
    pulse_width_s = pulse_width_ms / 1000.0
    # The input width to the injection system is the FWHM, but simpulse takes
    # the Gaussian sigma as the width, so divide by the appropraite factor
    pulse_width_s = pulse_width_s / (2. * np.sqrt(2 * np.log(2)))
    fluence_jy_ms = data["fluence_jy_ms"]
    # Simpulse expects fluence in (units) s
    fluence_jy_s = fluence_jy_ms / 1000.0
    t_injection = data["injection_time"]
    fpga0 = data["timestamp_fpga_injection"]
    spectral_model = data["spectral_model"]
    beam_model = data["beam_model"]
    plot = data["plot"]
    log.info("Data successfully unpacked!\n")

    # Make some basic assertions
    log.info("Asserting beam number is valid...")
    assert (
        0 < beam_no < 256
        or 1000 < beam_no < 1256
        or 2000 < beam_no < 2256
        or 3000 < beam_no < 3256
    ), "ERROR: invalid beam_no!"
    log.info("Beam number is valid!\n")
    log.info("Asserting injection time is not in the past...")
    t_injection_dt = datetime.strptime(t_injection, "%Y-%m-%dT%H:%M:%S.%fZ")
    assert (
        t_injection_dt > datetime.utcnow()
    ), "ERROR: t_injection must not be in the past!"
    log.info("Injection time is not in the past!")

    # Generate a pulse using simpulse and the specified params
    # (Note: the value set to 0.0 is for spectral index. I set it to 0 because
    # the spectral modulation folds in the spectral index already)
    t_inf = 0
    log.info("Generating pulse with simpulse...")
    pulse = simpulse.single_pulse(
        nt_per_sample,
        nfreq,
        freq_lo_MHz,
        freq_hi_MHz,
        dm,
        tau_1_ghz_ms,
        pulse_width_s,
        fluence_jy_s,
        0.0,
        t_inf,
    )
    log.info("Successfully generated pulse!")

    # Make the request to the RpcClient to inject the pulse at fpga0
    # NOTE: offsetting beam_no by 10000 as that is currently how we verify
    # injection beams
    beam_no_offset = beam_no + 10000
    log.info("Injecting into beam {}...".format(beam_no_offset))
    log.info("Mean of beam_model: {}".format(np.mean(beam_model)))
    log.info("Sum of spectral model: {}".format(np.sum(spectral_model)))
    log.info("Simpulse pulse: {}".format(pulse))
    resp = client.inject_single_pulse(
        beam_no_offset,
        pulse,
        fpga0,
        wait=True,
        nfreq=nfreq,
        spectral_model=spectral_model,
        beam_model=beam_model,
    )
    log.info("Inject single pulse response: {}".format(resp))
    log.info("Injection completed!")
    if plot:
        log.info("Request to plot pulse given, plotting...")
        fn = os.path.join(
            plot_dir, "injection_{}".format(t_injection_dt.strftime("%Y%m%d_%H%M%S"))
        )
        make_pulse_plot(
            pulse, spectral_model=spectral_model, beam_model=beam_model, fn=fn,
        )
        log.info("Plot generated at {}.png!".format(fn))
    return {"injection": True}


if __name__ == "__main__":
    SIMPULSE_HOST = os.environ.get("SIMPULSE_HOST", "0.0.0.0")
    SIMPULSE_PORT = os.environ.get("SIMPULSE_PORT", "8221")
    app.run(host=SIMPULSE_HOST, debug=True, port=SIMPULSE_PORT)
