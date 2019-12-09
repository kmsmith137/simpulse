import click
from datetime import datetime
import json
import logging
import numpy as np
import os
import requests
import time
from flask import Flask, request

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
    sm = data["sm"]
    width = data["width"]
    fluence = data["fluence"]
    spindex = data["spindex"]
    t_injection = data["t_injection"]
    gaussian_central_freq = data["gaussian_central_freq"]
    gaussian_fwhm = data["gaussian_fwhm"]
    frb_master_base_url = data["frb_master_base_url"]
    fpga0 = data["fpga0"]
    
    # Make some basic assertions
    assert 0 < beam_no < 255 or 1000 < beam_no < 1255 or 2000 < beam_no < 2255 or 3000 < beam_no < 3255, "ERROR: invalid beam_no!"
    t_injection_dt = datetime.strptime(t_injection, "%Y-%m-%dT%H:%M:%S.%fZ")
    assert t_injection_dt > datetime.utcnow(), "ERROR: t_injection must not be in the past!"

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
    #log.info("Retrieving rpc_server address for beam {}...".format(beam_no))
    #frb_master_url = frb_master_base_url + "/v1/parameters/get-beam-info/{}".format(
    #    beam_no
    #)
    #resp = requests.get(frb_master_url)
    #rpc_server = resp.json()["rpc_server"]
    # Use the "heavy" rpc port for data injection
    #if rpc_server.endswith("5555"):
    #    rpc_server = rpc_server.replace("5555", "5556")
    #log.info("rpc_server: {}".format(rpc_server))

    # Create an RpcClient object which will send a pulse injection request to L1
    # TODO: better server name (do these need unique names?)
    #servers = {"injections-test": rpc_server}
    #client = RpcClient(servers=servers)

    # Make the request to the RpcClient to inject the pulse at fpga0
    # NOTE: offsetting beam_no by 10000 as that is currently how we verify
    # injection beams
    # TODO: modify RpcClient to take in a frequency modulation array when
    # injecting a pulse
    beam_no_offset = beam_no + 10000
    log.info("Injecting into beam {}...".format(beam_no_offset))
    resp = client.inject_single_pulse(beam_no_offset, pulse, fpga0, wait=True, nfreq=nfreq)
    log.info("Injection completed!")
    return {"injection": True}


if __name__ == "__main__":
    app.run(host="0.0.0.0", debug=True, port=8004)
