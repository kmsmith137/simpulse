import click
import numpy as np
import requests

import simpulse


def gaussian(arr, x_0, f_low, f_hi, fwhm):
    assert len(arr.shape) == 2, "Error: Input array must be 2D!"
    x = np.linspace(f_low, f_hi, arr.shape[0])
    sigma = fwhm/(2*np.sqrt(2*np.log(2)))
    exp = np.exp(-1/2 * ((x-x_0)/sigma)**2)
    # TODO: properly normalizing?
    # a = 1/(np.sum(exp))
    a = 1
    # Do some transposing so multiplication works along frequency axis
    return (arr.T * a * exp).T

def distribute(distributor_base_url, 

@click.command()
@click.option(
    "--nt",
    "nt",
    required=True,
    type=int,
    default=1024,
    help="Number of time samples in the generated pulse array"
)
@click.option(
    "--nfreq",
    "nfreq",
    required=True,
    type=int,
    default=16384,
    help="Number of (equally spaced) frequency samples in the generated pulse array"
)
@click.option(
    "--freq_lo_MHz",
    "freq_lo_MHz",
    required=True,
    type=float,
    default=400,
    help="Frequency value of the bottom channel in the generated array"
)
@click.option(
    "--freq_hi_MHz",
    "freq_hi_MHz",
    required=True,
    type=float,
    default=800,
    help="Frequency value of the bottom channel in the generated array"
)
@click.option(
    "--dm",
    "dm",
    required=True,
    type=float,
    default=0,
    help="DM of the pulse to be generated"
)
@click.option(
    "--sm",
    "sm",
    required=True,
    type=float,
    default=0,
    help="SM (scattering time in ms at 1GHz) of the pulse to be generated"
)
@click.option(
    "--intrinsic_width",
    "width",
    required=True,
    type=float,
    default=0.001,
    help="Width of the pulse to be generated, in seconds"
)
@click.option(
    "--fluence",
    "fluence",
    required=True,
    type=float,
    # TODO: figure out fluence units
    help="Fluence of the pulse to be generated in units * s"
)
@click.option(
    "--spectral_index",
    "spindex",
    required=True,
    type=float,
    default=0,
    help="Spectral index of the pulse to be generated"
)
@click.option(
    "--t_arrival",
    "t_arrival",
    required=True,
    type=int,
    default=0,
    help="Arrival time of the pulse in seconds as freq -> infinity relative to an arbitrary origin"
)
# @click.option(
#     "--model",
#     "model",
#     required=True,
#     type=click.Choice(["gaussian", "powerlaw"]),
#     default="powerlaw"
# )
@click.option(
    "--gaussian_central_freq",
    "gaussian_central_freq",
    required=False,
    type=float,
    default=None,
    help="The frequency at which the gaussian peaks (in MHz) for a gaussian spectrum modulation"
)
# TODO: check if more intuitive to use sigmas with others
@click.option(
    "--gaussian_fwhm",
    "gaussian_fwhm",
    required=False,
    type=float,
    default=None,
    help="The FWHM (in MHz) for a gaussian spectrum modulation"
)

def main(
    nt,
    nfreq,
    freq_lo_MHz,
    freq_hi_MHz,
    dm,
    sm,
    width,
    fluence,
    spindex,
    t_arrival,
    gaussian_central_frequency,
    gaussian_fwhm
):
    # Generate a pulse using simpulse and the specified params
    pulse = simpulse.single_pulse(
        nt,
        nfreq,
        freq_lo_MHz,
        freq_hi_MHz,
        dm,
        sm,
        width,
        fluence,
        spindex,
        t_arrival
    )

    # Figure out what the array size needs to be to contain the dispersed pulse
    # TODO: should sampling time be 0.0009xx?
    dt_sample = 0.001 # in s
    t_end = pulse.get_endpoints()[1]
    n_chunks = int(t_end/(nt*dt_sample) + 1)
    data = np.zeros((nfreq, nt * n_chunks), dtype=np.float32)
    # TODO: should I make n chunks, or just one big intensity array?
    chunk_t0 = 0
    chunk_t1 = n_chunks * nt * dt_sample
    pulse.add_to_timestream(data, chunk_t0, chunk_t1, freq_hi_to_lo=True)
    if gaussian_central_frequency and gaussian_fwhm:
        # Modulate the frequency spectrum bye a gaussian profile
        data = gaussian(
            data, 
            gaussian_central_frequency, 
            freq_lo_MHz, 
            freq_hi_MHz, 
            gaussian_fwhm
        )
    # TODO: make a request to distributor service and put array into a bucket
    #np.save("test_pulse", data)

if __name__ == "__main__":
    main()
