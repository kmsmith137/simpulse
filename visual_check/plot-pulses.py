#!/usr/bin/env python
"""
This is a throwaway script from when I was visually debugging the pulse simulation code.
Some day, I may replace it with a more rigorous test.  The 'reference_*.png' files show 
what the results should look like.

to run, simply do $./plot-pulses.py
"""
import numpy as np
import matplotlib; matplotlib.use('Agg') # 'Agg' enables silent output to file
import matplotlib.pyplot as plt
import simpulse # contains objects and modules for the simulation

# this is the main function for plotting a simulated pulse object, which is obtained by
# sp = simpulse.single_pulse(pulse_nt, nfreq, freq_lo_MHz, freq_hi_MHz, dm, sm, 
#      intrinsic_width, fluence, spectral_index, undispersed_arrival_time)
def make_plot(sp, ifreq_list, color_list, label_list, filename):

    # check whether inputs are valid objects
    assert isinstance(sp, simpulse.single_pulse)
    assert len(ifreq_list) == len(color_list) == len(label_list)

    # given sp, get the end points of the time axis
    (plot_t0, plot_t1) = sp.get_endpoints()

    # loop over various plot types, e.g., with(out) labels
    for (plot_nt, ls, lflag) in [ (100,'-',True), (1000,':',False) ]:
        
        # make an empty array of (freq, time)
        ts = np.zeros((sp.nfreq, plot_nt))
        
        # add the pulse to the array, constrained by the end points
        sp.add_to_timestream(ts, plot_t0, plot_t1)

        # initialize empty vectors for holding plot data
        # upsample by 2 for a histogram look
        x = np.zeros(2*plot_nt)
        y = np.zeros(2*plot_nt)
        
        # assign time values
        t = np.linspace(plot_t0, plot_t1, plot_nt+1)
        x[0::2] = t[:-1]
        x[1::2] = t[1:]
        
        # looping over a set of freq channels, assign freq values to the y-axis vector
        for (ifreq, color, label) in zip(ifreq_list, color_list, label_list):
            
            #ISSUE: making the code more general, maybe ifreq need not to be an integer?
            y[0::2] = ts[ifreq,:]
            y[1::2] = ts[ifreq,:]
        
            # choose a good font
            font = {'family':'serif', 'size': 14}
            matplotlib.rc('font', **font)
            
            # plot with(out) label
            if lflag:
                plt.plot(x, y, color=color, ls=ls, label=label)
            else:
                plt.plot(x, y, color=color, ls=ls)

    plt.legend(loc='upper right').draw_frame(False)
    plt.savefig(filename)
    plt.clf()
    print 'wrote', filename

def plot1():
    """
    Some Gaussians with a little bit of dispersion and scattering.

    Visual checks:
       - arrival times of pulses should be 141.5, 118.4, 110.4 ms
       - fluences should be { 13.3, 30, 53.3 }
       - low frequency pulse should show scattering
       - no dispersion broadening visible
    """
    
    pulse_nt = 1024
    nfreq = 512
    freq_lo_MHz = 1000.0
    freq_hi_MHz = 2000.0
    dm = 10.0
    sm = 4.0
    intrinsic_width = 0.005
    fluence = 30.0
    spectral_index = 2.0
    undispersed_arrival_time = 0.100
    
    sp = simpulse.single_pulse(pulse_nt, nfreq, freq_lo_MHz, freq_hi_MHz, dm, sm, intrinsic_width, fluence, spectral_index, undispersed_arrival_time)

    make_plot(sp, 
              ifreq_list = [0, 256, 511], 
              color_list = ['r', 'b', 'm'],
              label_list = [ '1 GHz', '1.5 GHz', '2.0 GHz' ],
              filename = 'plot1.png')

def plot2():
    """Boxcars labeled by time intervals"""
    
    pulse_nt = 1024
    nfreq = 7
    freq_lo_MHz = 1000.0
    freq_hi_MHz = 2000.0
    dm = 10.0
    sm = 0.0
    intrinsic_width = 0.0
    fluence = 1.0
    spectral_index = 0.0
    undispersed_arrival_time = 0.100

    sp = simpulse.single_pulse(pulse_nt, nfreq, freq_lo_MHz, freq_hi_MHz, dm, sm, intrinsic_width, fluence, spectral_index, undispersed_arrival_time)

    make_plot(sp,
              ifreq_list = [0, 2, 4, 6], 
              color_list = ['r', 'b', 'm', 'k'],
              label_list = [ '131.8-141.5', '116.8-120.3', '114.1-116.8', '110.4-112.0' ],
              filename = 'plot2.png')

def plot3():
    """
    Same as plot1, but with intrisic_width set to zero, so that pulse width is dominated by scattering,
    and spectral_inex set to zero.
    """
    
    pulse_nt = 1024
    nfreq = 512
    freq_lo_MHz = 1000.0
    freq_hi_MHz = 2000.0
    dm = 10.0
    sm = 4.0
    intrinsic_width = 0.0
    fluence = 30.0
    spectral_index = 0.0
    undispersed_arrival_time = 0.100
    
    sp = simpulse.single_pulse(pulse_nt, nfreq, freq_lo_MHz, freq_hi_MHz, dm, sm, intrinsic_width, fluence, spectral_index, undispersed_arrival_time)

    make_plot(sp, 
              ifreq_list = [0, 256, 511], 
              color_list = ['r', 'b', 'm'],
              label_list = [ '1 GHz', '1.5 GHz', '2.0 GHz' ],
              filename = 'plot3.png')

plot1()
plot2()
plot3()
