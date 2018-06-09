#ifndef _SIMPULSE_SINGLE_PULSE_HPP
#define _SIMPULSE_SINGLE_PULSE_HPP

#if (__cplusplus < 201103) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
#error "This source file needs to be compiled with C++0x support (g++ -std=c++0x)"
#endif

#include <iostream>
#include <string>
#include <vector>
#include <memory>

namespace simpulse {
#if 0
}; // pacify emacs c-mode
#endif


// -------------------------------------------------------------------------------------------------
//
// class single_pulse: represents one dispersed, scattered pulse, in a fixed frequency channelization.


class single_pulse {
public:
    // These parameters of the pulse cannot be changed after construction.

    const int    pulse_nt;          // number of samples used "under the hood" to represent pulse (1024 is usually a good default)
    const int    nfreq;             // number of frequency channels, assumed equally spaced
    const double freq_lo_MHz;       // lower limit of frequency band
    const double freq_hi_MHz;       // upper limit of frequency band
    const double dm;                // dispersion measure in its standard units (pc cm^{-3})
    const double sm;                // we define the "scattering measure" SM to be scattering time in milliseconds at 1 GHz
    const double intrinsic_width;   // frequency-independent Gaussian width in seconds

    // These parameters can be changed after construction.

    double fluence;                   // F(nu_0) = int I(t) dt evaluated at central frequency
    double spectral_index;            // F(nu) = F(nu_0) (nu/nu_0)^alpha
    double undispersed_arrival_time;  // arrival time at high frequency, in seconds, and relative to the same origin used in add_to_timestream()
    
    //
    // Constructor arguments:
    //
    //   pulse_nt = number of samples used "under the hood" to represent the pulse (suggest power of two; 1024 is usually a good choice).
    //   nfreq, freq_lo_MHz, freq_hi_MHz = specification of frequency channels (assumed equally spaced).
    //   dm = dispersion measure in its standard units (pc cm^{-3}).
    //   sm = scattering measure, which we define to be scattering time in milliseconds (not seconds!) at 1 GHz.
    //   intrinsic_width = frequency-independent Gaussian width in seconds (not milliseconds).
    //   fluence = integrated flux in units Jy-s (where "Jy" really means "whatever units the output of add_to_timestream() has")
    //   undispersed_arrival_time = arrival time of pulse as freq->infty, in seconds relative to an arbitrary origin
    //
    single_pulse(int pulse_nt, int nfreq, double freq_lo_MHz, double freq_hi_MHz, 
		 double dm, double sm, double intrinsic_width, double fluence, 
		 double spectral_index, double undispersed_arrival_time);

    void set_fluence(double fluence);
    void set_spectral_index(double spectral_index);
    void set_undispersed_arrival_time(double undispersed_arrival_time);

    // Returns the earliest and latest arrival times in the band [freq_lo_MHz, freq_hi_MHz].
    // (Note that both of these will be larger than single_pulse::undispersed_arrival_time, unless the intrinsic width is very large).
    void get_endpoints(double &t0, double &t1) const;

    // This routine adds the pulse to a "block" of (frequency, time) samples.
    // It is sometimes called incrementally, as a stream of blocks generated.
    //
    // The 'out_t0' and 'out_t1' args are the endpoints of the sampled region, in seconds relative 
    // to the undispersed_arrival_time.
    //
    // The 'out' arg is a 2d array with shape (nfreq, out_nt).
    // The 'stride' is the pointer offset between frequencies.  If stride=0 then it defaults to out_nt.
    //
    // The frequencies are assumed ordered from lowest to highest.
    // WARNING: this is the opposite of the ordering used in rf_pipelines and bonsai!!
    //
    // If the highest-to-lowest ordering is desired, it can be obtained by using a negative stride
    // (and taking 'out' to point to the last row of the array rather than the first row)

    template<typename T> 
    void add_to_timestream(T *out, double out_t0, double out_t1, int out_nt, int stride=0, double weight=1.) const;
    
    //
    // Returns total signal-to-noise for all frequency channels and time samples combined.
    // The signal-to-noise of a sampled pulse depends on 'sample_dt', the length of a sample in seconds.
    //
    // In principle, it also depends on 'sample_t0', the starting time of an arbitrarily chosen sample,
    // although this dependence will be weak in realistic cases!
    //
    double get_signal_to_noise(double sample_dt, double sample_rms=1.0, double sample_t0=0.0) const;
    
    // Here is a version of get_signal_to_noise() which takes length-nfreq arrays 'sample_rms' and 'channel_weights'.
    // If 'channel_weights' are unspecified, then 1/sample_rms^2 weighting will be assumed.
    double get_signal_to_noise(double sample_dt, const double *sample_rms, const double *channel_weights=NULL, double sample_t0=0.0) const;

    // String representation
    void print(std::ostream &os) const;
    std::string str() const;

    // We make the single_pulse noncopyable, even though the default copy constructor is a sensible "deep" copy.  
    // We do this to catch performance bugs, since a deep copy of a single_pulse is probably unintentional.
    single_pulse(const single_pulse &) = delete;
    single_pulse& operator=(const single_pulse &) = delete;

protected:

    // Under-the-hood representation of the pulse.
    std::vector<double> pulse_t0;        // start time of pulse in i-th channel, relative to undispersed_arrival_time
    std::vector<double> pulse_t1;        // end time of pulse in i-th channel, relative to undispersed_arrival_time
    std::vector<double> pulse_freq_wt;   // weight of pulse in i-th channel, determined by spectral index.
    std::vector<double> pulse_cumsum;    // shape (nfreq, pulse_nt+1) array, containing cumulative sum of pulse, normalized to sum=1
    double min_t0;                       // minimum of all pulse_t0 values
    double max_t1;                       // maximum of all pulse_t1 values
    double max_dt;                       // maximum of all (pulse_t1-pulse_t0) values

    // Internal helper functions follow.

    void _compute_freq_wt();

    template<typename T>
    void _add_pulse_to_frequency_channel(T *out, double out_t0, double out_t1, int out_nt, int ifreq, double weight) const;    
};


inline std::ostream &operator<<(std::ostream &os, const single_pulse &sp) { sp.print(os); return os; }


}  // namespace simpulse

#endif // _SIMPULSE_SINGLE_PULSE_HPP
