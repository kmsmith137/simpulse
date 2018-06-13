#ifndef _SIMPULSE_PULSAR_PROFILES_HPP
#define _SIMPULSE_PULSAR_PROFILES_HPP

#if (__cplusplus < 201103) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
#error "This source file needs to be compiled with C++0x support (g++ -std=c++0x)"
#endif

#include <vector>
#include <memory>


namespace simpulse {
#if 0
}; // pacify emacs c-mode
#endif


// Defined in pulsar_phase_models.hpp
struct phase_model_base;


// -------------------------------------------------------------------------------------------------
//
// von_mises_profile
//
// This class represents a von Mises profile, which can be used to simulate pulsars.
// (Currently, the von Mises profile is the only periodic pulse profile implemented in simpulse,
// although in the future we might implement more profiles and define a 'pulsar_profile' base class.)
//
// In order to simulate a pulsar, you need two things: a phase model and a pulse profile.  Then, to do
// the simulation, you can call either profile.eval_integrated_samples() or profile.add_integrated_samples().
// These methods take the phase model as a parameter.
//
// By default, the profile is normalized so that its peak flux is equal to 1 (before applying the detrending subtraction).
// To get/set the amplitude, see below (methods get_peak_flux(), ..., set_multi_pulse_signal_to_noise()).
//
// Mathematically, a profile is just a function rho(Phi) which gives the flux 'rho' as a function of pulse phase 'phi'.
// The von Mises profile is the functional form:
//
//       rho(phi) = exp[ -2 kappa sin(pi*phi)^2 ]
//
// where kappa is a narrowness parameter, related to the duty cycle D by kappa = log(2) / (2 sin^2(pi*D/2)).


class von_mises_profile {
public:
    
    // The 'duty_cycle' is defined as D = (pulse full width at half maximum) / (pulse period).
    // A reasonable choice is D = 0.1 or so.
    //
    // If the boolean variable 'detrend' is true, then the mean will be subtracted from the profile.
    //
    // It is unlikely that you'll need to set the 'min_internal_nphi' constructor argument, which
    // changes the number of phase bins used internally to represent the pulse.  If set to zero, then
    // a reasonable default value will be chosen.
    //
    // It is unlikely that you'll need to set the 'internal_phi_block_size' argument, which determines
    // the block size for calls to phase_model_base::eval_phi_sequence().  If set to zero, then a
    // reasonable default value will be chosen.

    von_mises_profile(double duty_cycle, bool detrend, int min_internal_nphi=0, int internal_phi_block_size=0);

    const double duty_cycle;
    const bool detrend;
    const int internal_nphi;
    const int internal_phi_block_size;
    const double kappa;      // von Mises profile is exp(-2 kappa sin(pi*phi)^2)

    // These are the main routines used to simulate a pulsar in a regularly spaced sequence of time samples.
    // The 'out' argument should point to an array of length nt.
    //
    // The 't0' argument should be the _beginning_ of the first time sample, and 't1' should be the _end_
    // of the last time sample.  Thus, t1=t0+nt*dt, where dt is the length of a sample (not t1=t0+(nt-1)*dt).
    //
    // Reminder: if the 'detrend' flag was specified at construction, then the simulated flux will be detrended
    //
    // FIXME: add_integrated_samples() is currently not unit-tested (or python-wrapped).

    template<typename T> void eval_integrated_samples(T *out, double t0, double t1, ssize_t nt, const phase_model_base &pm) const;
    template<typename T> void add_integrated_samples(T *out, double t0, double t1, ssize_t nt, const phase_model_base &pm) const;


    // The overall amplitude of the pulses can be specified by either the 'peak_flux' or the 'mean_flux'.
    // Note: these refer to the peak/mean flux without any detrending offset subtracted! (even if the 'detrend' flag is set)

    double get_peak_flux() const;
    double get_mean_flux() const;
    void set_peak_flux(double peak_flux);
    void set_mean_flux(double mean_flux);

    // The "multi-pulse signal-to-noise" is the total SNR for a train of pulses with duration 'total_time'.
    //
    // The SNR calculations account for finite time resolution (and detrending, if detrend=True was specified
    // at construction).  Strictly speaking, the calculated SNR is an approximation to the true SNR, which may 
    // slightly depend on the exact arrival times of the pulses.
    //
    // In functions which calculate SNR:
    //   - The 'total_time' argument is the total duration of the pulse train.
    //   - The 'dt_sample' argument is the length of each time sample.
    //   - The 'pulse_freq' argument is the pulse frequency.
    //   - The 'sample_rms' argument is the RMS noise fluctuation in each time sample.

    double get_single_pulse_signal_to_noise(double dt_sample, double pulse_freq, double sample_rms=1.0) const;
    double get_multi_pulse_signal_to_noise(double total_time, double dt_sample, double pulse_freq, double sample_rms=1.0) const;
    void set_single_pulse_signal_to_noise(double snr, double dt_sample, double pulse_freq, double sample_rms=1.0);
    void set_multi_pulse_signal_to_noise(double snr, double total_time, double dt_sample, double pulse_freq, double sample_rms=1.0);
    

    // Returns the instantaneous flux evaluated at pulse phase 'phi'.
    // Reminder: if the 'detrend' flag was specified at construction, then the simulated flux will be detrended.
    double eval_instantaneous(double phi) const;    

    
    // Returns the Fourier transform of the profile
    //     rho_m = int_0^1 dphi rho(phi) e^{i m phi}.
    //
    // Note that rho_m is real, and rho_m = rho_{-m}, since the von Mises profile is symmetric.
    // The DC mode rho_0 will equal 'mean_flux' if detrend=False, or 0 if detrend=True.
    //
    // The return value is a 1D array of length 'nout'.  If nout=0 (the default), then it defaults to
    // (internal_nphi/2+1), the number of Fourier coefficients which are computed internally.  (If 'nout'
    // is larger than this, then the returned array is zero-padded.
    //
    // The returned FFT does not include a time-sample window function.  You may want to multiply the
    // output by j_0(pi*dphi/2), where dphi = (pulse_freq * dt_sample).

    template<typename T> void get_profile_fft(T *out, int nout) const;


    // This method is intended for debugging (hence the "_slow"!)
    // Returns the average flux over phase (not time) interval [phi0, phi1].
    double eval_integrated_sample_slow(double phi0, double phi1) const;

protected:
    const int internal_nphi2;

    double _peak_flux = 1.0;      // used internally to represent the amplitude
    double _mf_multiplier = 0.0;  // ratio (mean_flux / peak_flux), constant after construction
    
    // Note: padded to length (internal_nphi+1), for internal convenience interpolating.
    // These are always detrended (even if the 'detrend' flag is not set)
    std::vector<double> detrended_profile;
    std::vector<double> detrended_profile_antider;
    
    // Length internal_nphi2, normalized to peak_flux=1.
    // This is detrended if the 'detrend' flag is set.
    std::vector<double> profile_fft;

    mutable std::vector<double> phi_tmp;   // length (internal_phi_block_size + 1)

    double _get_rho2(double dphi) const;
};


}  // namespace simpulse

#endif // _SIMPULSE_PULSAR_PROFILES_HPP
