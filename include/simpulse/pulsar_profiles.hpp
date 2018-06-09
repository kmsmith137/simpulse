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
// Pulsar profiles.
//
// For now, the only pulse profile we implement is the von Mises profile.
// Eventually, we might implement more profiles and define a 'pulsar_profile' base class.
//
// The von Mises profile is
//
//   f(phi) = exp[ -2 kappa sin(pi*phi)^2 ]
//
// where 'phi' is the pulse phase (defined so that peak intensity occurs at integer values of phi),
// and 'kappa' is a width parameter (larger kappa corresponds to a narrower pulse).  The parameter
// 'kappa' is related to the duty cycle D (defined as the pulse FWHM divided by the pulse period) by:
//
//   kappa = log(2) / (2 sin^2(pi*D/2))


class von_mises_profile {
public:
    // If the 'detrend' flag is set, an offset will be subtracted from the pulse profile, so that the mean is zero.
    // If 'min_internal_nphi' is zero, then a reasonable default value of 'internal_nphi' will be chosen (recommended).
    von_mises_profile(double duty_cycle, bool detrend, int min_internal_nphi=0);

    const double duty_cycle;
    const bool detrend;
    const int internal_nphi;

    // These are the main routines used to simulate a pulsar.
    //
    // In eval_integrated_samples(), 't0' should be the _beginning_ of the first sample, 
    // and 't1' should be the _end_ of the last sample.  Thus, t1=t0+nt*dt, where dt is 
    // the length of a sample (not t1=t0+(nt-1)*dt).

    void eval_integrated_samples(double *out, double t0, double t1, ssize_t nt, const phase_model_base &pm, double amplitude=1.0) const;

    // void add_integrated_samples(double *out, double t0, double t1, ssize_t nt, const phase_model_base &pm, double amplitude=1.0) const;
    
    // By default, the profile is normalized so that its mean (not peak) value is 1 (before detrending!)
    double point_eval(double phi, double amplitude=1.0) const;    

    // Returns the total SNR of a single pulse, as computed by integrate_samples() with amplitude=1.
    double get_single_pulse_signal_to_noise(double dt_sample, double pulse_freq, double sample_rms=1.0) const;
    double get_multi_pulse_signal_to_noise(double total_time, double dt_sample, double pulse_freq, double sample_rms=1.0) const;

    double get_mean_flux() const { return mean_flux; }
    
    // get_profile_fft(): outputs values of rho_m = int_0^1 dphi rho(phi) e^{2 pi i m phi}.
    // Since we normalize the pulse to mean 1, rho_0 will either be 1 or 0, depending on whether the 'detrend' flag is set.
    // Note that rho_m is real, and rho_m = rho_{-m}, since the von Mises profile is symmetric under phi -> (1-phi).
    // Note: we only compute rho_m for m < internal_nphi2.  If nout > internal_nphi2, then 'out' will be zero-padded.
    // Note: out[0] = xxx

    template<typename T> void get_profile_fft(T *out, int nout) const;

    // For debugging.
    double eval_integrated_sample_slow(double phi0, double phi1, double amplitude=1.0) const;

protected:
    const int internal_nphi2;
    const double kappa;
    double mean_flux = 0.0;
    
    // Note: padded to length (internal_nphi+1), for internal convenience interpolating.
    std::vector<double> detrended_profile;
    std::vector<double> detrended_profile_antider;
    
    // Length internal_nphi2, normalized to profile_fft[0]=1.
    std::vector<double> profile_fft;

    mutable std::vector<double> phi_tmp;   // length (tmp_block_size + 1)
    const ssize_t phi_block_size = 1024;

    double _get_rho2(double dphi) const;
};


}  // namespace simpulse

#endif // _SIMPULSE_PULSAR_PROFILES_HPP
