#ifndef _SIMPULSE_PULSARS_HPP
#define _SIMPULSE_PULSARS_HPP

#if (__cplusplus < 201103) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
#error "This source file needs to be compiled with C++0x support (g++ -std=c++0x)"
#endif

#include <vector>
#include <memory>

namespace simpulse {
#if 0
}; // pacify emacs c-mode
#endif


// -------------------------------------------------------------------------------------------------
//
// Pulsar stuff


struct phase_model {
    // Can be called with nderivs=0 to get the phase model Phi(t), or nderivs > 0 to get derivatives of phi
    virtual double eval(double t, int nderivs=0) const = 0;

    // Returns mean phi over a given range.
    virtual double eval_mean_phi(double t1, double t2) const = 0;

    // Vectorized version of eval().
    // No 'nderivs' argument for now, but I may add it later!
    virtual void eval_phi(int nt, double *phi_out, const double *t_in) const = 0;

    static std::shared_ptr<phase_model> make_constant_acceleration(double phi0, double omega0, double omega_dot, double t0=0.0);
};


// For now, the only pulse_profile we consider is the von Mises profile.
// We define the duty cycle to be (pulse FWHM) / (pulse period).
struct von_mises_profile {
    double duty_cycle = 0.0;
    int nphi = 0;
    int nphi2 = 0;   // = nphi/2+1

    std::vector<double> profile_fft;      // length nphi2, normalized to profile_fft[0]=1.
    std::vector<double> profile_antider;  // padded to length (nphi+1)

    // An interpolation used to store the function N(x) = \sum_{n\ne 0} | rho_n j_0(nx/2) |^2 as a function of x.
    std::vector<double> ntab;
    double ntab_xmax;
    int ntab_size;

    // If min_nphi=0, then a reasonable default will be chosen (recommended).
    von_mises_profile(double duty_cycle, int min_nphi=0);

    // Note: 'out' is an array of length nt, but 'phi' is an array of length (nt+1).
    // The i-th sample is assumed to span phase range (phi[i], phi[i+1]).
    // This routine is pretty fast (~5e7 samples/sec on my laptop) but could be vectorized to run faster.
    // The amplitude has units of intensity.
    void eval(int nt, double *out, const double *phi, bool detrend, double amplitude=1.0) const;

    //
    // Returns the amplitude which gives signal-to-noise ratio 1.  Strictly speaking, the returned
    // amplitude is an approximation which is accurate in the limit where the observation contains
    // many pulse periods, and the pulse frequency doesn't change much over the observation.
    //
    //    eta = a noise parameter with units intensity-time^(1/2)
    //    omega = angular pulse frequency
    //    tsamp = timestream sample length
    //    T = total timestream length
    //
    // Note that 'omega' and 'tsamp' are only used for determining the finite-sample correction.
    // This correction can be disabled by either setting omega=0 or tsamp=0.
    //
    // FIXME: it may be convenient later to implement a vectorized version which evaluates
    // on a grid of omega values.
    //
    double get_normalization(double eta, double omega, double tsamp, double T, bool detrend) const;

    static int default_nphi(double duty_cycle);
};


}  // namespace simpulse

#endif // _SIMPULSE_PULSARS_HPP
