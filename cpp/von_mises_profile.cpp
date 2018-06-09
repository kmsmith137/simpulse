#include "../include/simpulse.hpp"
#include "../include/simpulse/internals.hpp"

#define DEBUG 1

using namespace std;


namespace simpulse {
#if 0
};  // pacify emacs c-mode!
#endif


// Von Mises profile with no detrending, normalized to peak (not mean!) value 1.
inline double _vm_profile(double kappa, double phi)
{
    return exp(-2. * kappa * square(sin(M_PI * phi)));
}

// Returns the linear interpolation f(u), where f(0)=f0 and f(1)=f1.
inline double _linterp(double u, double f0, double f1)
{
    return (1-u)*f0 + u*f1;
}


// -------------------------------------------------------------------------------------------------
//
// Constructor and helper functions.


static int choose_internal_nphi(double duty_cycle, int min_internal_nphi=0)
{
    if (duty_cycle <= 0.0)
	throw runtime_error("simpulse::von_mises_profile: duty_cycle must be >= 0");
    if (duty_cycle >= 0.5)
	throw runtime_error("simpulse::von_mises_profile: duty_cycle must be < 0.5");

    // No fundamental reason for these restrictions, but violating them is probably unintentional!
    // Note that in the limit duty_cyle->0, this implementation is inefficient, and it would be better
    // to supply an alternate implementation similar to single_pulse.

    if (duty_cycle < 1.0e-4)
	throw runtime_error("simpulse::von_mises_profile: we currently don't support duty cycles < 1.0e-4");
    if (min_internal_nphi < 0)
	throw runtime_error("simpulse::von_mises_profile: 'min_internal_nphi' argument must be >= 0");
    if (min_internal_nphi > 65536)
	throw runtime_error("simpulse::von_mises_profile: we currently don't support nphi >= 65536");    

    // Ultra-conservative value
    // FIXME I'd like to understand better how to choose this
    int default_nphi = round_up(30./duty_cycle);

    return max(default_nphi, min_internal_nphi);
}


von_mises_profile::von_mises_profile(double duty_cycle_, bool detrend_, int min_internal_nphi)
    : duty_cycle(duty_cycle_),
      detrend(detrend_),
      internal_nphi(choose_internal_nphi(duty_cycle_, min_internal_nphi)),   // note: constructor args sanity-checked here
      internal_nphi2(internal_nphi/2 + 1),
      kappa(log(2.0) / (2 * square(sin(M_PI*duty_cycle/2.)))),
      mean_flux(0.0)   // to be initialized below
{
    this->detrended_profile.resize(internal_nphi+1, 0.0);
    this->detrended_profile_antider.resize(internal_nphi+1, 0.0);
    this->profile_fft.resize(internal_nphi2, 0.0);

    double *rho = &detrended_profile[0];
    double *rho_a = &detrended_profile_antider[0];

    for (int iphi = 0; iphi < internal_nphi; iphi++) {
	rho[iphi] = _vm_profile(kappa, iphi / double(internal_nphi));
	mean_flux += rho[iphi] / internal_nphi;
    }

    // Compute 'profile_fft'.

    vector<complex<double>> ctmp(internal_nphi2, 0.0);

    fftw_plan plan = fftw_plan_dft_r2c_1d(internal_nphi, rho, reinterpret_cast<fftw_complex *> (&ctmp[0]), FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    for (int iphi = 0; iphi < internal_nphi2; iphi++)
	profile_fft[iphi] = ctmp[iphi].real() / internal_nphi;

    simpulse_assert(fabs(profile_fft[0] - mean_flux) < 1.0e-13);

    if (detrend)
	profile_fft[0] = 0.0;

    // Complete calculation of 'detrended_profile' (by subtracting mean flux and padding by one), and 'detrended_profile_antider'.

    for (int iphi = 0; iphi < internal_nphi; iphi++)
	rho[iphi] -= mean_flux;

    rho[internal_nphi] = rho[0];  // pad by one (with periodicity)

    for (int iphi = 0; iphi < internal_nphi; iphi++) {
	double x = 0.5 * (rho[iphi] + rho[iphi+1]);
	rho_a[iphi+1] = rho_a[iphi] + x;
    }

    simpulse_assert(fabs(rho_a[internal_nphi] - 0.0) < 1.0e-12);
}


// -------------------------------------------------------------------------------------------------
//
// eval_integrated_samples() and friends


// Inline helper class for periodic interpolation.  We write:
//
//   p = phi*n
//     = j*n + i + u
//
// where i,j are integers, 0 <= i < nphi, and 0 <= u <= 1 (within roundoff error)

struct ihelper {
    double p;
    double u;
    int i;

    double rho0;
    double rho1;
    
    ihelper(double phi, int n, const double *rho)
    {
	ssize_t j = (phi >= 0.0) ? ssize_t(phi) : (-ssize_t(-phi)-1);
	double x = (phi-j) * n;   // in interval [0,n], within roundoff error

	p = phi * n;
	i = int(x);
	i = max(i, 0);     // account for possibility that (i == -1) due to roundoff error
	i = min(i, n-1);   // account for possibility that (i == n) due to roundoff error
	u = x - i;

#if DEBUG
	simpulse_assert((i >= 0) && (i < n));
	simpulse_assert((u >= -1.0e-12) && (u <= 1+1.0e-12));
	simpulse_assert(fabs(p - j*n - i - u) < (1.0e-13 * (1 + fabs(p))));
#endif

	rho0 = rho[i];
	rho1 = rho[i+1];
    }

    inline double _linterp(double uu)
    {
	return (1-uu)*rho0 + uu*rho1;
    }
};


void von_mises_profile::eval_integrated_samples(double *out, double t0, double t1, ssize_t nt, const phase_model_base &pm, double amplitude) const
{
    if (_unlikely(nt <= 0))
	throw runtime_error("von_mises_profile::eval_integrated_samples(): expected nt > 0");
    if (_unlikely(t0 >= t1))
	throw runtime_error("von_mises_profile::eval_integrated_samples(): expected t0 < t1");
    if (_unlikely(!out))
	throw runtime_error("von_mises_profile::eval_integrated_samples(): 'out' is a null pointer");

    const double *rho = &detrended_profile[0];
    const double *rho_a = &detrended_profile_antider[0];
    const double mf = detrend ? 0.0 : mean_flux;

    if (!phi_tmp.size())
	phi_tmp.resize(phi_block_size+1, 0.0);   // note "+1" here

    // Number of time samples simulated so far.
    ssize_t nt0 = 0;

    while (nt0 < nt) {
	// Simulate one "block" of time samples i = nt0, ..., (nt1-1).
	ssize_t nt1 = min(nt, nt0 + phi_block_size);
	double t0_block = _linterp(nt0/double(nt), t0, t1);
	double t1_block = _linterp(nt1/double(nt), t0, t1);

	// Evaluate phase model in this block.
	// Note that the number of 'phi' values is (nt1-nt0+1), not (nt1-nt0).
	pm.eval_phi_sequence(t0_block, t1_block, nt1-nt0+1, &phi_tmp[0]);
	
	ihelper h0(phi_tmp[0], internal_nphi, rho);
	
	for (ssize_t it = nt0; it < nt1; it++) {
	    // Loop invariant: at top, 'h0' corresponds to phase phi_tmp[it].
	    ihelper h1(phi_tmp[it+1], internal_nphi, rho);
	    double dp = h1.p - h0.p;
	    double integral;

	    if (_unlikely(dp <= 0))
		throw runtime_error("von_mises_profile::eval_integrated_samples(): phase model is not monotone increasing?!");

	    // I put way too much thought into writing the logic below in a way which minimizes roundoff error!

	    if (h0.i == h1.i)
		integral = (h1.u - h0.u) * h0._linterp(0.5 * (h0.u + h1.u));
	    else {
		integral = (1.0 - h0.u) * h0._linterp(0.5*h0.u + 0.5);
		integral += rho_a[h1.i] - rho_a[h0.i+1];
		integral += h1.u * h1._linterp(0.5*h1.u);
	    }

	    out[it] = amplitude * (integral/dp + mf);
	    h0 = h1;
	}

	nt0 = nt1;
    }
}


double von_mises_profile::eval_integrated_sample_slow(double phi0, double phi1, double amplitude) const
{
    if (_unlikely(phi0 >= phi1))
	throw runtime_error("von_mises_profile::integrate_sample_slow(): expected phi0 < phi1");

    const double *rho = &detrended_profile[0];
    const double mf = detrend ? 0.0 : mean_flux;

    double p0 = phi0 * internal_nphi;
    double p1 = phi1 * internal_nphi;

    ssize_t i0 = round_down(p0);
    ssize_t i1 = round_up(p1);

    double num = 0.0;
    double den = 0.0;

    while (i0 < i1) {
	ssize_t j = xmod(i0, internal_nphi);
	double u0 = max(p0-i0, 0.0);
	double u1 = min(p1-i0, 1.0);
	double u = 0.5 * (u0+u1);

	num += (u1-u0) * _linterp(u, rho[j], rho[j+1]);
	den += (u1-u0);
	i0++;
    }

    simpulse_assert(den > 0.0);

    return amplitude * (num/den + mf);
}


double von_mises_profile::point_eval(double phi, double amplitude) const
{
    const double offset = detrend ? mean_flux : 0.0;
    return amplitude * (_vm_profile(kappa, phi) - offset);
}


// -------------------------------------------------------------------------------------------------
//
// Signal-to-noise calculation


// Helper method: returns expectation value <rho^2>, where expectation value is taken over phi,
// and boxcar filtering is applied.  The 'dphi' argument is the phase advance per boxcar.

double von_mises_profile::_get_rho2(double dphi) const
{
    double ret = profile_fft[0];  // can be zero, if detrend=true

    for (int m = 1; m < internal_nphi2; m++) {
	double rho_m = profile_fft[m] * bessj0(M_PI * m * dphi);
	ret += 2 * square(rho_m);
    }

    return ret;
}


double von_mises_profile::get_single_pulse_signal_to_noise(double dt_sample, double pulse_freq, double sample_rms) const
{
    if (_unlikely(dt_sample <= 0.0))
	throw runtime_error("von_mises_profile::get_single_pulse_signal_to_noise(): expected dt_sample > 0.0");
    if (_unlikely(pulse_freq <= 0.0))
	throw runtime_error("von_mises_profile::get_single_pulse_signal_to_noise(): expected pulse_freq > 0.0");
    if (_unlikely(sample_rms <= 0.0))
	throw runtime_error("von_mises_profile::get_single_pulse_signal_to_noise(): expected pulse_freq > 0.0");

    double rho2 = _get_rho2(pulse_freq * dt_sample);
    return sqrt(rho2 / (pulse_freq * dt_sample)) / sample_rms;
}


double von_mises_profile::get_multi_pulse_signal_to_noise(double total_time, double dt_sample, double pulse_freq, double sample_rms) const
{
    if (_unlikely(total_time <= 0.0))
	throw runtime_error("von_mises_profile::get_multi_pulse_signal_to_noise(): expected total_time > 0.0");
    if (_unlikely(dt_sample <= 0.0))
	throw runtime_error("von_mises_profile::get_multi_pulse_signal_to_noise(): expected dt_sample > 0.0");
    if (_unlikely(pulse_freq <= 0.0))
	throw runtime_error("von_mises_profile::get_multi_pulse_signal_to_noise(): expected pulse_freq > 0.0");
    if (_unlikely(sample_rms <= 0.0))
	throw runtime_error("von_mises_profile::get_multi_pulse_signal_to_noise(): expected pulse_freq > 0.0");

    double rho2 = _get_rho2(pulse_freq * dt_sample);
    return sqrt((total_time * rho2) / dt_sample) / sample_rms;
}


template<typename T>
void von_mises_profile::get_profile_fft(T *out, int nout) const
{
    if (!out)
	throw runtime_error("simpulse::von_mises_profile::get_profile_fft(): NULL pointer passed as 'out' argument");
    if (nout <= 0)
	throw runtime_error("simpulse::von_mises_profile::get_profile_fft(): the 'nout' argument must be > 0");

    int m = min(nout, internal_nphi2);

    for (int i = 0; i < m; i++)
	out[i] = T(profile_fft[i]);

    for (int i = m; i < nout; i++)
	out[i] = T(0);
}


// Instantiate von_mises_profile::get_profile_fft() for T=float,double.
template void von_mises_profile::get_profile_fft(float *out, int nout) const;
template void von_mises_profile::get_profile_fft(double *out, int nout) const;


}  // namespace simpulse
