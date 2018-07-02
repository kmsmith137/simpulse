#include <iostream>
#include "../include/simpulse/pulsar_phase_models.hpp"
#include "../include/simpulse/pulsar_profiles.hpp"
#include "../include/simpulse/internals.hpp"

#define DEBUG 0

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


static ssize_t round_up_to_fft_friendly_value(ssize_t n)
{
    sp_assert(n >= 2);
    
    // Try (2^m) and (2^m 3).
    ssize_t n1 = 1 << round_up(log2(n-0.1));
    ssize_t n2 = 3 * (1 << round_up(log2(n-0.1) - log2(3.)));
    
    sp_assert(n1 >= n && n1 < 2*n);
    sp_assert(n2 >= n && n2 < 2*n);

    return min(n1, n2);
}


static int choose_internal_nphi(double duty_cycle, int min_internal_nphi=0)
{
    sp_assert2(duty_cycle > 0.0, "simpulse::von_mises_profile: duty_cycle must be > 0");
    sp_assert2(duty_cycle < 0.5, "simpulse::von_mises_profile: duty_cycle must be < 0.5");

    // There is no fundamental reason for these restrictions, but violating them is probably unintentional
    // and a symptom of a bug!  Note that in the limit duty_cyle->0, our implementation of the periodic pulse
    // is inefficient, so if simulating lower duty cycles is of interest, it might be best to revisit the code.

    sp_assert2(duty_cycle >= 1.0e-4, "simpulse::von_mises_profile: we currently don't support duty cycles < 1.0e-4");
    sp_assert2(min_internal_nphi >= 0, "simpulse::von_mises_profile: 'min_internal_nphi' argument must be >= 0");
    sp_assert2(min_internal_nphi <= 65536, "simpulse::von_mises_profile: we currently don't support nphi >= 65536");

    // This value suffices to simulate the pulsar to ~1% accuracy or better.
    int nphi = round_up(30./duty_cycle);

    nphi = round_up_to_fft_friendly_value(nphi);
    nphi = max(nphi, min_internal_nphi);
    return nphi;
}


static int choose_internal_phi_block_size(int requested_block_size)
{
    if (requested_block_size == 0)
	return 256;  // default value

    sp_assert2(requested_block_size > 0, "simpulse::von_mises_profile: internal_phi_block_size must be >= 0");
    sp_assert2(requested_block_size <= 4096, "simpulse::von_mises_profile: internal_phi_block_size must be <= 4096");

    return requested_block_size;
}    


von_mises_profile::von_mises_profile(double duty_cycle_, bool detrend_, int min_internal_nphi, int internal_phi_block_size_)
    : duty_cycle(duty_cycle_),
      detrend(detrend_),
      internal_nphi(choose_internal_nphi(duty_cycle_, min_internal_nphi)),   // 'duty_cycle', 'min_internal_nphi' sanity-checked here
      internal_phi_block_size(choose_internal_phi_block_size(internal_phi_block_size_)),
      kappa(log(2.0) / (2 * square(sin(M_PI*duty_cycle/2.)))),
      internal_nphi2(internal_nphi/2 + 1),
      _peak_flux(1.0),
      _mf_multiplier(0.0)   // to be initialized below
{
    this->detrended_profile.resize(internal_nphi+1, 0.0);
    this->detrended_profile_antider.resize(internal_nphi+1, 0.0);
    this->profile_fft.resize(internal_nphi2, 0.0);

    double *rho = &detrended_profile[0];
    double *rho_a = &detrended_profile_antider[0];

    for (int iphi = 0; iphi < internal_nphi; iphi++) {
	rho[iphi] = _vm_profile(kappa, iphi / double(internal_nphi));
	_mf_multiplier += rho[iphi] / internal_nphi;
    }

    // Compute 'profile_fft'.  (Don't forget detrending correction!)

    vector<complex<double>> ctmp(internal_nphi2, 0.0);

    fftw_plan plan = fftw_plan_dft_r2c_1d(internal_nphi, rho, reinterpret_cast<fftw_complex *> (&ctmp[0]), FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    for (int iphi = 0; iphi < internal_nphi2; iphi++)
	profile_fft[iphi] = ctmp[iphi].real() / internal_nphi;

    sp_assert(fabs(profile_fft[0] - _mf_multiplier) < 1.0e-13);

    if (detrend)
	profile_fft[0] = 0.0;

    // Complete calculation of 'detrended_profile' (by subtracting mean flux and padding by one), and 'detrended_profile_antider'.

    for (int iphi = 0; iphi < internal_nphi; iphi++)
	rho[iphi] -= _mf_multiplier;

    rho[internal_nphi] = rho[0];  // pad by one (with periodicity)

    for (int iphi = 0; iphi < internal_nphi; iphi++) {
	double x = 0.5 * (rho[iphi] + rho[iphi+1]);
	rho_a[iphi+1] = rho_a[iphi] + x;
    }

    sp_assert(fabs(rho_a[internal_nphi] - 0.0) < 1.0e-11);
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
	sp_assert((i >= 0) && (i < n));
	sp_assert((u >= -1.0e-12) && (u <= 1+1.0e-12));
	sp_assert(fabs(p - j*n - i - u) < (1.0e-12 * (1 + fabs(p))));
#endif

	rho0 = rho[i];
	rho1 = rho[i+1];
    }

    inline double _linterp(double uu)
    {
	return (1-uu)*rho0 + uu*rho1;
    }
};


// This is really ugly, is there a better way?
//
// The inline function _integrate_samples() forms the "body" of eval_integrated_samples() 
// or add_integrated_samples().  It needs a lot of arguments since it is not a member function!
//
// The template argument Tout is an inline "output" class which is called in the innermost loop, 
// and updates the output array.

template<typename Tout>
inline void _integrate_samples(Tout &out, double t0, double t1, ssize_t nt, const phase_model_base &pm,
			       double peak_flux, double mean_flux, ssize_t internal_nphi, const double *rho,
			       const double *rho_a, ssize_t phi_block_size, double *phi_tmp)
{
    // Number of time samples simulated so far.
    ssize_t nt0 = 0;

    while (nt0 < nt) {
	// Simulate one "block" of time samples i = nt0, ..., (nt1-1).
	ssize_t nt1 = min(nt, nt0 + phi_block_size);
	double t0_block = _linterp(nt0/double(nt), t0, t1);
	double t1_block = _linterp(nt1/double(nt), t0, t1);

	// Evaluate phase model in this block.
	// Note that the number of 'phi' values is nphi=(nt1-nt0+1), not nphi=(nt1-nt0).
	// As a consequence, phase_model_base::eval_phi_sequence() is always called with nphi > 1 as expected.
	pm.eval_phi_sequence(t0_block, t1_block, nt1-nt0+1, &phi_tmp[0]);
	
	ihelper h0(phi_tmp[0], internal_nphi, rho);
	
	for (ssize_t it = nt0; it < nt1; it++) {
	    // Loop invariant: at top, h0 corresponds to time it, and h1 corresponds to time (it+1).
	    ihelper h1(phi_tmp[it-nt0+1], internal_nphi, rho);
	    double dp = h1.p - h0.p;
	    double integral;

	    sp_assert2(dp > 0.0, "von_mises_profile::eval_integrated_samples(): phase model is not monotone increasing?!");

	    // I put way too much thought into writing the logic below in a way which minimizes roundoff error!

	    if (h0.i == h1.i)
		integral = (h1.u - h0.u) * h0._linterp(0.5 * (h0.u + h1.u));
	    else {
		integral = (1.0 - h0.u) * h0._linterp(0.5*h0.u + 0.5);
		integral += rho_a[h1.i] - rho_a[h0.i+1];
		integral += h1.u * h1._linterp(0.5*h1.u);
	    }

	    out.process_sample(it, peak_flux * (integral/dp) + mean_flux);
	    h0 = h1;
	}

	nt0 = nt1;
    }    
}


// This helper class is used as the "output class" in _integrate_samples(),
// when _integrate_samples() is called from eval_integrated_samples().

template<typename T>
struct _array_overwriter {
    T *out;
    _array_overwriter(T *out_) : out(out_) { }
    inline void process_sample(ssize_t it, double val) { out[it] = T(val); }
};


// This helper class is used as the "output class" in _integrate_samples(),
// when _integrate_samples() is called from eval_integrated_samples().

template<typename T>
struct _array_accumulator {
    T *out;
    _array_accumulator(T *out_) : out(out_) { }
    inline void process_sample(ssize_t it, double val) { out[it] += T(val); }
};


template<typename T>
void von_mises_profile::eval_integrated_samples(T *out, double t0, double t1, ssize_t nt, const phase_model_base &pm) const
{
    sp_assert2(nt > 0, "simpulse::von_mises_profile::eval_integrated_samples(): expected nt > 0");
    sp_assert2(t0 < t1, "simpulse::von_mises_profile::eval_integrated_samples(): expected t0 < t1");
    sp_assert2(out, "simpulse::von_mises_profile::eval_integrated_samples(): 'out' is a null pointer");

    if (!phi_tmp.size())
	phi_tmp.resize(internal_phi_block_size+1, 0.0);   // note "+1" here

    _array_overwriter<T> hout(out);

    _integrate_samples(hout, t0, t1, nt, pm, _peak_flux,
		       detrend ? 0.0 : (_mf_multiplier * _peak_flux),
		       internal_nphi, &detrended_profile[0], 
		       &detrended_profile_antider[0], 
		       internal_phi_block_size, &phi_tmp[0]);
}


template<typename T>
void von_mises_profile::add_integrated_samples(T *out, double t0, double t1, ssize_t nt, const phase_model_base &pm) const
{
    sp_assert2(nt > 0, "simpulse::von_mises_profile::add_integrated_samples(): expected nt > 0");
    sp_assert2(t0 < t1, "simpulse::von_mises_profile::add_integrated_samples(): expected t0 < t1");
    sp_assert2(out, "simpulse::von_mises_profile::add_integrated_samples(): 'out' is a null pointer");

    if (!phi_tmp.size())
	phi_tmp.resize(internal_phi_block_size+1, 0.0);   // note "+1" here

    _array_accumulator<T> hout(out);

    _integrate_samples(hout, t0, t1, nt, pm, _peak_flux,
		       detrend ? 0.0 : (_mf_multiplier * _peak_flux),
		       internal_nphi, &detrended_profile[0], 
		       &detrended_profile_antider[0], 
		       internal_phi_block_size, &phi_tmp[0]);
}


double von_mises_profile::eval_integrated_sample_slow(double phi0, double phi1) const
{
    sp_assert2(phi0 < phi1, "simpulse::von_mises_profile::integrate_sample_slow(): expected phi0 < phi1");

    const double *rho = &detrended_profile[0];
    const double pf = _peak_flux;
    const double mf = detrend ? 0.0 : (_mf_multiplier * _peak_flux);

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

    sp_assert(den > 0.0);

    return pf * (num/den) + mf;
}


double von_mises_profile::eval_instantaneous(double phi) const
{
    double offset = detrend ? (_mf_multiplier * _peak_flux) : 0.0;
    return _peak_flux * _vm_profile(kappa, phi) - offset;
}


// -------------------------------------------------------------------------------------------------
//
// Member functions which get/set the normalization.


// Helper method: returns expectation value <rho^2>, where expectation value is taken over phi,
// and boxcar filtering is applied.  The 'dphi' argument is the phase advance per boxcar.

double von_mises_profile::_get_rho2(double dphi) const
{
    double ret = square(profile_fft[0]);  // can be zero, if detrend=true

    for (int m = 1; m < internal_nphi2; m++) {
	double rho_m = profile_fft[m] * bessj0(M_PI * m * dphi);
	ret += 2 * square(rho_m);
    }

    return square(_peak_flux) * ret;
}


double von_mises_profile::get_peak_flux() const
{
    return _peak_flux;
}

double von_mises_profile::get_mean_flux() const
{
    return _mf_multiplier * _peak_flux;
}


double von_mises_profile::get_single_pulse_signal_to_noise(double dt_sample, double pulse_freq, double sample_rms) const
{
    sp_assert2(dt_sample > 0.0, "simpulse::von_mises_profile::get_single_pulse_signal_to_noise(): expected dt_sample > 0.0");
    sp_assert2(pulse_freq > 0.0, "simpulse::von_mises_profile::get_single_pulse_signal_to_noise(): expected pulse_freq > 0.0");
    sp_assert2(sample_rms > 0.0, "simpulse::von_mises_profile::get_single_pulse_signal_to_noise(): expected pulse_freq > 0.0");

    double rho2 = _get_rho2(pulse_freq * dt_sample);
    return sqrt(rho2 / (pulse_freq * dt_sample)) / sample_rms;
}


double von_mises_profile::get_multi_pulse_signal_to_noise(double total_time, double dt_sample, double pulse_freq, double sample_rms) const
{
    sp_assert2(total_time > 0.0, "simpulse::von_mises_profile::get_single_pulse_signal_to_noise(): expected total_time > 0.0");
    sp_assert2(dt_sample > 0.0, "simpulse::von_mises_profile::get_single_pulse_signal_to_noise(): expected dt_sample > 0.0");
    sp_assert2(pulse_freq > 0.0, "simpulse::von_mises_profile::get_single_pulse_signal_to_noise(): expected pulse_freq > 0.0");
    sp_assert2(sample_rms > 0.0, "simpulse::von_mises_profile::get_single_pulse_signal_to_noise(): expected pulse_freq > 0.0");
    sp_assert2(total_time > dt_sample, "simpulse::von_mises_profile::get_single_pulse_signal_to_noise(): expected total_time > dt_sample");

    double rho2 = _get_rho2(pulse_freq * dt_sample);
    return sqrt((total_time * rho2) / dt_sample) / sample_rms;
}


void von_mises_profile::set_peak_flux(double peak_flux)
{
    this->_peak_flux = peak_flux;
}

void von_mises_profile::set_mean_flux(double mean_flux)
{
    this->_peak_flux = mean_flux / _mf_multiplier;
}

void von_mises_profile::set_single_pulse_signal_to_noise(double snr, double dt_sample, double pulse_freq, double sample_rms)
{
    double current_snr = get_single_pulse_signal_to_noise(dt_sample, pulse_freq, sample_rms);
    this->_peak_flux *= (snr / current_snr);
}

void von_mises_profile::set_multi_pulse_signal_to_noise(double snr, double total_time, double dt_sample, double pulse_freq, double sample_rms)
{
    double current_snr = get_multi_pulse_signal_to_noise(total_time, dt_sample, pulse_freq, sample_rms);
    this->_peak_flux *= (snr / current_snr);
}


template<typename T>
void von_mises_profile::get_profile_fft(T *out, int nout) const
{
    sp_assert2(nout > 0, "simpulse::von_mises_profile::get_profile_fft(): the 'nout' argument must be > 0");
    sp_assert2(out, "simpulse::von_mises_profile::get_profile_fft(): NULL pointer passed as 'out' argument");

    int m = min(nout, internal_nphi2);

    for (int i = 0; i < m; i++)
	out[i] = T(_peak_flux * profile_fft[i]);

    for (int i = m; i < nout; i++)
	out[i] = T(0);
}


string von_mises_profile::str() const
{
    stringstream ss;
    ss << "simpulse.von_mises_profile(duty_cycle=" << duty_cycle << ", detrend=" 
       << (detrend ? "true" : "false") << ", peak_flux=" << _peak_flux << ")";
    return ss.str();
}


// -------------------------------------------------------------------------------------------------


#define INSTANTIATE_TEMPLATES(T) \
    template void von_mises_profile::eval_integrated_samples(T *out, double t0, double t1, ssize_t nt, const phase_model_base &pm) const; \
    template void von_mises_profile::add_integrated_samples(T *out, double t0, double t1, ssize_t nt, const phase_model_base &pm) const; \
    template void von_mises_profile::get_profile_fft(T *out, int nout) const

INSTANTIATE_TEMPLATES(float);
INSTANTIATE_TEMPLATES(double);


}  // namespace simpulse
