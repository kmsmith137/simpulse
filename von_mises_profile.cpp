#include "simpulse_internals.hpp"

#define DEBUG 0

using namespace std;


namespace simpulse {
#if 0
};  // pacify emacs c-mode!
#endif


von_mises_profile::von_mises_profile(double duty_cycle_, int min_nphi)
    : duty_cycle(duty_cycle_)
{
    if (duty_cycle <= 0.0)
	throw runtime_error("simpulse::von_mises_profile constructor: 'duty_cycle' argument must be >= 0");

    // Note: the von Mises profile actually makes mathematical sense for D < 1
    if (duty_cycle >= 0.5)
	throw runtime_error("simpulse::von_mises_profile constructor: 'duty_cycle' argument must be < 0.5");

    // No fundamental reason for these restrictions, but violating them is probably unintentional!
    // Note that in the limit duty_cyle->0, this implementation is inefficient, and it would be better
    // to supply an alternate implementation similar to single_pulse.
    if (duty_cycle < 1.0e-4)
	throw runtime_error("simpulse::von_mises_profile constructor: we currently don't support duty cycles < 1.0e-4");
    if (min_nphi > 65536)
	throw runtime_error("simpulse::von_mises_profile constructor: we currently don't support nphi >= 65536");

    this->nphi = max(default_nphi(duty_cycle), min_nphi);
    this->nphi2 = nphi/2 + 1;
    this->profile_fft.resize(nphi2, 0.0);
    this->profile_antider.resize(nphi+1, 0.0);

    vector<double> profile(nphi);
    vector<complex<double> > tmp_fft(nphi2, 0.0);

    //
    // The von Mises profile is 
    //   f(phi) = exp[ -2 kappa sin(phi/2)^2 ]
    //
    // where kappa is related to the duty cycle D by
    //   kappa = log(2) / (2 sin^2(pi*D/2))
    //
    double kappa = log(2.0) / (2 * square(sin(M_PI*duty_cycle/2.)));

    for (int iphi = 0; iphi < nphi; iphi++) {
        double phi = 2*M_PI * iphi / double(nphi);
	profile[iphi] = exp(-2. * kappa * square(sin(phi/2.)));
    }

    // We normalize so that the mean intensity is equal to one.
    double mean = 0.0;
    for (int iphi = 0; iphi < nphi; iphi++)
	mean += profile[iphi] / nphi;
    for (int iphi = 0; iphi < nphi; iphi++)
	profile[iphi] /= mean;

    // Cumulative sum.  Note that profile_antider[0] was initialized to zero above.
    for (int iphi = 0; iphi < nphi; iphi++) {
	double p1 = profile[(iphi)];
	double p2 = profile[(iphi+1) % nphi];
	profile_antider[iphi+1] = profile_antider[iphi] + (p1+p2)/2. * (2*M_PI)/nphi;
    }

    fftw_plan plan = fftw_plan_dft_r2c_1d(nphi, &profile[0], reinterpret_cast<fftw_complex *> (&tmp_fft[0]), FFTW_ESTIMATE);
    fftw_execute(plan);

    for (int iphi = 0; iphi < nphi2; iphi++)
	profile_fft[iphi] = tmp_fft[iphi].real() / nphi;

    simpulse_assert(fabs(profile_antider[nphi] - 2*M_PI) < 1.0e-12);
    simpulse_assert(fabs(profile_fft[0] - 1.0) < 1.0e-12);

    // Convervative cutoff (pulsar frequency must be sub-Nyquist)
    this->ntab_xmax = M_PI + 1.0e-9;

    // Also a conservative choice
    this->ntab_size = round_up(16. / duty_cycle);
    this->ntab.resize(ntab_size, 0.0);

    for (int ix = 0; ix < ntab_size; ix++) {
	double x = ix * ntab_xmax / double(ntab_size-1);

	for (int n = 1; n < nphi2; n++) {
	    double w = (2*n == nphi) ? 1.0 : 2.0;
	    ntab[ix] += w * square(profile_fft[n] * bessj0(n*x/2.));
	}
    }
}
    

// Helper for von_mises_profile::eval().
// We write phi = 2*pi/nphi * (i+t), where 0 <= i < nphi and 0 <= t <= 1.
inline void dissect_phi(double phi, int nphi, const double *antider, int &i, double &t, double &a)
{
    t = fmod(phi,2*M_PI) / (2*M_PI) * nphi;
    t += 2 * nphi;   // force positive

    i = int(t);
    t -= double(i);
    i %= nphi;
    
#if DEBUG
    double dphi = 2*M_PI/nphi * (i+t) - phi;
    simpulse_assert(fabs(sin(dphi/2.)) < 1.0e-10);   // a way of asserting that dphi is divisible by 2pi
    simpulse_assert(i >= 0 && i < nphi);
    simpulse_assert(t > -1.0e-10);
    simpulse_assert(t < 1.0 + 1.0e-10);
#endif

    a = (1.0-t) * antider[i] + t * antider[i+1];
}


void von_mises_profile::eval(int nt, double *out, const double *phi, bool detrend, double amplitude) const
{
    const double detrend_corr = detrend ? 1.0 : 0.0;
    const double *antider = &profile_antider[0];

    int i0;
    double t0, a0;
    dissect_phi(phi[0], nphi, antider, i0, t0, a0);

    for (int it = 0; it < nt; it++) {
	// Loop invariant: at top, 
	//   phi[it] = 2*pi/nphi * (i0+t0),
	//   a0 = value of profile_antider array, interpolated at phi[it]

	int i1;
	double t1, a1;
	dissect_phi(phi[it+1], nphi, antider, i1, t1, a1);

	double dphi = phi[it+1] - phi[it];
	double wraparounds = dphi - (i1+t1-i0-t0) * (2*M_PI/nphi);   // multiple of (2pi)
	double integral = a1 - a0 + wraparounds;

	out[it] = amplitude * (integral/dphi - detrend_corr);

	// Maintain loop invariant
	i0 = i1;
	t0 = t1;
	a0 = a1;
    }
}


double von_mises_profile::get_normalization(double eta, double omega, double tsamp, double T, bool detrend) const
{
    double n0 = detrend ? 0.0 : 1.0;
    double x = omega * tsamp;
    double y = x / ntab_xmax * (ntab_size-1);
    
    if ((y < 0) || (y > ntab_size-1))
	throw runtime_error("von_mises_profile::get_normalization: product omega*tsamp is negative or too large");
    
    int iy = int(y);

    // just being paranoid
    iy = max(iy,0);
    iy = min(iy,ntab_size-2);

    double t = y - double(iy);
    double n_interp = (1-t)*ntab[iy] + t*ntab[iy+1];
    return eta / sqrt(T * (n_interp + n0));
}


// static member function
int von_mises_profile::default_nphi(double duty_cycle)
{
    if ((duty_cycle <= 0) || (duty_cycle > 0.5))
	throw runtime_error("von_mises_profile::default_nphi(): invalid duty cycle");
    
    // See comment in constructor
    if (duty_cycle < 1.0e-4)
	throw runtime_error("simpulse::von_mises_profile constructor: we currently don't support duty cycles < 1.0e-4");

    // Ultra-conservative value
    // FIXME I'd like to understand better how to choose this
    return round_up(30./duty_cycle);
}


}  // namespace simpulse
