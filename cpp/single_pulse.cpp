#include "../include/simpulse.hpp"
#include "../include/simpulse/internals.hpp"

using namespace std;


namespace simpulse {
#if 0
};  // pacify emacs c-mode!
#endif


single_pulse::single_pulse(int pulse_nt_, int nfreq_, double freq_lo_MHz_, double freq_hi_MHz_, 
			   double dm_, double sm_, double intrinsic_width_, double fluence_, 
			   double spectral_index_, double undispersed_arrival_time_)
    : pulse_nt(pulse_nt_), nfreq(nfreq_), freq_lo_MHz(freq_lo_MHz_), freq_hi_MHz(freq_hi_MHz_), 
      dm(dm_), sm(sm_), intrinsic_width(intrinsic_width_), fluence(fluence_), 
      spectral_index(spectral_index_), undispersed_arrival_time(undispersed_arrival_time_)
{
    simpulse_assert(pulse_nt >= 64);   // using fewer time samples than this is probably a mistake
    simpulse_assert(nfreq > 0);
    simpulse_assert(freq_lo_MHz > 0.0);
    simpulse_assert(freq_hi_MHz > freq_lo_MHz);

    simpulse_assert(dm >= 0.0);
    simpulse_assert(sm >= 0.0);
    simpulse_assert(intrinsic_width >= 0.0);
    simpulse_assert(fluence >= 0.0);

    // Implementing delta function pulses wouldn't be a big deal, but creates corner cases
    // and so far I haven't seen a strong reason to implement it.
    if ((dm == 0.0) && (sm == 0.0) && (intrinsic_width == 0.0))
	throw runtime_error("single_pulse: delta function pulse (dm=sm=width=0) is currently not allowed");

    this->pulse_t0.resize(nfreq, 0.0);
    this->pulse_t1.resize(nfreq, 0.0);
    this->pulse_freq_wt.resize(nfreq, 0.0);
    this->pulse_cumsum.resize(nfreq * (pulse_nt+1), 0.0);

    this->_compute_freq_wt();

    int nfft = 2*pulse_nt;
    int nfft2 = nfft/2 + 1;

    double *bufr = checked_fftw_malloc<double> (nfft);
    complex<double> *bufc = checked_fftw_malloc<complex<double> > (nfft2);
    fftw_plan plan = fftw_plan_dft_c2r_1d(nfft, reinterpret_cast<fftw_complex *> (bufc), bufr, FFTW_ESTIMATE);

    // The following loop synthesizes the pulse.
    // We sample the pulse at time t_i = t0 + (i+0.5)*(t1-t0)/pulse_nt, where 0 <= i < pulse_nt.

    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	double nu_lo = freq_lo_MHz + (ifreq) * (freq_hi_MHz - freq_lo_MHz) / (double)nfreq;
	double nu_hi = freq_lo_MHz + (ifreq+1) * (freq_hi_MHz - freq_lo_MHz) / (double)nfreq;
	double nu_c = (nu_lo + nu_hi) / 2.;
	
	double dm_delay0 = dispersion_delay(dm, nu_hi);
	double dm_delay1 = dispersion_delay(dm, nu_lo);
	double dm_width = dm_delay1 - dm_delay0;
	double tscatt = scattering_time(sm, nu_c);

	double t0 = dm_delay0 - 0.1*dm_width - 4.*intrinsic_width - tscatt;
	double t1 = dm_delay1 + 0.1*dm_width + 4.*intrinsic_width + 10.*tscatt;
	double tc = (dm_delay0 + dm_delay1) / 2.;         // pulse center in channel
	double dt = tc - (t0 + (t1-t0)/(2.*pulse_nt));    // pulse center relative to first sample

	simpulse_assert(t0 < t1);
	this->pulse_t0[ifreq] = t0;
	this->pulse_t1[ifreq] = t1;

	double *p = &this->pulse_cumsum[ifreq * (pulse_nt+1)];
	double omega0 = 2*M_PI * (double)pulse_nt / (double)nfft / (t1-t0);

	for (int j = 0; j < nfft2; j++) {
	    double omega = omega0 * j;
	    
	    // Fourier transform of pulse
	    // Note: FFTW sign convention is T(x) = sum_k T(k) e^{ik.x}

	    bufc[j] = exp(-square(intrinsic_width*omega)/2.);            // Gaussian pulse
	    bufc[j] *= bessj0(dm_width * omega / 2.);                    // dispersion broadening
	    bufc[j] /= complex<double>(1., tscatt * omega);              // scatter broadening
	    bufc[j] *= complex<double>(cos(dt*omega), -sin(dt*omega));   // phase shift to center
	}

	// FFT bufc -> bufr
	fftw_execute(plan);
	
	// Evaluate cumsum, cleaning up samples which are negative due to discretization effects
	for (int it = 0; it < pulse_nt; it++)
	    p[it+1] = p[it] + max(bufr[it], 0.0);

	// Normalize to sum=1
	for (int it = 0; it < pulse_nt+1; it++)
	    p[it] = p[it] / p[pulse_nt];
    }

    fftw_free(bufr);
    fftw_free(bufc);
    fftw_destroy_plan(plan);

    // Initialize min_t0, max_t1, max_dt

    this->min_t0 = pulse_t0[0];
    this->max_t1 = pulse_t1[0];
    this->max_dt = pulse_t1[0] - pulse_t0[0];

    for (int ifreq = 1; ifreq < nfreq; ifreq++) {
	this->min_t0 = min(min_t0, pulse_t0[ifreq]);
	this->max_t1 = max(max_t1, pulse_t1[ifreq]);
	this->max_dt = max(max_dt, pulse_t1[ifreq] - pulse_t0[ifreq]);
    }
}


void single_pulse::_compute_freq_wt()
{
    if ((spectral_index < -20.1) || (spectral_index > 20.1))
	throw runtime_error("single_pulse::spectral_index set to 'extreme' value " + xto_string(spectral_index) + ", this is currently disallowed");

    double nu0 = (freq_lo_MHz + freq_hi_MHz) / 2.0;

    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	double nu = freq_lo_MHz + (ifreq+0.5) * (freq_hi_MHz - freq_lo_MHz) / (double)nfreq;
	this->pulse_freq_wt[ifreq] = pow(nu/nu0, spectral_index);
    }
}


void single_pulse::set_fluence(double fluence_)
{
    this->fluence = fluence_;
    simpulse_assert(fluence_ >= 0.0);
}

void single_pulse::set_spectral_index(double spectral_index_)
{
    this->spectral_index = spectral_index_;
    this->_compute_freq_wt();
}

void single_pulse::set_undispersed_arrival_time(double undispersed_arrival_time_)
{
    this->undispersed_arrival_time = undispersed_arrival_time_;
}

void single_pulse::get_endpoints(double &t0, double &t1) const
{
    t0 = this->undispersed_arrival_time + this->min_t0;
    t1 = this->undispersed_arrival_time + this->max_t1;
}


// Helper function called by _add_pulse_to_frequency_channel().  
// The 'arr' arg is an array of length (pulse_nt+1).
// The 's' arg is time in "sample coords", i.e. elements of 'arr' correspond to times s=0,1,...,pulse_nt.
inline double _interpolate_cumsum(int pulse_nt, const double *arr, double s)
{
    if (s < 1.0e-10)
	return 0.0;
    if (s > pulse_nt - 1.0e-10)
	return arr[pulse_nt];
    
    int is = (int)s;
    double ds = s - is;
    simpulse_assert(is >= 0 && is < pulse_nt);
    
    return (1-ds)*arr[is] + (ds)*arr[is+1];
}


//
// Helper function called by single_pulse::add_to_timestream() and single_pulse::get_signal_to_noise().
//
// Note to self: when I implement pulsars, I'm imagining that the "under-the-hood" parts of 
// struct single_pulse will be factored into a new class, with a member functions which are
// similar to _interpolate_cumsum() and _add_pulse_to_frequency_channel().
//
template<typename T>
void single_pulse::_add_pulse_to_frequency_channel(T *out, double out_t0, double out_t1, int out_nt, int ifreq, double weight) const
{
    simpulse_assert(out);
    simpulse_assert(out_nt > 0);
    simpulse_assert(out_t0 < out_t1);
    simpulse_assert(ifreq >= 0 && ifreq < nfreq);

    // Convert input times to "sample coords"
    double s0 = pulse_nt * (out_t0 - undispersed_arrival_time - pulse_t0[ifreq]) / (pulse_t1[ifreq] - pulse_t0[ifreq]);
    double s1 = pulse_nt * (out_t1 - undispersed_arrival_time - pulse_t0[ifreq]) / (pulse_t1[ifreq] - pulse_t0[ifreq]);

    if ((s0 >= pulse_nt) || (s1 <= 0))
	return;
    
    double out_dt = (out_t1 - out_t0) / out_nt;
    double w = weight * fluence * pulse_freq_wt[ifreq] / out_dt;
    const double *cs = &pulse_cumsum[ifreq*(pulse_nt+1)];

    for (int it = 0; it < out_nt; it++) {
	double a = _interpolate_cumsum(pulse_nt, cs, s0 + (it)*(s1-s0)/(double)out_nt);
	double b = _interpolate_cumsum(pulse_nt, cs, s0 + (it+1)*(s1-s0)/(double)out_nt);
	out[it] += w * (b-a);
    }
}


template<typename T>
void single_pulse::add_to_timestream(T *out, double out_t0, double out_t1, int out_nt, int stride, double weight) const
{
    if (stride == 0)
	stride = out_nt;

    simpulse_assert(out);
    simpulse_assert(out_nt > 0);
    simpulse_assert(out_t0 < out_t1);
    simpulse_assert(abs(stride) >= out_nt);   // allow negative stride as explained in simpulse.hpp
    
    // Return early if data does not overlap pulse
    if (out_t0 > undispersed_arrival_time + max_t1)
	return;
    if (out_t1 < undispersed_arrival_time + min_t0)
	return;

    for (int ifreq = 0; ifreq < nfreq; ifreq++)
	_add_pulse_to_frequency_channel(out + ifreq*stride, out_t0, out_t1, out_nt, ifreq, weight);
}


double single_pulse::get_signal_to_noise(double sample_dt, double sample_rms, double sample_t0) const
{
    simpulse_assert(sample_dt > 0.0);
    simpulse_assert(sample_rms > 0.0);

    int nsamp_max = (int)(max_dt/sample_dt) + 3;
    vector<double> buf(nsamp_max, 0.0);

    double acc = 0.0;

    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	// Range of samples spanned by pulse
	double s0 = (undispersed_arrival_time + pulse_t0[ifreq] - sample_t0) / sample_dt;
	double s1 = (undispersed_arrival_time + pulse_t1[ifreq] - sample_t0) / sample_dt;

	ssize_t j = round_down(s0);
	ssize_t k = round_up(s1);
	simpulse_assert(k-j <= nsamp_max);

	memset(&buf[0], 0, nsamp_max * sizeof(double));
	_add_pulse_to_frequency_channel(&buf[0], sample_t0 + j*sample_dt, sample_t0 + k*sample_dt, k-j, ifreq, 1.0);
	
	for (ssize_t i = 0; i < k-j; i++)
	    acc += buf[i]*buf[i];
    }

    return sqrt(acc) / sample_rms;
}


double single_pulse::get_signal_to_noise(double sample_dt, const double *sample_rms, const double *channel_weights, double sample_t0) const
{
    if (sample_rms == nullptr)
	throw runtime_error("simpulse::single_pulse::get_signal_to_noise(): sample_rms pointer was NULL");
    if (sample_dt <= 0.0)
	throw runtime_error("simpulse::single_pulse::get_signal_to_noise(): sample_dt must be positive");

    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	if (sample_rms[ifreq] < 0.0)
	    throw runtime_error("simpulse::single_pulse::get_signal_to_noise(): all values of 'sample_rms' array must be nonnegative");
	if (channel_weights && (channel_weights[ifreq] < 0.0))
	    throw runtime_error("simpulse::single_pulse::get_signal_to_noise(): all values of 'channel_weights' array must be nonnegative");
    }

    vector<double> wtmp;

    if (channel_weights == nullptr) {
	wtmp.resize(nfreq);
	channel_weights = &wtmp[0];

	for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	    if (sample_rms[ifreq] <= 0.0)
		throw runtime_error("simpulse::single_pulse::get_signal_to_noise(): if 'channel_weights' pointer is unspecified, then all sample_rms values must be positive");
	    wtmp[ifreq] = 1.0 / (sample_rms[ifreq] * sample_rms[ifreq]);
	}
    }

    int nsamp_max = (int)(max_dt/sample_dt) + 3;
    vector<double> buf(nsamp_max, 0.0);

    double sig_ampl = 0.0;
    double noise_var = 0.0;

    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	// Range of samples spanned by pulse
	double s0 = (undispersed_arrival_time + pulse_t0[ifreq] - sample_t0) / sample_dt;
	double s1 = (undispersed_arrival_time + pulse_t1[ifreq] - sample_t0) / sample_dt;

	ssize_t j = round_down(s0);
	ssize_t k = round_up(s1);
	simpulse_assert(k-j <= nsamp_max);

	memset(&buf[0], 0, nsamp_max * sizeof(double));
	_add_pulse_to_frequency_channel(&buf[0], sample_t0 + j*sample_dt, sample_t0 + k*sample_dt, k-j, ifreq, 1.0);

	double t = 0.0;
	for (ssize_t i = 0; i < k-j; i++)
	    t += square(buf[i]);
	
	sig_ampl += channel_weights[ifreq] * t;
	noise_var += square(channel_weights[ifreq] * sample_rms[ifreq]) * t;
    }

    if (noise_var <= 0.0)
	throw runtime_error("simpulse::single_pulse::get_signal_to_noise(): computed noise variance is zero, this means that too many sample_rms (or channel_weights) values were zero");

    return sig_ampl / sqrt(noise_var);
}


void single_pulse::print(ostream &os) const
{
    os << "single_pulse(pulse_nt=" << pulse_nt << ",nfreq=" << nfreq << ",freq_lo_MHz=" << freq_lo_MHz << ",freq_hi_MHz=" << freq_hi_MHz
       << ",dm=" << dm << ",sm=" << sm << ",intrinsic_width=" << intrinsic_width << ",fluence=" << fluence 
       << ",spectral_index=" << spectral_index << ",undispersed_arrival_time=" << undispersed_arrival_time << ")";
}


string single_pulse::str() const
{
    stringstream ss;
    this->print(ss);
    return ss.str();
}


#define INSTANTIATE(T) \
    template void single_pulse::add_to_timestream(T *out, double out_t0, double out_t1, int out_nt, int stride, double weight) const; \
    template void single_pulse::_add_pulse_to_frequency_channel(T *out, double out_t0, double out_t1, int out_nt, int ifreq, double weight) const

INSTANTIATE(float);
INSTANTIATE(double);


}   // namespace simpulse
