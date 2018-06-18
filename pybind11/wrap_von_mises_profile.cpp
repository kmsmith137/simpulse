#include "simpulse_pybind11.hpp"
#include "../include/simpulse/pulsar_phase_models.hpp"
#include "../include/simpulse/pulsar_profiles.hpp"
#include "../include/simpulse/internals.hpp"   // sp_assert()

#include <iostream>
#include <pybind11/numpy.h>

namespace py = pybind11;
using namespace pybind11::literals;
using namespace simpulse;
using namespace std;

namespace simpulse_pybind11 {
#if 0
}  // pacify emacs c-mode
#endif


void wrap_von_mises_profile(py::module &m)
{
    py::options options;
    options.disable_function_signatures();


    const char *doc =
	"This class represents a von Mises pulsar profile.  (Currently, the von Mises profile is the only\n"
	"periodic profile implemented in simpulse, although in the future we might implement more profiles\n"
	"and define a 'pulsar_profile' base class.)\n"
	"\n"
	"Constructor syntax::\n"
	"\n"
	"    p = simpulse.von_mises_profile(duty_cycle, detrend, min_internal_nphi=0,\n"
	"                                   internal_phi_block_size=0)\n"
	"\n"
	"where:\n"
	"\n"
	"   - The 'duty_cycle' is defined as D = (pulse full width at half maximum) / (pulse period).\n"
	"     A reasonable choice is D = 0.1 or so.\n"
	"\n"
	"   - If the boolean variable 'detrend' is true, then the mean will be subtracted from the profile.\n"
	"\n"
	"   - It is unlikely that you'll need to set the 'min_internal_nphi' constructor argument, which\n"
	"     changes the number of phase bins used internally to represent the pulse.  If set to zero, then\n"
	"     a reasonable default value will be chosen.  (The default suffices to simulate the pulsar to ~1%\n"
	"     accuracy or better.)\n"
	"\n"
	"   - It is unlikely that you'll need to set the 'internal_phi_block_size' argument, which determines\n"
	"     the block size for calls to phase_model_base.eval_phi_sequence().  If set to zero, then a\n"
	"     reasonable default value will be chosen.\n"
	"\n"
	"In order to simulate a pulsar, you need two things: a phase model and a pulse profile.  Then, to do\n"
	"the simulation, you can call either profile.eval_integrated_samples() or profile.add_integrated_samples().\n"
	"These methods take the phase model as a parameter.\n"
	"\n"
	"By default, the profile is normalized so that its peak flux is equal to 1 (before applying the detrending subtraction).\n"
	"The normalization can be changed by setting the 'peak_flux' property, setting the 'mean_flux' propery, calling\n"
	"set_single_pulse_signal_to_noise(), or calling set_multi_pulse_signal_to_noise().\n"
	"\n"
	"Mathematically, a profile is just a function rho(phi) which gives the flux 'rho' as a function of pulse phase 'phi'.\n"
	"The von Mises profile is the functional form:\n"
	"\n"
	"   rho(phi) = exp[ -2 kappa sin(pi*phi)^2 ]\n"
	"\n"
	"where kappa is a narrowness parameter, related to the duty cycle D by kappa = log(2) / (2 sin^2(pi*D/2)).";


    auto eval_integrated_samples = [](von_mises_profile &self, double t0, double t1, ssize_t nt, const phase_model_base &pm)
    {
	sp_assert2(nt > 0, "simpulse::von_mises_profile::eval_integrated_sample(): expected nt > 0");
	py::array_t<double> ret{ size_t(nt) };

	self.eval_integrated_samples(ret.mutable_data(), t0, t1, nt, pm);
	return ret;
    };


    auto get_profile_fft = [](von_mises_profile &self, int nout)
    {
	sp_assert2(nout >= 0, "simpulse::von_mises_profile::get_profile_fft(): expected nout >= 0");

	if (nout == 0)
	    nout = self.internal_nphi/2 + 1;

	py::array_t<double> ret{ size_t(nout) };

	self.get_profile_fft(ret.mutable_data(), nout);
	return ret;
    };


    py::class_<von_mises_profile>(m, "von_mises_profile", doc)
	.def(py::init<double,bool,int,int>(), "duty_cycle"_a, "detrend"_a, "min_internal_nphi"_a = 0, "internal_phi_block_size"_a = 0)

	.def_readonly("duty_cycle", &von_mises_profile::duty_cycle, 
		      "Duty cycle of the pulsar, defined as D = (pulse FWHM) / (pulse period).")

	.def_readonly("detrend", &von_mises_profile::detrend,
		      "If this boolean variable is true, then the mean will be subtracted from the pulse profile.")

	.def_readonly("internal_nphi", &von_mises_profile::internal_nphi,
		      "Number of phase bins used internally to represent the pulse.")

	.def_readonly("internal_phi_block_size", &von_mises_profile::internal_phi_block_size,
		      "Block size used internally, in calls to phase_model_base.eval_phi_sequqnece().")

	.def_readonly("kappa", &von_mises_profile::kappa,
		      "Parameter appearing in von Mises profile: rho(phi) = exp(-2 kappa sin(pi*phi)^2)")

	.def("__repr__", &von_mises_profile::str)

	.def("eval_integrated_samples", eval_integrated_samples, "t0"_a, "t1"_a, "nt"_a, "phase_model"_a,
	     "eval_integrated_samples(t0, t1, nt, phase_model) -> (1D array)\n"
	     "\n"
	     "This method can be used to simulate a pulsar in a regularly spaced sequence of time samples.\n"
	     "\n"
	     "  - The 't0' argument should be the *beginning* of the first time sample, and 't1' should be the *end*\n"
	     "    of the last time sample.  Thus, t1=t0+nt*dt, where dt is the length of a sample (not t1=t0+(nt-1)*dt).\n"
	     "\n"
	     "  - The 'phase_model' argument should be an object of type 'simpulse.phase_model_base'.\n"
	     "\n"
	     "  - The output is a 1D array of length 'nsamples', containing simulated flux values.  Note that each\n"
	     "    flux value is obtained by averaging (in time) the profile over the time sample, not instantanously\n"
	     "    sampling at the midpoint.\n"
	     "\n"
	     "Reminder: if the 'detrend' flag was specified at construction, then the simulated flux will be detrended.")

	.def("eval_integrated_sample_slow", &von_mises_profile::eval_integrated_sample_slow, "phi0"_a, "phi1"_a,
	     "eval_integrated_sample_slow(phi0, phi1) -> float\n"
	     "\n"
	     "This method is intended for debugging (hence the \"_slow\"!)\n"
	     "Returns the average flux over phase (not time) interval [phi0, phi1].")

	.def("eval_instantaneous", &von_mises_profile::eval_instantaneous, "phi"_a,
	     "eval_instantaneous(phi) -> float\n"
	     "\n"
	     "Returns the instantaneous flux evaluated at pulse phase 'phi'.\n"
	     "Reminder: if the 'detrend' flag was specified at construction, then the simulated flux will be detrended.")

	.def_property("peak_flux", 
		      &von_mises_profile::get_peak_flux, 
		      &von_mises_profile::set_peak_flux,
		      "Peak flux of the pulse (settable property).\n\n"
		      "Note: we define the peak flux with no detrending correction applied (even if the 'detrend' flag is True)")

	.def_property("mean_flux", 
		      &von_mises_profile::get_mean_flux, 
		      &von_mises_profile::set_mean_flux,
		      "Mean flux of the pulse (settable property).\n\n"
		      "Note: we define the mean flux with no detrending correction applied (even if the 'detrend' flag is True)")

	.def("get_single_pulse_signal_to_noise", 
	     &von_mises_profile::get_single_pulse_signal_to_noise,
	     "dt_sample"_a, "pulse_freq"_a, "sample_rms"_a = 1.0,

	     "get_single_pulse_signal_to_noise(dt_sample, pulse_freq, sample_rms=1.0) -> float\n"
	     "\n"
	     "Returns the SNR of a single pulse, accounting for finite time resolution (and detrending,\n"
	     "if detrend=True was specified at construction).\n"
	     "\n"
	     "Strictly speaking, the return value is an approximation to the true SNR, which may slightly depend\n"
	     "on the exact arrival time of the pulse.\n"
	     "\n"
	     "  - The 'dt_sample' argument is the length of each time sample.\n\n"
	     "  - The 'pulse_freq' argument is the pulse frequency.\n\n"
	     "  - The 'sample_rms' argument is the RMS noise fluctuation in each time sample.")

	.def("get_multi_pulse_signal_to_noise", 
	     &von_mises_profile::get_multi_pulse_signal_to_noise,
	     "total_time"_a, "dt_sample"_a, "pulse_freq"_a, "sample_rms"_a = 1.0,

	     "get_multi_pulse_signal_to_noise(total_time, dt_sample, pulse_freq, sample_rms=1.0) -> float\n"
	     "\n"
	     "Returns the total SNR of a train of pulses, accounting for finite time resolution (and detrending,\n"
	     "if detrend=True was specified at construction).\n"
	     "\n"
	     "Strictly speaking, this is an approximation to the true SNR, which may slightly depend\n"
	     "on the exact arrival times of the pulses.\n"
	     "\n"
	     "  - The 'total_time' argument is the total duration of the pulse train.\n\n"
	     "  - The 'dt_sample' argument is the length of each time sample.\n\n"
	     "  - The 'pulse_freq' argument is the pulse frequency.\n\n"
	     "  - The 'sample_rms' argument is the RMS noise fluctuation in each time sample.")

	.def("set_single_pulse_signal_to_noise", 
	     &von_mises_profile::set_single_pulse_signal_to_noise,
	     "snr"_a, "dt_sample"_a, "pulse_freq"_a, "sample_rms"_a = 1.0,

	     "set_single_pulse_signal_to_noise(snr, dt_sample, pulse_freq, sample_rms=1.0)\n"
	     "\n"
	     "Changes the pulse normalization to give the specified single-pulse signal-to-noise ratio.\n"
	     "For more information and description of the arguments, see the get_single_pulse_signal_to_noise() docstring.")

	.def("set_multi_pulse_signal_to_noise", 
	     &von_mises_profile::set_multi_pulse_signal_to_noise,
	     "snr"_a, "total_time"_a, "dt_sample"_a, "pulse_freq"_a, "sample_rms"_a = 1.0,

	     "set_multi_pulse_signal_to_noise(snr, total_time, dt_sample, pulse_freq, sample_rms=1.0)\n"
	     "\n"
	     "Changes the pulse normalization to give the specified multi-pulse signal-to-noise ratio.\n"
	     "For more information and description of the arguments, see the get_multi_pulse_signal_to_noise() docstring.")

	.def("get_profile_fft", get_profile_fft, "nout"_a = 0,
	     "get_profile_fft(nout=0) -> (1D array)\n"
	     "\n"
	     "Returns the Fourier transform of the profile:\n\n"
	     "  rho_m = int_0^1 dphi rho(phi) e^{i m phi}.\n"
	     "\n"
	     "Note that rho_m is real, and rho_m = rho_{-m}, since the von Mises profile is symmetric.\n\n"
	     "The DC mode rho_0 will equal 'mean_flux' if detrend=False, or 0 if detrend=True.\n"
	     "\n"
	     "The return value is a 1D array of length 'nout'.  If nout=0 (the default), then it defaults to\n"
	     "(internal_nphi/2+1), the number of Fourier coefficients which are computed internally.  (If 'nout'\n"
	     "is larger than this, then the returned array is zero-padded.\n"
	     "\n"
	     "The returned FFT does not include a time-sample window function.  You may want to multiply the\n"
	     "output by j_0(pi*dphi/2), where dphi = (pulse_freq * dt_sample).")
    ;
}


}  // namespace simpulse_pybind11
