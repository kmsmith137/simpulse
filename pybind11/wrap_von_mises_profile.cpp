#include "simpulse_pybind11.hpp"
#include "../include/simpulse/pulsar_phase_models.hpp"
#include "../include/simpulse/pulsar_profiles.hpp"
#include "../include/simpulse/internals.hpp"   // simpulse_assert()

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
    auto eval_integrated_samples = [](von_mises_profile &self, double t0, double t1, ssize_t nt, const phase_model_base &pm, double amplitude)
    {
	if (nt <= 0)
	    throw runtime_error("simpulse::von_mises_profile::eval_integrated_sample(): expected nt > 0");

	py::array_t<double> ret{ size_t(nt) };

	self.eval_integrated_samples(ret.mutable_data(), t0, t1, nt, pm, amplitude);
	return ret;
    };


    auto get_profile_fft = [](von_mises_profile &self, int nout)
    {
	if (nout < 0)
	    throw runtime_error("simpulse::von_mises_profile::get_profile_fft(): expected nout >= 0");

	if (nout == 0)
	    nout = self.internal_nphi/2 + 1;

	py::array_t<double> ret{ size_t(nout) };

	self.get_profile_fft(ret.mutable_data(), nout);
	return ret;
    };


    py::class_<von_mises_profile>(m, "von_mises_profile")
	.def(py::init<double,bool,int>(), "duty_cycle"_a, "detrend"_a, "min_internal_nphi"_a = 0,
	     "We define the duty cycle to be (pulse FWHM) / (pulse period).\n"
	     "If min_nphi=0, then a reasonable default will be chosen (recommended).")

	// Read-only properties (constant after construction)
	.def_readonly("duty_cycle", &von_mises_profile::duty_cycle, "We define the duty cycle to be (pulse FWHM) / (pulse period).")
	.def_readonly("detrend", &von_mises_profile::detrend)
	.def_readonly("internal_nphi", &von_mises_profile::internal_nphi)
	.def_property_readonly("mean_flux", &von_mises_profile::get_mean_flux)

	.def("eval_integrated_samples", eval_integrated_samples, "t0"_a, "t1"_a, "nt"_a, "phase_model"_a, "amplitude"_a = 1.0)
	.def("eval_integrated_sample_slow", &von_mises_profile::eval_integrated_sample_slow, "phi0"_a, "phi1"_a, "amplitude"_a = 1.0)
	.def("point_eval", &von_mises_profile::point_eval, "phi"_a, "amplitude"_a = 1.0)

	.def("get_single_pulse_signal_to_noise", &von_mises_profile::get_single_pulse_signal_to_noise,
	     "dt_sample"_a, "pulse_freq"_a, "sample_rms"_a = 1.0)

	.def("get_multi_pulse_signal_to_noise", &von_mises_profile::get_multi_pulse_signal_to_noise,
	     "total_time"_a, "dt_sample"_a, "pulse_freq"_a, "sample_rms"_a = 1.0)

	.def("get_profile_fft", get_profile_fft, "nout"_a = 0)
    ;
}


}  // namespace simpulse_pybind11
