#include "simpulse_pybind11.hpp"
#include "../include/simpulse/pulsar_phase_models.hpp"

#include <pybind11/numpy.h>

namespace py = pybind11;
using namespace pybind11::literals;
using namespace simpulse;
using namespace std;

namespace simpulse_pybind11 {
#if 0
}  // pacify emacs c-mode
#endif



void wrap_phase_model_base(py::module &m)
{
    auto eval_phi_sequence = [](phase_model_base &self, double t0, double t1, ssize_t nsamples, int nderivs)
    {
	if (nsamples <= 0)
	    throw runtime_error("simpulse::phase_model::eval_phi_sequence(): expected nsamples > 0");

	py::array_t<double> ret{ size_t(nsamples) };

	self.eval_phi_sequence(t0, t1, nsamples, ret.mutable_data(), nderivs);
	return ret;
    };

    py::class_<phase_model_base>(m, "phase_model_base")
	.def("eval_phi", &phase_model_base::eval_phi, "t"_a, "nderivs"_a = 0)
	.def("eval_phi_sequence", eval_phi_sequence, "t0"_a, "t1"_a, "nsamples"_a, "nderivs"_a = 0)
    ;
}


void wrap_constant_acceleration_phase_model(py::module &m)
{
    py::class_<constant_acceleration_phase_model, phase_model_base>(m, "constant_acceleration_phase_model")
	.def(py::init<double,double,double,double>(), "phi0"_a, "f0"_a, "fdot"_a, "t0"_a)

	.def_readonly("phi0", &constant_acceleration_phase_model::phi0)
	.def_readonly("f0", &constant_acceleration_phase_model::f0)
	.def_readonly("fdot", &constant_acceleration_phase_model::fdot)
	.def_readonly("t0", &constant_acceleration_phase_model::t0)
    ;
}


}  // namespace simpulse_pybind11
