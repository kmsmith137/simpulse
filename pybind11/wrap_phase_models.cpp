#include "simpulse_pybind11.hpp"
#include "../include/simpulse/internals.hpp"  // sp_assert()
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
    py::options options;
    options.disable_function_signatures();

    const char *doc =
	"This virtual base class represents a pulsar phase model Phi(t).  By definition, this is a function which takes\n"
	"integer values at times when pulses are observed.  The specific form of Phi(t) is determined by the subclass.\n"
	"\n"
	"Currently, only one subclass of phase_model_base is implemented: constant_acceleration_phase_model (see below).\n";


    auto eval_phi_sequence = [](phase_model_base &self, double t0, double t1, ssize_t nsamples, int nderivs)
    {
	sp_assert2(nsamples > 0, "simpulse::phase_model::eval_phi_sequence(): expected nsamples > 0");
	py::array_t<double> ret{ size_t(nsamples) };

	self.eval_phi_sequence(t0, t1, nsamples, ret.mutable_data(), nderivs);
	return ret;
    };


    py::class_<phase_model_base>(m, "phase_model_base", doc)
	.def("eval_phi", &phase_model_base::eval_phi, "t"_a, "nderivs"_a = 0,
	     "eval_phi(t, nderivs=0) -> float\n\n"
	     "Evaluates the phase model Phi(t) at a single time 't'.\n"
	     "Can be called with nderivs=0 to get Phi(t), or nderivs > 0 to get derivatives (d/dt)^n Phi.")

	.def("eval_phi_sequence", eval_phi_sequence, "t0"_a, "t1"_a, "nsamples"_a, "nderivs"_a = 0,
	     "eval_phi_sequence(t0, t1, nsamples, nderivs=0) -> (1D array)\n\n"
	     "Evaluates the phase model at a sequence of equally spaced time samples.\n"
	     "The 't0' and 't1' args are the starting/ending times of the sampled region.\n\n"
	     "The return value is an array of length 'nsamples', whose i-th entry is the phase model Phi(t)\n"
	     "evaluated at t = ((nsamples-i-1)*t0 + i*t1) / (nsamples-1).\n\n"
	     "Can be called with nderivs=0 to get Phi(t), or nderivs > 0 to get derivatives (d/dt)^n Phi.")
    ;
}


void wrap_constant_acceleration_phase_model(py::module &m)
{
    py::options options;
    options.disable_function_signatures();

    const char *doc =
	"This class represents a pulsar whose acceleration (i.e. frequency derivative df/dt)\n"
	"is constant in time.  (Subclass of phase_model_base.)\n"
	"\n"
	"Constructor syntax::\n"
	"\n"
	"    pm = simpulse.constant_acceleration_phase_model(phi0, f0, fdot, t0)\n"
	"\n"
	"where:\n\n"
	"    phi0 = value of the phase at reference time 't0'\n\n"
	"    f0 = value of the frequency f=dphi/dt at reference time 't0'\n\n"
	"    fdot = frequency derivative df/dt (note that this is independent of time)\n\n"
	"    t0 = a reference time (can be zero, but a nonzero value may be convenient)";
    

    py::class_<constant_acceleration_phase_model, phase_model_base>(m, "constant_acceleration_phase_model", doc)
	.def(py::init<double,double,double,double>(), "phi0"_a, "f0"_a, "fdot"_a, "t0"_a)

	.def_readonly("phi0", &constant_acceleration_phase_model::phi0, 
		      "Value of the phase at reference time 't0' (same as constructor argument with same name)")

	.def_readonly("f0", &constant_acceleration_phase_model::f0, 
		      "Value of the frequency at reference time 't0' (same as constructor argument with same name)")

	.def_readonly("fdot", &constant_acceleration_phase_model::fdot,
		      "Frequency derivative df/dt (same as constructor argument with same name)")

	.def_readonly("t0", &constant_acceleration_phase_model::t0,
		      "Reference time where 'phi0' and 'f0' are defined (same as constructor argument with same name)")
    ;
}


}  // namespace simpulse_pybind11
