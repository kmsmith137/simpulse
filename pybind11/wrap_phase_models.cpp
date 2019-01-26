#include "simpulse_pybind11.hpp"
#include "simpulse_pybind11_array_helpers.hpp"

#include "../include/simpulse/internals.hpp"  // sp_assert()
#include "../include/simpulse/pulsar_phase_models.hpp"


namespace py = pybind11;
using namespace pybind11::literals;
using namespace simpulse;
using namespace std;

namespace simpulse_pybind11 {
#if 0
}  // pacify emacs c-mode
#endif


struct phase_model_trampoline : public phase_model_base 
{
    virtual double eval_phi(double t, int nderivs) const override
    {
	// FIXME replace this by something which prints a more helpful error message if function is not defined on python side.
	PYBIND11_OVERLOAD_PURE(double, phase_model_base, eval_phi, t, nderivs);
    }


    virtual void eval_phi_sequence(double t0, double t1, ssize_t nsamples, double *phi_out, int nderivs) const override
    {
	py::gil_scoped_acquire gil;
	py::function overload = pybind11::get_overload(this, "eval_phi_sequence");

	if (!overload)
	    return phase_model_base::eval_phi_sequence(t0, t1, nsamples, phi_out, nderivs);

	auto o = overload(t0, t1, nsamples, nderivs);
	auto a = py::cast<in_carray<double>> (o);

	if (!is_contiguous(a))
	    throw runtime_error("simpulse internal error: is_contiguous(in_carray<T>) returned false, this should never happen");

	if (!shape_equals(a, 1, &nsamples)) {
	    stringstream ss;
	    ss << "simpulse: eval_phi_sequence() returned shape-" << shape_string(a) << " array, expected shape (" << nsamples << ",)";
	    throw runtime_error(ss.str());
	}

	memmove(phi_out, a.data(), nsamples * sizeof(double));
    }


    virtual std::string str() const override
    {
	py::gil_scoped_acquire gil;
	py::function overload = pybind11::get_overload(this, "__repr__");

	if (!overload) {
	    string class_name = py::cast(this).get_type().attr("__name__").cast<string>();
	    return "<" + class_name + " instance>";
	}
	
	auto o = overload();
	return py::cast<string> (o);
    }
};


void wrap_phase_model_base(py::module &m)
{
    py::options options;
    options.disable_function_signatures();

    const char *doc =
	"This virtual base class represents a pulsar phase model Phi(t).  By definition, this is a function which takes\n"
	"integer values at times when pulses are observed.  The specific form of Phi(t) is determined by the subclass.\n"
	"\n"
	"Currently, only one subclass is implemented: constant_acceleration_phase_model.\n"
	"\n"
	"It is possible to define subclasses from python.  A python subclass must define the eval_phi() method,\n"
	"and may optionally define eval_phi_sequence().  See :ref:`Example 3: python phase model` for an example.";


    auto eval_phi_sequence = [](phase_model_base &self, double t0, double t1, ssize_t nsamples, int nderivs)
    {
	sp_assert2(nsamples > 0, "simpulse::phase_model::eval_phi_sequence(): expected nsamples > 0");
	py::array_t<double> ret{ size_t(nsamples) };

	self.eval_phi_sequence(t0, t1, nsamples, ret.mutable_data(), nderivs);
	return ret;
    };


    py::class_<phase_model_base, phase_model_trampoline>(m, "phase_model_base", doc)
	.def(py::init<>())

	.def("eval_phi", &phase_model_base::eval_phi, "t"_a, "nderivs"_a = 0,
	     "eval_phi(t, nderivs=0) -> float\n\n"
	     "Evaluates the phase model Phi(t) at a single time 't'.\n"
	     "\n"
	     "The ``nderivs`` argument can be used to compute the n-th derivative (d^n phi) / dt^n.\n"
	     "Currently we don't actually use 'nderivs' anywhere, and I may remove it from the\n"
	     "simpulse API eventually!\n")

	.def("eval_phi_sequence", eval_phi_sequence, "t0"_a, "t1"_a, "nsamples"_a, "nderivs"_a = 0,
	     "eval_phi_sequence(t0, t1, nsamples, nderivs=0) -> (1D array)\n\n"
	     "Evaluates the phase model at a sequence of equally spaced time samples.\n"
	     "The 't0' and 't1' args are the starting/ending times of the sampled region.\n\n"
	     "The return value is an array of length 'nsamples', whose i-th entry is the phase model Phi(t)\n"
	     "evaluated at t = ((nsamples-i-1)*t0 + i*t1) / (nsamples-1).\n\n"
	     "Can be called with nderivs=0 to get Phi(t), or nderivs > 0 to get derivatives (d/dt)^n Phi.")

	.def("__repr__", &phase_model_base::str)
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
	"    t0 = a reference time (can be zero, but a nonzero value may be convenient)\n"
	"\n"
	"WARNING: constant_acceleration_phase_models are pickleable, but *only* if you\n"
	"use the cPickle (not pickle) module, and specify pickle protocol version 2, e.g.::\n"
	"\n"
	"    x = simpulse.constant_acceleration_phase_model(...)\n"
	"    s = cPickle.dumps(x, 2)    # note protocol version 2 here!\n"
	"    print cPickle.loads(x)\n"
	"\n"
	"Variants of this will result in either exceptions or mystery segfaults!!\n";

    auto __getstate__ = [](const constant_acceleration_phase_model &p)
    {
	return py::make_tuple(p.phi0, p.f0, p.fdot, p.t0);
    };

    auto __setstate__ = [](py::tuple t)
    {
	if (t.size() != 4)
	    throw runtime_error("Internal error in simpule::constant_acceleration_phase_model::__setstate__(): len(t) != 4");

	double phi0 = t[0].cast<double> ();
	double f0 = t[1].cast<double> ();
	double fdot = t[2].cast<double> ();
	double t0 = t[3].cast<double> ();
	
	return constant_acceleration_phase_model(phi0, f0, fdot, t0);
    };

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

	.def(py::pickle(__getstate__, __setstate__))
    ;
}


void wrap_keplerian_binary_phase_model(py::module &m)
{
    py::options options;
    options.disable_function_signatures();

    // This string becomes the python docstring, and also goes into the sphinx documentation.
    // If you change it and do 'make install', the python docstring will be updated.
    //
    // However, the sphinx documentation won't be updated until you run the 'run-sphinx.py'
    // script in the simpulse_docs repo.  (This updates the "local" sphinx documentation,
    // to update the web-published version, you need to do 'git push' in the simpulse_docs repo.)

    const char *doc =
	"This class represents a binary pulsar with relativistic effects neglected.\n"
	"(Subclass of phase_model_base.)\n"
	"\n"
	"Constructor syntax::\n"
	"\n"
	"    pm = simpulse.keplerian_binary_phase_model(e, a, Porb, nx, ny, P, t0, phi0)\n"
	"\n"
	"where:\n\n"
	"    e = eccentricity\n\n"
	"    a = semimajor axis\n\n"
	"    Porb = orbital period\n\n"
	"    nx, ny = unit vector in the direction of Earth\n\n"
	"    P = pulse period\n\n"
	"    t0 = time delay parameter between Earth and binary center of mass (mod Porb)\n\n"
	"    phi0 = initial phase (mod 1)\n"
	"\n"
	"WARNING: keplerian_binary_phase_models are pickleable, but *only* if you\n"
	"use the cPickle (not pickle) module, and specify pickle protocol version 2, e.g.::\n"
	"\n"
	"    x = simpulse.keplerian_binary_phase_model(...)\n"
	"    s = cPickle.dumps(x, 2)    # note protocol version 2 here!\n"
	"    print cPickle.loads(x)\n"
	"\n"
	"Variants of this will result in either exceptions or mystery segfaults!!\n";
    

    auto __getstate__ = [](const keplerian_binary_phase_model &p)
    {
	return py::make_tuple(p.e, p.a, p.Porb, p.nx, p.ny, p.P, p.t0, p.phi0);
    };

    auto __setstate__ = [](py::tuple t)
    {
	if (t.size() != 8)
	    throw runtime_error("Internal error in simpule::keplerian_binary_phase_model::__setstate__(): len(t) != 8");

	double e = t[0].cast<double> ();
	double a = t[1].cast<double> ();
	double Porb = t[2].cast<double> ();
	double nx = t[3].cast<double> ();
	double ny = t[4].cast<double> ();
	double P = t[5].cast<double> ();
	double t0 = t[6].cast<double> ();
	double phi0 = t[7].cast<double> ();

	return keplerian_binary_phase_model(e, a, Porb, nx, ny, P, t0, phi0);
    };


    
    py::class_<keplerian_binary_phase_model, phase_model_base>(m, "keplerian_binary_phase_model", doc)
	// The following boilerplate exports the C++ constructor to python.
	// Note that you need to declare the constructor argument types explicitly (double, double, ...)
	// and provide names for the arguments.  
	// Warning: if something is wrong here (for example, if the number of arguments doesn't match the 
	// C++ constructor) you'll get an incomprehensible compiler error!)
	.def(py::init<double,double,double,double,double,double,double,double>(), 
	     "e"_a, "a"_a, "Porb"_a, "nx"_a, "ny"_a, "P"_a, "t0"_a, "phi0"_a)

	// Make const public members visible from python.  (The last argument is the docstring.)
	.def_readonly("e", &keplerian_binary_phase_model::e, "Eccentricity")
	.def_readonly("a", &keplerian_binary_phase_model::a, "Semimajor axis")
	.def_readonly("b", &keplerian_binary_phase_model::b, "Semiminor axis")
	.def_readonly("Porb", &keplerian_binary_phase_model::Porb, "Orbital period")
	.def_readonly("nx", &keplerian_binary_phase_model::nx, "X-component of unit vector in direction of Earth")
	.def_readonly("ny", &keplerian_binary_phase_model::ny, "Y-component of unit vector in direction of Earth")
	.def_readonly("P", &keplerian_binary_phase_model::P, "Pulse period")
	.def_readonly("t0", &keplerian_binary_phase_model::t0, "Time delay parameter between Earth and binary center of mass (mod Porb)")
	.def_readonly("phi0", &keplerian_binary_phase_model::phi0, "Initial phase (mod 1)")

	.def(py::pickle(__getstate__, __setstate__))
    ;
}


}  // namespace simpulse_pybind11
