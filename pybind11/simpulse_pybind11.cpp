#include "simpulse_pybind11.hpp"
#include "../include/simpulse/inlines.hpp"   // dispersion_delay(), scattering_time()

#include <pybind11/numpy.h>

namespace py = pybind11;
using namespace pybind11::literals;
using namespace simpulse;
using namespace std;


static void wrap_inlines(py::module &m)
{
    py::options options;
    options.disable_function_signatures();

    m.def("dispersion_delay", &dispersion_delay, "dm"_a, "freq_MHz"_a,
	  "dispersion_delay(dm, freq_MHz) -> float\n"
	  "\n"
	  "Returns the dispersion delay in seconds, given:\n"
	  "\n"
	  "   - The dispersion measure (in its usual units, pc cm^{-3}).\n\n"
	  "   - The frequency (in MHz).");

    m.def("scattering_time", &scattering_time, "sm"_a, "freq_MHz"_a,
	  "scattering_time(sm, freq_MHz) -> float\n"
	  "\n"
	  "Returns the scattering timescale in seconds (not milliseconds!), given\n"
	  "\n"
	  "   - The scattering measure 'sm', which we define to be the scattering time\n"
	  "     in milliseconds (not seconds!) at 1 GHz.\n"
	  "\n"
	  "   - The frequency (in MHz).");
}


PYBIND11_MODULE(simpulse_pybind11, m)
{
    m.doc() = "simpulse: C++/python library for simulating pulses in radio astronomy";

    wrap_inlines(m);

    simpulse_pybind11::wrap_single_pulse(m);
    simpulse_pybind11::wrap_phase_model_base(m);
    simpulse_pybind11::wrap_constant_acceleration_phase_model(m);
    simpulse_pybind11::wrap_sinusoidal_phase_model(m);
    simpulse_pybind11::wrap_keplerian_binary_phase_model(m);
    simpulse_pybind11::wrap_von_mises_profile(m);
}
