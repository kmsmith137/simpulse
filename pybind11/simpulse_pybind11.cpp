#include "../include/simpulse.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
using namespace pybind11::literals;
using namespace simpulse;
using namespace std;


PYBIND11_MODULE(simpulse_pybind11, m)
{
    m.doc() = "simpulse: C++/python library for simulating pulses in radio astronomy";

    m.def("dispersion_delay", &dispersion_delay, "dm"_a, "freq_MHz"_a,
	  "dispersion_delay(dm, freq_MHz) -> dispersion delay in seconds");

    m.def("scattering_time", &scattering_time, "sm"_a, "freq_MHz"_a,
	  "scattering_time(sm, freq_MHz) -> scattering time in seconds");
}
