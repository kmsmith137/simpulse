#ifndef _SIMPULSE_PYBIND11_HPP
#define _SIMPULSE_PYBIND11_HPP

#include <iostream>
#include <pybind11/pybind11.h>


namespace simpulse_pybind11 {
#if 0
}  // pacify emacs c-mode
#endif


// wrap_single_pulse.cpp
void wrap_single_pulse(pybind11::module &m);

// wrap_phase_models.cpp
void wrap_phase_model_base(pybind11::module &m);
void wrap_constant_acceleration_phase_model(pybind11::module &m);
void wrap_keplerian_binary_phase_model(pybind11::module &m);

// wrap_von_mises_profile.cpp
void wrap_von_mises_profile(pybind11::module &m);


}  // namespace simpulse_pybind11

#endif // _SIMPULSE_PYBIND11_HPP
