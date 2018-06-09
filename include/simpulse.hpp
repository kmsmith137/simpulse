#ifndef _SIMPULSE_HPP
#define _SIMPULSE_HPP

#if (__cplusplus < 201103) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
#error "This source file needs to be compiled with C++0x support (g++ -std=c++0x)"
#endif

// Currently only contains two functions, dispersion_delay() and scattering_time()
#include "simpulse/inlines.hpp"

// single_pulse: this class represents one dispersed, scattered pulse, in a fixed frequency channelization.
#include "simpulse/single_pulse.hpp"

// Pulsar phase models + 'class von_mises_profile'.
// Note that in the pulsar case, we only simulate 1D timestreams,
// whereas in the single_pulse case, we simulate 2D (frequency,time) arrays).
#include "simpulse/pulsar_phase_models.hpp"
#include "simpulse/pulsar_profiles.hpp"

// Not included by default
// #include "simpulse/internals.hpp"

#endif // _SIMPULSE_HPP
