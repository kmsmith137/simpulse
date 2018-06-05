#ifndef _SIMPULSE_INLINES_HPP
#define _SIMPULSE_INLINES_HPP

#if (__cplusplus < 201103) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
#error "This source file needs to be compiled with C++0x support (g++ -std=c++0x)"
#endif

#include <cmath>

namespace simpulse {
#if 0
}; // pacify emacs c-mode
#endif


// Returns dispersion delay in seconds
// The 'DM' parameter is the dispersion measure in its standard units (pc cm^{-3})
inline double dispersion_delay(double dm, double freq_MHz)
{
    return 4.148806e3 * dm / (freq_MHz * freq_MHz);
}

// Returns scattering time in seconds
// We define the 'SM' to be the scattering time in _milliseconds_ (not seconds) at 1 GHz
inline double scattering_time(double sm, double freq_MHz)
{
    return 1.0e-3 * sm / pow(freq_MHz/1000.0, 4.4);
}


}  // namespace simpulse

#endif // _SIMPULSE_INLINES_HPP
