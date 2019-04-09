#ifndef _SIMPULSE_INLINES_HPP
#define _SIMPULSE_INLINES_HPP

#if (__cplusplus < 201103) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
#error "This source file needs to be compiled with C++0x support (g++ -std=c++0x)"
#endif

#include <cmath>
#include <mutex>

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


// Returns scattering time in seconds (not milliseconds!)
// We define the 'SM' to be the scattering time in _milliseconds_ (not seconds!) at 1 GHz
inline double scattering_time(double sm, double freq_MHz)
{
    return 1.0e-3 * sm / pow(freq_MHz/1000.0, 4.4);
}


// Global lock protecting the FFTW planner.  This should really be in libfftw, not simpulse!
// However:
//
//   - A global planning lock was only added to FFTW in 3.3.5
//   - It didn't really work until 3.3.6, so users need to check the fftw_version() and report an error if ==3.3.5.
//   - FFTW needs to be configured with --enable-threads, and needs to be linked with -lfftw3_threads
//   - Currently (April 2019), most machines "in the wild" don't have a version of FFTW which satisfies these constraints.
//
// As a result, I decided that using FFTW's lock would be inconvenient for users, and I just put a lock in simpulse
// instead.  This works as long as simpulse is the "bottom-level" library which uses libfftw (e.g. libproton also uses
// libfftw, but it can "see" and use the simpulse lock).  Eventually I'll phase out simpulse::fftw_global_planning_lock,
// add the requirement that FFTW >= 3.3.6, and call fftw_make_planner_thread_safe().

extern std::mutex fftw_global_planning_lock;


}  // namespace simpulse

#endif // _SIMPULSE_INLINES_HPP
