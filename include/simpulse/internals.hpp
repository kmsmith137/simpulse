#ifndef _SIMPULSE_INTERNALS_HPP
#define _SIMPULSE_INTERNALS_HPP

#if (__cplusplus < 201103) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
#error "This source file needs to be compiled with C++0x support (g++ -std=c++0x)"
#endif

#include <cmath>
#include <cstring>
#include <complex>
#include <sstream>
#include <stdexcept>
#include <fftw3.h>

#include "simpulse.hpp"


// Branch predictor hint
#ifndef _unlikely
#define _unlikely(cond)  (__builtin_expect(cond,0))
#endif

// simpulse_assert(): like assert, but throws an exception in order to work smoothly with python.
#define simpulse_assert(cond) simpulse_assert2(cond, __LINE__)

#define simpulse_assert2(cond,line) \
    do { \
        if (_unlikely(!(cond))) { \
	    const char *msg = "simpulse: assertion '" __STRING(cond) "' failed (" __FILE__ ":" __STRING(line) ")\n"; \
	    throw std::runtime_error(msg); \
	} \
    } while (0)


namespace simpulse {
#if 0
}; // pacify emacs c-mode
#endif


// Catering to old versions of gcc
// FIXME was this really necessary?
template<typename T> inline std::string xto_string(const T &x)
{
    std::stringstream ss;
    ss << x;
    return ss.str();
}

template<typename T> inline T *checked_fftw_malloc(size_t nelts)
{
    size_t nbytes = nelts * sizeof(T);

    void *ret = fftw_malloc(nbytes);
    if (!ret)
	throw std::runtime_error("fftw_malloc couldn't allocate " + xto_string(nbytes) + " bytes");

    memset(ret, 0, nbytes);
    return reinterpret_cast<T *> (ret);
}

inline double square(double x)
{
    return x*x;
}

// Returns j0(x) = sin(x)/x
inline double bessj0(double x)
{
    return (x*x > 1.0e-100) ? (sin(x)/x) : 1.0;
}

// These implementations of round_down() and round_up() are guaranteed correct for negative arguments
inline int round_down(double x)
{
    int i = (int)x;
    return (i <= x) ? i : (i-1);
}

inline int round_up(double x)
{
    int i = (int)x;
    return (i >= x) ? i : (i+1);
}


}  // namespace simpulse

#endif // _SIMPULSE_INTERNALS_HPP
