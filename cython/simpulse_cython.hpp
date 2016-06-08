#include "simpulse.hpp"

namespace simpulse {
#if 0
}; // pacify emacs c-mode
#endif

inline void _add_to_timestream_float(simpulse::single_pulse *sp, float *out, double out_t0, double out_t1, int out_nt, int stride)
{
    sp->add_to_timestream(out, out_t0, out_t1, out_nt, stride);
}

inline void _add_to_timestream_double(simpulse::single_pulse *sp, double *out, double out_t0, double out_t1, int out_nt, int stride)
{
    sp->add_to_timestream(out, out_t0, out_t1, out_nt, stride);
}

}  // namespace simpulse
