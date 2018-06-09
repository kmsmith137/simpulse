#include "../include/simpulse.hpp"
#include "../include/simpulse/internals.hpp"

using namespace std;


namespace simpulse {
#if 0
};  // pacify emacs c-mode!
#endif


constant_acceleration_phase_model::constant_acceleration_phase_model(double phi0_, double f0_, double fdot_, double t0_)
    : phi0(phi0_), f0(f0_), fdot(fdot_), t0(t0_)
{
    if (f0 <= 0.0)
	throw std::runtime_error("simpulse::constant_acceleration_phase_model constructor: expected f0 > 0");
}


// virtual override
double constant_acceleration_phase_model::eval_phi(double t, int nderivs) const
{
    t -= t0;
    
    if (nderivs == 0)
	return phi0 + f0*t + fdot*t*t/2.;
    if (nderivs == 1)
	return f0 + fdot*t;
    if (nderivs == 2)
	return fdot;

    if (_unlikely(nderivs < 0))
	throw runtime_error("simpulse::phase_model::eval() was called with negative 'nderivs'");
    
    return 0.0;  // derivs > 2
}


}  // namespace simpulse
