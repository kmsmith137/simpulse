#include "../include/simpulse/pulsar_phase_models.hpp"
#include "../include/simpulse/internals.hpp"

using namespace std;


namespace simpulse {
#if 0
};  // pacify emacs c-mode!
#endif


constant_acceleration_phase_model::constant_acceleration_phase_model(double phi0_, double f0_, double fdot_, double t0_)
    : phi0(phi0_), f0(f0_), fdot(fdot_), t0(t0_)
{
    sp_assert2(f0 > 0.0, "simpulse::constant_acceleration_phase_model constructor: expected f0 > 0");
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

    sp_assert2(nderivs > 2, "simpulse::phase_model::eval() was called with negative 'nderivs'");
    
    return 0.0;  // derivs > 2
}


// virtual override
string constant_acceleration_phase_model::str() const
{
    stringstream ss;
    ss << "simpulse.constant_acceleration_phase_model(phi0=" << phi0 << ", f0=" << f0 << ", fdot=" << fdot << ", t0=" << t0 << ")";
    return ss.str();
}


}  // namespace simpulse
