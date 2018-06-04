#include "../include/simpulse.hpp"
#include "../include/simpulse/internals.hpp"

using namespace std;


namespace simpulse {
#if 0
};  // pacify emacs c-mode!
#endif


struct constant_acceleration_phase_model : public phase_model {
    double phi0 = 0.0;
    double omega0 = 0.0;
    double omega_dot = 0.0;
    double t0 = 0.0;

    constant_acceleration_phase_model(double phi0_, double omega0_, double omega_dot_, double t0_)
        : phi0(phi0_), omega0(omega0_), omega_dot(omega_dot_), t0(t0_)
    {
	if (omega0 <= 0.0)
            throw std::runtime_error("simpulse::constant_acceleration_phase_model constructor: expected omega0 > 0");
    }

    virtual double eval(double t, int nderivs) const override
    {
	t -= t0;

	if (nderivs == 0)
	    return phi0 + omega0*t + omega_dot*t*t/2.;
	if (nderivs == 1)
	    return omega0 + omega_dot*t;
	if (nderivs == 2)
	    return omega_dot;
	if (nderivs < 0)
	    throw runtime_error("simpulse::phase_model::eval() was called with negative 'nderivs'");

	return 0.0;
    }

    virtual double eval_mean_phi(double t1, double t2) const override
    {
	double tmid = (t1+t2)/2. - t0;
	double dt = t2-t1;
	return phi0 + omega0*tmid + omega_dot*tmid*tmid/2. + omega_dot*dt*dt/24.;
    }

    virtual void eval_phi(int nt, double *phi_out, const double *t_in) const override
    {
	if (nt <= 0)
	    throw runtime_error("simpulse::phase_model::eval_phi(): expected nt > 0");

	for (int it = 0; it < nt; it++) {
	    double t = t_in[it] - t0;
	    phi_out[it] = phi0 + omega0*t + 0.5*omega_dot*t*t;
	}
    }
};


// static member function
shared_ptr<phase_model> phase_model::make_constant_acceleration(double phi0, double omega0, double omega_dot, double t0)
{
    return make_shared<constant_acceleration_phase_model> (phi0, omega0, omega_dot, t0);
}


}  // namespace simpulse
