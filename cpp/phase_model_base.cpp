#include "../include/simpulse/pulsar_phase_models.hpp"
#include "../include/simpulse/internals.hpp"

using namespace std;


namespace simpulse {
#if 0
};  // pacify emacs c-mode!
#endif


void phase_model_base::eval_phi_sequence(double t0, double t1, ssize_t nsamples, double *phi_out, int nderivs) const
{
    sp_assert2(nsamples > 0, "phase_model_base::eval_phi_sequence(): expected nsamples > 0");
    sp_assert2(phi_out, "phase_model_base::eval_phi_sequence(): 'phi_out' is a null pointer");

    if (nsamples == 1)
	phi_out[0] = this->eval_phi(0.5 * (t0+t1));

    // nsamples > 1
    for (ssize_t i = 0; i < nsamples; i++) {
	double t = ((nsamples-1-i)*t0 + (i)*t1) / double(nsamples-1);
	phi_out[i] = this->eval_phi(t, nderivs);
    }
}


}  // namespace simpulse
