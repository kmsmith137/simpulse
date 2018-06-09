#ifndef _SIMPULSE_PULSAR_PHASE_MODELS_HPP
#define _SIMPULSE_PULSAR_PHASE_MODELS_HPP

#if (__cplusplus < 201103) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
#error "This source file needs to be compiled with C++0x support (g++ -std=c++0x)"
#endif

#include <vector>
#include <memory>

namespace simpulse {
#if 0
}; // pacify emacs c-mode
#endif


struct phase_model_base 
{
    // Evaluates the phase model at a single time 't'.
    // Can be called with nderivs=0 to get Phi(t), or nderivs > 0 to get derivatives Phi^{(n)}(t).

    virtual double eval_phi(double t, int nderivs=0) const = 0;


    // Evaluates the phase model at a sequence of equally spaced time samples:
    // The 't0' and 't1' args are the starting/ending times of the sampled region.
    //
    // eval_phi_sequence() is not pure virtual: there is a default implementation which loops over eval_phi().
    // Subclasses may optionally override eval_sequence() as an optimization to improve performace.  
    // (FIXME: I'm planning to time this, to see how much it actually matters!)

    virtual void eval_phi_sequence(double t0, double t1, ssize_t nsamples, double *phi_out, int nderivs=0) const;
};


struct constant_acceleration_phase_model : public phase_model_base 
{
    const double phi0;
    const double f0;
    const double fdot;
    const double t0;

    constant_acceleration_phase_model(double phi0, double f0, double fdot, double t0);

    virtual double eval_phi(double t, int nderivs) const override;
};
    


}  // namespace simpulse

#endif // _SIMPULSE_PULSAR_PHASE_MODELS_HPP
