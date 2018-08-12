// When implementing 'class keplerian_binary_phase_model', you may find it useful to look at
// constant_acceleration_phase_model.cpp, for a simple example of a phase_model class in C++.

#include "../include/simpulse/pulsar_phase_models.hpp"
#include "../include/simpulse/internals.hpp"

using namespace std;


namespace simpulse {
#if 0
};  // pacify emacs c-mode!
#endif


// Note: if you add more constructor arguments, you'll need to change the python-linkage
// boilerplate in pybind11/wrap_phase_models.cpp

keplerian_binary_phase_model::keplerian_binary_phase_model(double e, double a, double b, double Porb, double nx, double ny, double P, double t0, double phi0)
{
    // Placeholder
}


// virtual override
double keplerian_binary_phase_model::eval_phi(double t, int nderivs) const
{
    // Placeholder
    throw runtime_error("keplerian_binary_phase_model::eval_phi() not implemented yet");
}


// virtual override
string keplerian_binary_phase_model::str() const
{
    // Placeholder.
    // Returns human-readable string representation of the phase_model object.
    // (To call it, try 'print phase_model' from python.)

    throw runtime_error("keplerian_binary_phase_model::str() not implemented yet");
}


}  // namespace simpulse
