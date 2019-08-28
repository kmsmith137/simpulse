#include "../include/simpulse/pulsar_phase_models.hpp"
#include "../include/simpulse/internals.hpp"

using namespace std;


namespace simpulse {
#if 0
};  // pacify emacs c-mode!
#endif


sinusoidal_phase_model::sinusoidal_phase_model(double pulse_phase_, double pulse_freq_, double orbital_phase_, 
					       double orbital_freq_, double beta_, double t0_)
  : pulse_phase(pulse_phase_), pulse_freq(pulse_freq_), orbital_phase(orbital_phase_), orbital_freq(orbital_freq_), beta(beta_), t0(t0_)
{}


// virtual override
double sinusoidal_phase_model::eval_phi(double t, int nderivs) const
{
  sp_assert2(nderivs <= 2, "sinusoidal_phase_model: nderivs > 2 not implemented yet!");

  t -= t0;
  double theta = orbital_freq * t + orbital_phase;
  double c = pulse_freq * beta / (2 * M_PI * orbital_freq);

  if (nderivs == 0)
    return (pulse_phase + (pulse_freq * t) + c * (sin(2 * M_PI * theta) - sin(2 * M_PI * orbital_phase)));
  if (nderivs == 1)
    return pulse_freq * (1 + beta * cos(2 * M_PI * theta));
  return -2 * M_PI * orbital_freq * pulse_freq * beta * sin(2 * M_PI * theta);
}


// virtual override
string sinusoidal_phase_model::str() const
{
  stringstream ss;
  ss << "simpulse.sinusoidal_phase_model(pulse_phase=" << pulse_phase << ", pulse_freq=" << pulse_freq << ", orbital_phase=" << orbital_phase 
     << ", orbital_freq=" << orbital_freq << ", beta=" << beta << ", t0=" << t0 << ")";
  return ss.str();
}


}  // namespace simpulse
