#include "../include/simpulse/pulsar_phase_models.hpp"
#include "../include/simpulse/internals.hpp"
#include <iostream>

using namespace std;


namespace simpulse {
#if 0
};  // pacify emacs c-mode!
#endif


keplerian_binary_phase_model::keplerian_binary_phase_model(double e_, double a_, double Porb_, double nx_, double ny_, double P_, double t0_, double phi0_)
    : e(e_), a(a_), b(xsqrt(a*a*(1-e*e))), Porb(Porb_), nx(nx_), ny(ny_), P(P_), t0(t0_), phi0(phi0_)
{
    sp_assert2((e >= 0) && (e < 1), "simpulse::keplerian_binary_phase_model: e must be in [0, 1)");
    sp_assert2(sqrt(nx * nx + ny * ny) <= 1, "simpulse::keplerian_binary_phase_model: nx and ny must be components of a unit vector");
    sp_assert2(max(a, b) < Porb / 10, "simpulse::keplerian_binary_phase_model: a and b must be much smaller than Porb");
}


// Newton's method to solve for E given a tobs
double keplerian_binary_phase_model::find_E(double tobs) const
{
    // We precompute coefficients a, b, c, d so that the equation being solved is
    // f(E) = C E - A cos(E) - B sin(E) + D = 0
    // This precomputation speeds things up a little.
    double C = Porb / 2 / M_PI;
    double D = t0 + a * nx * e - tobs;
    double A = a * nx;
    double B = b * ny + C * e;

    double E = -D / C;
    double thresh = 1.0e-7 * P;

    for (int i = 0; i < 20; i++)
    {
	double cosE = cos(E);
	double sinE = sin(E);

	// Left-hand side f(E)
	double f_E = C * E - A * cosE - B * sinE + D;
	if (abs(f_E) < thresh)
	    return E;
	
	// Derivative df/dE
	double df_dE = C + A * sinE - B * cosE;

	// Newton's method update
	E -= f_E / df_dE;
    }

    throw runtime_error("keplerian_binary_phase_model::find_E failed to converge (should never happen)");
}

// virtual override
double keplerian_binary_phase_model::eval_phi(double t, int nderivs) const
{
    double E = find_E(t);
    
    if (nderivs == 0)
	return Porb / 2 / M_PI * (E - e * sin(E)) / P + phi0;

    double dEdtobs = 1 / (Porb / 2 / M_PI * (1 - e * cos(E)) + a * nx * sin(E) - b * ny * cos(E));
    if (nderivs == 1)
	return Porb / 2 / M_PI / P * (1 - e * cos(E)) * dEdtobs;
    if (nderivs == 2)
    {
        return Porb / 2 / M_PI / P * pow(dEdtobs, 2) * 
	       (e * sin(E) - (1 - e * cos(E)) * (Porb / 2 / M_PI * e * sin(E)  + a * nx * cos(E) + b * ny * sin(E)) * dEdtobs);
    }
    throw runtime_error("keplerian_binary_phase_model::eval_phi() nderivs != 0, 1, 2 not implemented yet");
}


// virtual override
string keplerian_binary_phase_model::str() const
{
    stringstream ss;
    ss << "simpulse.keplerian_binary_phase_model(e=" << e << ", a=" << a << ", b=" << b << ", Porb=" <<
	Porb << ", nx=" << nx << ", ny=" << ny << ", P=" << P << ", t0=" << t0 << ", phi0=" << phi0 << ")";
    return ss.str();
}

}  // namespace simpulse
