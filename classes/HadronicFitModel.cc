#include "classes/HadronicFitModel.hh"

#include "Elegent/Constants.h"
#include "Elegent/Math.h"

//----------------------------------------------------------------------------------------------------

double HadronicFitModel::DsigmaDTHighT(double t) const
{
	double x = -t;

	// exp4+exp2
	if (modulusHighTVariant == 1)
	{
		const double P0 = 6.11442e+02;
		const double P1 = -2.07544e+01;
		const double P2 = 1.01559e+00;
		const double P3 = 2.23444e+01;
		const double P4 = -9.65895e+01;
		const double P5 = 2.89226e-04;
		const double P6 = 1.44707e+01;
		const double P7 = -1.09700e+01;

		return P0*exp(P1*x + P2*x*x + P3*x*x*x + P4*x*x*x*x) + P5*exp(P6*x + P7*x*x);
	}

	// p1*exp3+p1*exp1
	if (modulusHighTVariant == 2)
	{
		const double P0 = 6.24949e+02;
		const double P1 = -2.56314e+02;
		const double P2 = -2.04532e+01;
		const double P3 = 8.49336e+00;
		const double P4 = -1.60850e+01;
		const double P5 = -1.11034e+01;
		const double P6 = 2.25886e+01;
		const double P7 = -7.02090e+00;

		return (P0 + P1*x) * exp(P2*x + P3*x*x + P4*x*x*x) + (P5 + P6*x) * exp(P7*x);
	}

	// p1*exp3+p2*exp2
	if (modulusHighTVariant == 3)
	{
		const double P0 = 7.16305e+02;
		const double P1 = -2.37871e+02;
		const double P2 = -1.96623e+01;
		const double P3 = 9.34281e+00;
		const double P4 = -1.50302e+01;
		const double P5 = -1.02707e+02;
		const double P6 = 8.08324e+01;
		const double P7 = 2.20613e+02;
		const double P8 = -1.29148e+01;
		const double P9 = 3.09810e+00;

		return (P0 + P1*x) * exp(P2*x + P3*x*x + P4*x*x*x) + (P5 + P6*x + P7*x*x) * exp(P8*x + P9*x*x);
	}

	// exp3-intf-exp1
	if (modulusHighTVariant == 4)
	{
		const double P0 = 2.65511e+00;
		const double P1 = 2.55649e+01;
		const double P2 = -1.02703e+01;
		const double P3 = 4.42715e+00;
		const double P4 = -6.83600e+00;
		const double P5 = 9.00437e-01;
		const double P6 = -2.16005e+00;

		return P1*P1*exp(2*P2*x + 2*P3*x*x + 2*P4*x*x*x) + 2 * cos(P0) * P1*exp(P2*x + P3*x*x + P4*x*x*x) * P5*exp(P6*x) + P5*P5*exp(2*P6*x);
	}

	return 0;
}

//----------------------------------------------------------------------------------------------------

TComplex HadronicFitModel::Amp(double t) const
{
	// ---------- modulus ----------

	double m = 0.;

	if (modulusMode == mmExp)
	{
		// blending region
		// TODO: necessary ?
		/*
		double t_avg = (t2 + t1)/2., t_si = (t2 - t_avg) / 3.;
		double t_opt_m1 = -t_avg + 5. * t_si;
		double t_opt_m2 = -t_avg - 5. * t_si;
		*/

		// main modulus part (for low |t|)
		double bPol = 0., tPow = t;
		bPol += b1 * tPow; tPow *= t;
		bPol += b2 * tPow; tPow *= t;
		bPol += b3 * tPow; tPow *= t;
		bPol += b4 * tPow; tPow *= t;
		bPol += b5 * tPow; tPow *= t;
		bPol += b6 * tPow; tPow *= t;
		bPol += b7 * tPow; tPow *= t;
		bPol += b8 * tPow; tPow *= t;
		bPol += b9 * tPow;
		const double m1 = a * exp(bPol);

		// TODO: check if necessary
		m = m1;
		/*
		if (t > t_opt_m1)
		{
			// optimisation for very low |t|: only m1 component
			m = m1;
		} else {	
			// fixed part for higher |t|
			double m2 = hts * sqrt(DsigmaDTHighT(t) / cnts->sig_fac);

			if (t < t_opt_m2)
			{
				// optimisation for very high |t|: only m2 component
				m = m2;
			} else {
				// full modulus (low and high-|t| parts blended)
				double p = (TMath::Erf( (-t - t_avg) / t_si / sqrt(2.) ) + 1.) / 2.;
				m = m1*(1.-p) + m2*p;
			}
		}
		*/
	}

	// ---------- phase ----------

	double ph = 0.;

	if (phaseMode == pmUnknown)
		printf("ERROR in HadronicFitModel::Amp > phase mode is `unknown'.\n");

	if (phaseMode == pmConstant)
		ph = p0;

	if (phaseMode == pmBailly)
	{
		double rho = cos(p0)/sin(p0);

		double ze = 0.;
		if (rho == 0.)
			ze = 0.;
		else {
			double r = t / p_td;
			double ep = 1E-4;
			if (r < 1. - ep) ze = atan(rho / (1. - r));
			if (r > 1. - ep && r < 1. + ep) ze = M_PI/2. + (r-1.)/rho;
			if (r > 1. + ep) ze = atan(rho / (1. - r)) + M_PI;
		}

		ph = M_PI/2. - ze;
	}

	if (phaseMode == pmStandard)
	{
		ph = p0 + atan( (fabs(t) - fabs(p_t0)) / p_tau ) - atan(- fabs(p_t0) / p_tau);
	}

	if (phaseMode == pmPeripheral)
	{
		double r = fabs(t / p_tm);
		ph = p0 - p_A * exp(p_ka * (log(r) - r + 1.));
	}

	if (phaseMode == pmPolynomial)
	{
		double ph_pol = p1*t + p2*t*t + p3*t*t*t;

		double mt = -t;
		double mu = 1., si = 0.2;
		double w = (1. - TMath::Erf( (mt - mu)/(si * 1.414) )) / 2.;

		ph = p0 + ph_pol * w;
	}

	return m * TComplex::Exp(i * ph);
}

//----------------------------------------------------------------------------------------------------

TComplex HadronicFitModel::Amp_J0(double t, double *par, const void *vobj)
{
	const HadronicFitModel *obj = (HadronicFitModel *) vobj;
	const double &b = par[0];	// impact parameter in GeV^-1

	return obj->Amp(t) * TMath::BesselJ0(b * sqrt(-t));
}

//----------------------------------------------------------------------------------------------------

TComplex HadronicFitModel::Prf(double b_fm) const
{
	double b = b_fm / cnts->hbarc;	// b in GeV^-1
	double par[] = { b };

	TComplex I = ComplexIntegrate(Amp_J0, par, this, upper_bound_t, 0., 0., precision_t,
		integ_workspace_size_t, integ_workspace_t, "HadronicFitModel::Prf");

	return I / 4. / cnts->p_cms / cnts->sqrt_s;
}
