#include "classes/HadronicFitModel.hh"

#include "Elegent/Constants.h"
#include "Elegent/Math.h"

using namespace Elegent;

//----------------------------------------------------------------------------------------------------

double HadronicFitModel::DsigmaDTHighT(double t) const
{
	const double mt = -t;

	// exp1
	if (modulusHighTVariant == mhtVariant1)
	{
		const double A = 0.5E-2;
		const double B = 1.7;

		return A*exp(-B*mt);
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
		const double t_avg = (t2 + t1)/2., t_si = (t2 - t_avg) / 3.;
		const double t_opt_m1 = -t_avg + 5. * t_si;
		const double t_opt_m2 = -t_avg - 5. * t_si;

		// main modulus part (for low |t|)
		double bPol = 0., tPow = t;
		bPol += b1 * tPow; tPow *= t;
		bPol += b2 * tPow; tPow *= t;
		bPol += b3 * tPow; tPow *= t;
		bPol += b4 * tPow; tPow *= t;
		bPol += b5 * tPow;
		const double m1 = a * exp(bPol);

		if (t > t_opt_m1)
		{
			// optimisation for very low |t|: only m1 component
			m = m1;
		} else {	
			// fixed part for higher |t|
			const double m2 = hts * sqrt(DsigmaDTHighT(t) / cnts->sig_fac);

			if (t < t_opt_m2)
			{
				// optimisation for very high |t|: only m2 component
				m = m2;
			} else {
				// full modulus (low and high-|t| parts blended)
				const double p = (TMath::Erf( (-t - t_avg) / t_si / sqrt(2.) ) + 1.) / 2.;
				m = m1*(1.-p) + m2*p;
			}
		}
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
