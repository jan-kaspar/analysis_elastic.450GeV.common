#ifndef _HadronicFitModel_hh_
#define _HadronicFitModel_hh_

#include "Elegent/Model.h"

#include <string>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

//----------------------------------------------------------------------------------------------------

class HadronicFitModel : public Elegent::Model
{
  public:
	/// modulus parameters (low |t|)
	enum ModulusMode { mmUnknown, mmExp } modulusMode;
	double a, b1, b2, b3, b4, b5, b6, b7, b8, b9;

	/// modulus parameters (high |t|)
	unsigned int modulusHighTVariant;
	double hts; ///< scale factor for the high-|t| part

	/// parameters for blending the low-|t| (variable) and high-|t| (fixed) modulus
	/// the interval (t1, t2) corresponds to (-3 sigma, +3 sigma)
	double t1, t2;

	/// phase parameters
	enum PhaseMode { pmUnknown, pmConstant, pmBailly, pmStandard, pmPeripheral, pmPolynomial } phaseMode;
	double p0, p1, p2, p3, p_td, p_t0, p_tau, p_A, p_ka, p_tm;

    HadronicFitModel() :
		modulusMode(mmExp),
		a(0.), b1(0.), b2(0.), b3(0.), b4(0.), b5(0.), b6(0.), b7(0.), b8(0.), b9(0.),

		modulusHighTVariant(2),
		hts(1.),
		t1(0.2), t2(0.5),

		phaseMode(pmUnknown),
		p0(0.), p1(0.), p2(0.), p3(0.), p_td(0.), p_t0(0.), p_tau(0.), p_A(0.), p_ka(0.), p_tm(0.),
		integ_workspace_initialized(false)
	{
	}

    ~HadronicFitModel() override {}

	void Init() override
	{
		fullLabel.name = "hadronic fit model";
		shortLabel.name = "hfm";

		precision_t = 1E-4;
		upper_bound_t = -50.;

		if (!integ_workspace_initialized)
		{
			integ_workspace_size_t = 100;
			integ_workspace_t = gsl_integration_workspace_alloc(integ_workspace_size_t);
			integ_workspace_initialized = true;
		}
	}

    void Print() const override
	{
		printf(">> HadronicFitModel\n");

		printf("\tmodulus: %s\n", GetModulusModeString().c_str());
		if (modulusMode == mmExp)
			printf("\t\ta=%.3E, b1=%.3E, b2=%.3E, b3=%.3E, b4=%.3E, b5=%.3E, b6=%.3E, b7=%.3E, b8=%.3E, b9=%.3E\n",
				a, b1, b2, b3, b4, b5, b6, b7, b8, b9);

		printf("\tmodulus at high |t|: variant %u, hts=%.3f\n", modulusHighTVariant, hts);

		printf("\tblending parameters: t1=%.3E, t2=%.3E\n", t1, t2);

		printf("\tphase: %s\n", GetPhaseModeString().c_str());
		if (phaseMode == pmConstant)
			printf("\t\tp0=%.3E\n", p0);
		if (phaseMode == pmBailly)
			printf("\t\tp0=%.3E, p_td=%.3E\n", p0, p_td);
		if (phaseMode == pmStandard)
			printf("\t\tp0=%.3E, p_t0=%.3E, p_tau=%.3E\n", p0, p_t0, p_tau);
		if (phaseMode == pmPeripheral)
			printf("\t\tp0=%.3E, p_A=%.3E, p_ka=%.3E, p_tm=%.3E\n", p0, p_A, p_ka, p_tm);
		if (phaseMode == pmPolynomial)
			printf("\t\tp0=%.3E, p1=%.3E, p2=%.3E, p3=%.3E\n", p0, p1, p2, p3);
	}

	void PrintCode(const char *name = "hfm") const
	{
		if (modulusMode == mmExp)
			printf("%s->modulusMode = HadronicFitModel::mmExp;\n%s->a=%.5E;\n%s->b1=%.5E;\n%s->b2=%.5E;\n%s->b3=%.5E;\n%s->b4=%.5E;\n%s->b5=%.5E;\n%s->b6=%.5E;\n%s->b7=%.5E;\n%s->b8=%.5E;\n%s->b9=%.5E;\n",
				name, name, a, name, b1, name, b2, name, b3, name, b4, name, b5, name, b6, name, b7, name, b8, name, b9);

		printf("%s->modulusHighTVariant = %u\n", name, modulusHighTVariant);
		printf("%s->hts = %.1E\n", name, hts);

		if (modulusMode == mmExp)
			printf("%s->t1=%.5E;\n%s->t2=%.5E;\n", name, t1, name, t2);

		printf("\n");

		if (phaseMode == pmConstant)
			printf("%s->phaseMode = HadronicFitModel::pmConstant;\n%s->p0=%.5E;\n", name, name, p0);
		if (phaseMode == pmBailly)
			printf("%s->phaseMode = HadronicFitModel::pmBailly;\n%s->p0=%.5E;\n%s->p_td=%.5E;\n", name, name, p0, name, p_td);
		if (phaseMode == pmStandard)
			printf("%s->phaseMode = HadronicFitModel::pmStandard;\n%s->p0=%.5E;\n%s->p_t0=%.5E;\n%s->p_tau=%.5E;\n",
				name, name, p0, name, p_t0, name, p_tau);
		if (phaseMode == pmPeripheral)
			printf("%s->phaseMode = HadronicFitModel::pmPeripheral;\n%s->p0=%.5E;\n%s->p_A=%.5E;\n%s->p_ka=%.5E;\n%s->p_tm=%.5E;\n",
				name, name, p0, name, p_A, name, p_ka, name, p_tm);
		if (phaseMode == pmPolynomial)
			printf("%s->phaseMode = HadronicFitModel::pmPolynomial;\n%s->p0=%.5E;\n%s->p1=%.5E;\n%s->p2=%.5E;\n%s->p3=%.5E;\n",
				name, name, p0, name, p1, name, p2, name, p3);
	}

	std::string GetStateString() const
	{
		char buf_m[200];
		if (modulusMode == mmExp)
			sprintf(buf_m, "modulus: exp, a=%.3E, b1=%.3E, b2=%.3E, b3=%.3E, b4=%.3E, b5=%.3E, b6=%.3E, b7=%.3E, b8=%.3E, b9=%.3E",
				a, b1, b2, b3, b4, b5, b6, b7, b8, b9);

		char buf_p[200];
		if (phaseMode == pmConstant)
			sprintf(buf_p, "phase: constant, p0=%.3E", p0);
		if (phaseMode == pmBailly)
			sprintf(buf_p, "phase: Bailly, p0=%.3E, p_td=%.3E", p0, p_td);
		if (phaseMode == pmStandard)
			sprintf(buf_p, "phase: standard, p0=%.3E, p_t0=%.3E, p_tau=%.3E", p0, p_t0, p_tau);
		if (phaseMode == pmPeripheral)
			sprintf(buf_p, "phase: peripheral, p0=%.3E, p_A=%.3E, p_ka=%.3E, p_tm=%.3E", p0, p_A, p_ka, p_tm);
		if (phaseMode == pmPolynomial)
			sprintf(buf_p, "phase: polynomial, p0=%.3E, p1=%.3E, p2=%.3E, p3=%.3E", p0, p1, p2, p3);

		return std::string(buf_m) + "; " + buf_p;
	}

	std::string GetStateStringShort() const
	{
		char buf_m[200];
		if (modulusMode == mmExp)
			sprintf(buf_m, "a=%.3E, b1=%.3E, b2=%.3E, b3=%.3E, b4=%.3E, b5=%.3E, b6=%.3E, b7=%.3E, b8=%.3E, b9=%.3E | ", a, b1, b2, b3, b4, b5, b6, b7, b8, b9);

		char buf_p[200];
		if (phaseMode == pmConstant)
			sprintf(buf_p, "p0=%.3E", p0);
		if (phaseMode == pmBailly)
			sprintf(buf_p, "p0=%.3E, p_td=%.3E", p0, p_td);
		if (phaseMode == pmStandard)
			sprintf(buf_p, "p0=%.3E, p_t0=%.3E, p_tau=%.3E", p0, p_t0, p_tau);
		if (phaseMode == pmPeripheral)
			sprintf(buf_p, "p0=%.3E, p_A=%.3E, p_ka=%.3E, p_tm=%.3E", p0, p_A, p_ka, p_tm);
		if (phaseMode == pmPolynomial)
			sprintf(buf_p, "p0=%.3E, p1=%.3E, p2=%.3E, p3=%.3E", p0, p1, p2, p3);

		return std::string(buf_m) + buf_p;
	}


    std::string GetModulusModeString() const
	{
		if (modulusMode == mmExp) return "exp";

		return "unknown";
	}

    std::string GetPhaseModeString() const
	{
		if (phaseMode == pmConstant) return "constant";
		if (phaseMode == pmBailly) return "Bailly";
		if (phaseMode == pmStandard) return "standard";
		if (phaseMode == pmPeripheral) return "peripheral";
		if (phaseMode == pmPolynomial) return "polynomial";

		return "unknown";
	}

	unsigned int GetNumberOfPhaseParameters() const
	{
		if (phaseMode == pmConstant) return 1;
		if (phaseMode == pmBailly) return 1;
		if (phaseMode == pmStandard) return 1;
		//if (phaseMode == pmPeripheral) return 4;
		if (phaseMode == pmPeripheral) return 1;
		if (phaseMode == pmPolynomial) return 1;
		return 0;
	}

	// TODO: needed ?
	double DsigmaDTHighT(double t) const;

    TComplex Amp(double t) const override;

	double upper_bound_t, precision_t;

	bool integ_workspace_initialized;
	unsigned long integ_workspace_size_t;
	gsl_integration_workspace *integ_workspace_t;

	static TComplex Amp_J0(double t, double *par, const void *vobj);
    TComplex Prf(double b_fm) const override;
};

#endif
