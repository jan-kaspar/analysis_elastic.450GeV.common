#include "classes/HadronicFitModel.hh"

#include "Elegent/Constants.h"
#include "Elegent/CoulombInterference.h"
#include "Elegent/Model.h"

#include "TFile.h"
#include "TGraph.h"
#include "TComplex.h"
#include "TDirectory.h"
#include "TF1.h"

#include <cstdio>
#include <cstring>

using namespace std;
using namespace Elegent;

//----------------------------------------------------------------------------------------------------

void PrintUsage()
{
	printf("USAGE: program <option> <option>\n");
	printf("OPTIONS:\n");
}

//----------------------------------------------------------------------------------------------------

TF1 *ff = new TF1("ff", "[0] * exp(-[1]*x)");

void AnalyzeModel(Model *m, const string &cniMode)
{
	model = m;

	// make graph
	TGraph *g_dsdt_CH = new TGraph();

	for (double mt = 0.9E-4; mt <= 0.03; )
	{
		TComplex F_CH = 0;

		if (cniMode == "AmpSum")
			F_CH = coulomb->Amp_pure(-mt) + model->Amp(-mt);

		if (cniMode == "KL")
			F_CH = coulomb->Amp_KL(-mt);

		int idx = g_dsdt_CH->GetN();
		g_dsdt_CH->SetPoint(idx, mt, cnts->sig_fac * F_CH.Rho2());

		double dmt = 0.1;
		if (mt < 0.10) dmt = 5E-4;
		if (mt < 0.004) dmt = 1E-4;
		if (mt < 0.001) dmt = 5E-6;
		mt += dmt;
	}

	g_dsdt_CH->Write("g_dsdt_CH");

	// make fits
	vector<pair<double, double>> ranges = {
		{ 0.005, 0.010 },
		{ 0.010, 0.015 },
		{ 0.015, 0.020 },
		{ 0.020, 0.025 },
		{ 0.005, 0.020 },
	};

	for (const auto &r : ranges)
	{
		ff->SetParameters(200, 17);
		ff->SetRange(r.first, r.second);
		g_dsdt_CH->Fit(ff, "Q", "", r.first, r.second);

		char buf[100];
		sprintf(buf, "fit_%.3f_%.3f", r.first, r.second);
		ff->Write(buf);
	}
}

//----------------------------------------------------------------------------------------------------

int main()
{
	// init Elegent
	Constants::Init(2*450, cnts->mPP);
    cnts->Print();

	coulomb->ffType = coulomb->ffPuckett;
	coulomb->Print();

	// prepare output
	TFile *f_out = TFile::Open("coulomb_effect_on_B.root", "recreate");

	// sample models
	HadronicFitModel *hfm = new HadronicFitModel();

	hfm->modulusMode = HadronicFitModel::mmExp;
	hfm->a = 5.67E6;
	hfm->b1 = 8.5;
	hfm->b2 = 0.;
	hfm->b3 = 0.;

	hfm->t1 = 0.5;
	hfm->t2 = 1.5;

	hfm->phaseMode = HadronicFitModel::pmConstant;

	coulomb->mode = CoulombInterference::mKL;

	for (const string &cniMode : {"KL", "AmpSum"})
	{
		for (const double rho : {0.10, 0.15})
		{
			char buf[100];
			sprintf(buf, "exp1_con_rho%.2f_%s", rho, cniMode.c_str());
			gDirectory = f_out->mkdir(buf);

			hfm->p0 = M_PI/2. - atan(rho);

			AnalyzeModel(hfm, cniMode);
		}
	}

	// clean up
	delete f_out;

	return 0;
}
