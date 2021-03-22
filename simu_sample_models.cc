#include "classes/command_line_tools.hh"

#include "Elegent/Constants.h"
#include "Elegent/CoulombInterference.h"

#include "classes/HadronicFitModel.hh"

#include "TFile.h"
#include "TGraph.h"

#include <cstring>

using namespace std;
using namespace Elegent;

//----------------------------------------------------------------------------------------------------

void PrintUsage()
{
	printf("USAGE: program <option> <option>\n");
	printf("OPTIONS:\n");
	//printf("    -cfg <file>       config file\n");
}

//----------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------

void SampleModel(Model *m, const string &label)
{
	model = m;

	TDirectory *d_top = gDirectory;

	gDirectory = d_top->mkdir(label.c_str());

	TGraph *g_dsdt_H = new TGraph();
	TGraph *g_dsdt_CH = new TGraph();

	for (double mt = 1E-4; mt <= 0.11; )
	{
		coulomb->mode = CoulombInterference::mPH;
		TComplex F_H = coulomb->Amp(-mt);

		coulomb->mode = CoulombInterference::mKL;
		TComplex F_CH = coulomb->Amp(-mt);

		int idx = g_dsdt_H->GetN();
		g_dsdt_H->SetPoint(idx, mt, cnts->sig_fac * F_H.Rho2());
		g_dsdt_CH->SetPoint(idx, mt, cnts->sig_fac * F_CH.Rho2());

		double dmt = 5E-3;
		if (mt < 0.10) dmt = 1E-3;
		if (mt < 0.004) dmt = 1E-4;
		if (mt < 0.001) dmt = 1E-5;
		mt += dmt;
	}

	g_dsdt_H->Write("g_dsdt_H");
	g_dsdt_CH->Write("g_dsdt_CH");

	gDirectory = d_top;
}

//----------------------------------------------------------------------------------------------------

int main(int argc, const char **argv)
{
	// defaults
	//string cfg_file = "config.py";

	// parse command line
	for (int argi = 1; (argi < argc) && (cl_error == 0); ++argi)
	{
		if (strcmp(argv[argi], "-h") == 0 || strcmp(argv[argi], "--help") == 0)
		{
			cl_error = 1;
			continue;
		}

		//if (TestStringParameter(argc, argv, argi, "-cfg", cfg_file)) continue;

		printf("ERROR: unknown option '%s'.\n", argv[argi]);
		cl_error = 1;
	}

	if (cl_error)
	{
		PrintUsage();
		return 1;
	}

	// init Elegent
	Constants::Init(2*450, cnts->mPP);
    cnts->Print();

	coulomb->ffType = coulomb->ffPuckett;
	coulomb->Print();

	// prepare output
	TFile *f_out = TFile::Open("models.root", "recreate");

	// sample models
	HadronicFitModel *hfm = new HadronicFitModel();

	coulomb->mode = CoulombInterference::mKL;

	hfm->a = 5.67E6;
	hfm->b1 = 8.5;
	hfm->b2 = 0.;
	hfm->b3 = 0.;

	hfm->t1 = 0.2;
	hfm->t2 = 0.5;

	hfm->phaseMode = HadronicFitModel::pmConstant;

	hfm->p0 = M_PI/2. - atan(0.05);
	SampleModel(hfm, "exp1_con_rho0.05");

	hfm->p0 = M_PI/2. - atan(0.10);
	SampleModel(hfm, "exp1_con_rho0.10");

	hfm->p0 = M_PI/2. - atan(0.15);
	SampleModel(hfm, "exp1_con_rho0.15");

	// clean up
	delete f_out;

	return 0;
}
