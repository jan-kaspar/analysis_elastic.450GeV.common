#include "TFile.h"
#include "TGraph.h"
#include "TF1.h"

#include <vector>
#include <string>
#include <cstdio>

using namespace std;

//----------------------------------------------------------------------------------------------------

TGraph* MakeDiff(const TGraph *g1, const TGraph *g2)
{
	TGraph *g_res = new TGraph();

	for (int i = 0; i < g1->GetN(); ++i)
	{
		double x1, y1;
		g1->GetPoint(i, x1, y1);

		if (x1 > 0.1)
			continue;

		double y2 = g2->Eval(x1);

		g_res->SetPoint(g_res->GetN(), x1, y1 - y2);
	}

	return g_res;
}

//----------------------------------------------------------------------------------------------------

int main()
{
	// get input
	TFile *f_in = TFile::Open("../data/simu/iv2_TEST/models.root");

	TGraph *g_mod = (TGraph *) f_in->Get("exp1_con_rho0.15/g_dsdt_CH");
	TGraph *g_ref = (TGraph *) f_in->Get("exp1_con_rho0.10/g_dsdt_CH");

	TGraph *g_diff = MakeDiff(g_mod, g_ref);

	// prepare output
	TFile *f_out = new TFile("test2.root", "recreate");

	TF1 *ff = new TF1("ff", "[0] + [1]*x");

	g_diff->Fit(ff, "", "", 0.005, 0.020);
	g_diff->Write("fit_0.005_0.020");

	g_diff->Fit(ff, "", "", 0.015, 0.025);
	g_diff->Write("fit_0.015_0.025");

	delete f_out;

	return 0;
}
