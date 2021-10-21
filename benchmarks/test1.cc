#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TF1.h"

#include <vector>
#include <string>
#include <cstdio>

using namespace std;

//----------------------------------------------------------------------------------------------------

void NormalisationCheck(TH1D *h_dsdt, double t_min_in, double t_max_in)
{
	printf(">> NormalisationCheck\n");

	// get sum in Coulomb region
	int bi_min = h_dsdt->GetXaxis()->FindBin(t_min_in + 1E-6);
	int bi_max = h_dsdt->GetXaxis()->FindBin(t_max_in - 1E-6);
	printf("    bin min: %i, min max: %i --> %i bins\n", bi_min, bi_max, bi_max - bi_min + 1);

	double s = 0.;
	for (int bi = bi_min; bi <= bi_max; ++bi)
		s += h_dsdt->GetBinContent(bi) * h_dsdt->GetBinWidth(bi);

	printf("    sum of bin content: %.3E\n", s);

	// get reference Coulomb cross-section
	const double t_min_be = h_dsdt->GetBinLowEdge(bi_min);
	const double t_max_be = h_dsdt->GetXaxis()->GetBinUpEdge(bi_max);
	const double a_ref = 0.000260508; // ds/dt = a / t^2
	const double si_ref = a_ref * (1./t_min_be - 1./t_max_be);

	printf("    Coulomb integral: from = %.2E, to = %.2E, si_ref = %.3E\n", t_min_be, t_max_be, si_ref);

	// get scale
	printf("    scale = %.5f\n", si_ref / s);
}

//----------------------------------------------------------------------------------------------------

int main()
{
	// get input
	TFile *f_in = TFile::Open("../data/input/version_2/import.root");

	TH1D *h_dsdt_high = (TH1D *) f_in->Get("sb1/high_beta/h_dsdt_cen_stat");
	//TH1D *h_dsdt_low = (TH1D *) f_in->Get("sb1/low_beta/h_dsdt_cen_stat");

	// prepare output
	TFile *f_out = new TFile("test1.root", "recreate");

	h_dsdt_high->Write("h_dsdt_high");

	// check at very low |t|
	NormalisationCheck(h_dsdt_high, 3E-4, 5E-4);
	NormalisationCheck(h_dsdt_high, 3E-4, 4E-4);

	// fit function
	TF1 *ff = new TF1("ff", "[0] * exp(-[1]*x)");
	//TF1 *ff = new TF1("ff", "[0] * exp(-[1]*x + [2]*x*x)");
	//TF1 *ff = new TF1("ff", "[0] * exp(-[1]*x + [2]*x*x + [3]*x*x*x)");

	// slope checks
	{
		printf("\n>> slope checks\n");

		gDirectory = f_out->mkdir("slope checks");

		for (double t_min = 0.005; t_min <= 0.020; t_min += 0.005)
		{
			const double t_max = t_min + 0.005;
			h_dsdt_high->Fit(ff, "Q", "", t_min, t_max);
			const auto B = ff->GetParameter(1), B_unc = ff->GetParError(1);
			printf("    fit from %.1E to %.1E: B = %.2f +- %.2f\n", t_min, t_max, B, B_unc);

			char buf[100];
			sprintf(buf, "from_%.1E_to_%.1E", t_min, t_max);
			h_dsdt_high->Write(buf);
		}
	}

	// extrapolation/si_tot checks
	{
		printf("\n>> extrapolation checks\n");

		TDirectory *d_ext_checks = f_out->mkdir("extrapolation checks");

		const double sq_hbarc = 3.893790E-01;
		const double rho = 0.14;

		struct FitRange { double t_min, t_max; };
		vector<FitRange> fitRanges = {
			{ 0.005, 0.015 },
			{ 0.005, 0.020 },
		};

		for (const auto &fitRange : fitRanges)
		{
			printf("    from %.3f to %.3f\n", fitRange.t_min, fitRange.t_max);

			char buf[100];
			sprintf(buf, "from_%.3f_to_%.3f", fitRange.t_min, fitRange.t_max);
			gDirectory = d_ext_checks->mkdir(buf);

			h_dsdt_high->Fit(ff, "Q", "", fitRange.t_min, fitRange.t_max);
			h_dsdt_high->Write("h_dsdt_high_with_fit");

			const double A = ff->GetParameter(0), A_unc = ff->GetParError(0);

			const double K = sqrt( 16.*M_PI * sq_hbarc / (1. + rho * rho) );
			const double si_tot_high = K * sqrt(A);
			const double si_tot_high_unc = K / 2. / sqrt(A) * A_unc;
			printf("        si_tot = %.2f +- %.2f mb\n", si_tot_high, si_tot_high_unc);

			TH1D *h_results = new TH1D("h_results", "", 8, 0.5, 8.5);
			TAxis *xa = h_results->GetXaxis();
			xa->SetBinLabel(1, "A"); h_results->SetBinContent(1, ff->GetParameter(0));
			xa->SetBinLabel(2, "A_unc"); h_results->SetBinContent(2, ff->GetParError(0));
			xa->SetBinLabel(3, "B"); h_results->SetBinContent(3, ff->GetParameter(1));
			xa->SetBinLabel(4, "B_unc"); h_results->SetBinContent(4, ff->GetParError(1));
			xa->SetBinLabel(5, "rho"); h_results->SetBinContent(5, rho);
			xa->SetBinLabel(6, "rho_unc"); h_results->SetBinContent(6, 0.);
			xa->SetBinLabel(7, "si_tot"); h_results->SetBinContent(7, si_tot_high);
			xa->SetBinLabel(8, "si_tot_unc"); h_results->SetBinContent(8, si_tot_high_unc);
			h_results->Write();
		}
	}

	// low beta fit
	/*
	printf("\n----- low beta -----\n");
	h_low->Fit(ff, "", "", 0.017, 0.053);
	h_low->Write("h_low");

	const double si_tot_low = sqrt( 16.*M_PI * sq_hbarc / (1. + rho * rho) * ff->GetParameter(0) );
	printf("si_tot = %.2f\n", si_tot_low);
	*/

	delete f_out;

	return 0;
}
