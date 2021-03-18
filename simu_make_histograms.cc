#include "classes/command_line_tools.hh"

#include "TDirectory.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TRandom.h"
#include "TMatrixDSymEigen.h"

#include <cstdio>
#include <cstring>
#include <memory>
#include <sstream>
#include <map>

using namespace std;

//----------------------------------------------------------------------------------------------------

void BuildModelHistogram(const TGraph *g_model, const TH1D *h_input, int bi_min, int bi_max, TH1D *h_simu_ideal)
{
	for (int bi = bi_min; bi <= bi_max; ++bi)
	{
		const double l = h_input->GetBinLowEdge(bi);
		const double w = h_input->GetBinWidth(bi);
		const double input_c = h_input->GetBinContent(bi);
		const double input_u = h_input->GetBinError(bi);

		const unsigned int n_div = 10;
		double c = 0.;
		for (unsigned j = 0; j < n_div; ++j)
		{
			const double t = l + (double(j) + 0.5) * w / n_div;
			c += g_model->Eval(t) / n_div;
		}

		h_simu_ideal->SetBinContent(bi, c);
		h_simu_ideal->SetBinError(bi, input_u / input_c * c);
	}
}

//----------------------------------------------------------------------------------------------------

void AddStatisticalErrors(int bi_min, int bi_max, TH1D *h_simu)
{
	for (int bi = bi_min; bi <= bi_max; ++bi)
		h_simu->SetBinContent(bi, h_simu->GetBinContent(bi) + gRandom->Gaus() * h_simu->GetBinError(bi));
}

//----------------------------------------------------------------------------------------------------

void AddSystematicErrors(const TMatrixD &m_syst_unc, const TH1D *h_simu_ideal, int bi_min, int bi_max, TH1D *h_simu)
{
	// generate relative bin perturbations
	int n_bins = bi_max - bi_min + 1;

	TMatrixDSym m_red(n_bins);
	for (int i = 0; i < n_bins; ++i)
	{
		for (int j = 0; j < n_bins; ++j)
			m_red(i, j) = m_syst_unc(i + bi_min-1, j + bi_min-1);
	}

	TMatrixDSymEigen eig_decomp(m_red);
	const TVectorD& eig_values(eig_decomp.GetEigenValues());
	TMatrixDSym S(n_bins);
	for (int i = 0; i < n_bins; i++)
		S(i, i) = (eig_values(i) >= 0.) ? sqrt(eig_values(i)) : 0.;
	TMatrixD m_gen = eig_decomp.GetEigenVectors() * S;

	TVectorD g(n_bins);
	for (int i = 0; i < n_bins; ++i)
		g(i) = gRandom->Gaus();
	TVectorD rel_unc_red = m_gen * g;

	// apply relative bin perturbations
	for (int i = 0; i < n_bins; ++i)
	{
		const int bi = i + bi_min;
		h_simu->SetBinContent(bi, h_simu->GetBinContent(bi) + rel_unc_red(i) * h_simu_ideal->GetBinContent(bi));
	}
}

//----------------------------------------------------------------------------------------------------

void AddNormalisationError(double normalisation_error, int bi_min, int bi_max, TH1D *h_simu)
{
	for (int bi = bi_min; bi <= bi_max; ++bi)
	{
		h_simu->SetBinContent(bi, h_simu->GetBinContent(bi) * (1. + normalisation_error));
		h_simu->SetBinError(bi, h_simu->GetBinError(bi) * (1. + normalisation_error));
	}
}

//----------------------------------------------------------------------------------------------------

void PrintUsage()
{
	printf("USAGE: program <option> <option>\n");
	printf("OPTIONS:\n");
	printf("    -model-file <string>    file with model graphs");
	printf("    -model-obj <string>     sampled model object");
	printf("    -input-file <string>    file with real-data input");
	printf("    -binnings <string>      comma-separated list of binnings");
	printf("    -seed <integer>         random seed");
	printf("    -sim_stat_err <bool>    whether to simulate statistical errors");
	printf("    -sim_syst_err <bool>    whether to simulate systematic errors");
	printf("    -sim_norm_err <bool>    whether to simulate normalisation errors");
}

//----------------------------------------------------------------------------------------------------

int main(int argc, const char **argv)
{
	// defaults
	string model_file = "";
	string model_obj = "";
	string input_file = "";
	string binnings_str = "";
	unsigned int seed = 0;
	bool sim_stat_err = false;
	bool sim_syst_err = false;
	bool sim_norm_err = false;

	// parse command line
	for (int argi = 1; (argi < argc) && (cl_error == 0); ++argi)
	{
		if (strcmp(argv[argi], "-h") == 0 || strcmp(argv[argi], "--help") == 0)
		{
			cl_error = 1;
			continue;
		}

		if (TestStringParameter(argc, argv, argi, "-model-file", model_file)) continue;
		if (TestStringParameter(argc, argv, argi, "-model-obj", model_obj)) continue;

		if (TestStringParameter(argc, argv, argi, "-input-file", input_file)) continue;
		if (TestStringParameter(argc, argv, argi, "-binnings", binnings_str)) continue;

		if (TestUIntParameter(argc, argv, argi, "-seed", seed)) continue;

		if (TestBoolParameter(argc, argv, argi, "-sim-stat-err", sim_stat_err)) continue;
		if (TestBoolParameter(argc, argv, argi, "-sim-syst-err", sim_syst_err)) continue;
		if (TestBoolParameter(argc, argv, argi, "-sim-norm-err", sim_norm_err)) continue;

		printf("ERROR: unknown option '%s'.\n", argv[argi]);
		cl_error = 1;
	}

	if (cl_error)
	{
		PrintUsage();
		return 1;
	}

	// parse list of binnings
	vector<string> binnings;
	{
		stringstream ss(binnings_str);
		string binning;
		while (getline(ss, binning, ','))
			binnings.push_back(binning);
	}

	// datasets
	struct DatasetInfo
	{
		double t_min, t_max;
	};

	map<string, DatasetInfo> datasets = {
		{ "high_beta", {2.1E-4, 0.025} },
		{ "low_beta", {0.015, 0.09} },
	};

	// open input files
	unique_ptr<TFile> f_model(TFile::Open(model_file.c_str()));
	unique_ptr<TFile> f_input(TFile::Open(input_file.c_str()));

	unique_ptr<TFile> f_out(TFile::Open("histograms.root", "recreate"));

	if (!f_model || !f_input || !f_out)
	{
		printf("ERROR: cannot open some of the files.\n");
		return 2;
	}

	// load input
	TGraph *g_model = (TGraph *) f_model->Get(model_obj.c_str());

	// set seed
	gRandom->SetSeed(seed);

	// generate normalisation error
	const double normalisation_error = gRandom->Gaus() * 0.10;

	// loop over binnings
	for (const auto &binning : binnings)
	{
		printf("* %s\n", binning.c_str());

		// make binning directory
		TDirectory *d_binning = f_out->mkdir(binning.c_str());

		// loop over datasets
		for (const auto &dsit : datasets)
		{
			gDirectory = d_binning->mkdir(dsit.first.c_str());

			string input_d = binning + "/" + dsit.first;
			TH1D *h_input = (TH1D *) f_input->Get((input_d + "/h_dsdt_cen_stat").c_str());
			TMatrixD *m_syst_unc = (TMatrixD *) f_input->Get((input_d + "/m_dsdt_rel_syst").c_str());

			int bi_min = h_input->FindBin(dsit.second.t_min);
			int bi_max = h_input->FindBin(dsit.second.t_max);

			TH1D *h_simu_ideal = new TH1D(*h_input);
			h_simu_ideal->Reset();
			BuildModelHistogram(g_model, h_input, bi_min, bi_max, h_simu_ideal);

			TH1D *h_simu = new TH1D(*h_simu_ideal);

			if (sim_stat_err)
				AddStatisticalErrors(bi_min, bi_max, h_simu);

			if (sim_syst_err)
				AddSystematicErrors(*m_syst_unc, h_simu_ideal, bi_min, bi_max, h_simu);

			if (sim_norm_err)
				AddNormalisationError(normalisation_error, bi_min, bi_max, h_simu);

			h_simu->Write("h_dsdt_cen_stat");
			m_syst_unc->Write("m_dsdt_rel_syst");
		}
	}

	return 0;
}
