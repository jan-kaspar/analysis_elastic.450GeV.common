#include "classes/command_line_tools.hh"
#include "classes/Result.hh"
#include "classes/Stat.hh"

#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TProfile2D.h"

#include <cstdio>
#include <cstring>
#include <memory>
#include <string>

using namespace std;

//----------------------------------------------------------------------------------------------------

void PrintUsage()
{
	printf("USAGE: program <option> <option>\n");
	printf("OPTIONS:\n");

	printf("    -ref-file <string>         file with reference data\n");
	printf("    -ref-object <string>       object with reference data\n");

	printf("    -seed-min <integer>        first seed to process\n");
	printf("    -seed-max <integer>        last seed to process\n");

	printf("    -input-pattern <string>    pattern of input file paths\n");

	printf("    -output <string>           output file name\n");
}

//----------------------------------------------------------------------------------------------------

struct ReducedResult {
	double de_si_tot, de_B, de_rho;
};

//----------------------------------------------------------------------------------------------------

void MakeOneFit(const vector<ReducedResult> &res, double range_de_B, double range_de_si_tot)
{
	// make fit
	// convention: x = de_B, y = de_si_tot, z = de_rho

	double s_xx=0., s_xy=0., s_yy=0., s_zx=0., s_zy=0.;
	for (const auto &r : res)
	{
		s_xx += r.de_B * r.de_B;
		s_xy += r.de_B * r.de_si_tot;
		s_yy += r.de_si_tot * r.de_si_tot;
		s_zx += r.de_rho * r.de_B;
		s_zy += r.de_rho * r.de_si_tot;
	}

	const double det = s_xx*s_yy - s_xy*s_xy;
	const double al = ( s_yy*s_zx - s_xy*s_zy) / det;
	const double be = (-s_xy*s_zx + s_xx*s_zy) / det;

	printf("al = %.3E\n", al);
	printf("be = %.3E\n", be);

	// plot residuals
	TProfile2D *p2 = new TProfile2D("p2", ";de_si_tot;de_B", 25, -range_de_si_tot, +range_de_si_tot, 25, -range_de_B, +range_de_B, "s");
	for (const auto &r : res)
	{
		const double diff = r.de_rho - al * r.de_B - be * r.de_si_tot;
		p2->Fill(r.de_si_tot, r.de_B, diff);
	}
	p2->Write();

	TH2D *h2_rms = new TH2D("h2_rms", ";de_si_tot;de_B", 25, -range_de_si_tot, +range_de_si_tot, 25, -range_de_B, +range_de_B);
	for (int bi = 1; bi <= h2_rms->GetNbinsX(); ++bi)
	{
		for (int bj = 1; bj <= h2_rms->GetNbinsY(); ++bj)
			h2_rms->SetBinContent(bi, bj, p2->GetBinError(bi, bj));
	}
	h2_rms->Write();
}

//----------------------------------------------------------------------------------------------------

int main(int argc, const char **argv)
{
	// defaults
	string ref_file;
	string ref_object;

	unsigned int seed_min = 0;
	unsigned int seed_max = 0;

	string input_pattern = "seed_%u/fit.root";

	string output = "process_fits.root";

	// parse command line
	for (int argi = 1; (argi < argc) && (cl_error == 0); ++argi)
	{
		if (strcmp(argv[argi], "-h") == 0 || strcmp(argv[argi], "--help") == 0)
		{
			cl_error = 1;
			continue;
		}

		if (TestStringParameter(argc, argv, argi, "-ref-file", ref_file)) continue;
		if (TestStringParameter(argc, argv, argi, "-ref-object", ref_object)) continue;

		if (TestUIntParameter(argc, argv, argi, "-seed-min", seed_min)) continue;
		if (TestUIntParameter(argc, argv, argi, "-seed-max", seed_max)) continue;

		if (TestStringParameter(argc, argv, argi, "-input-pattern", input_pattern)) continue;

		if (TestStringParameter(argc, argv, argi, "-output", output)) continue;

		printf("ERROR: unknown option '%s'.\n", argv[argi]);
		cl_error = 1;
	}

	if (cl_error)
	{
		PrintUsage();
		return 1;
	}

	// list of parameters
	struct ParameterInfo
	{
		string tag;
		double range;

		shared_ptr<TH1D> h_diff;

		ParameterInfo(const string &_t, double _r) : tag(_t), range(_r), h_diff(new TH1D("", "", 40, -_r, +_r)) {}
	};

	vector<ParameterInfo> parameters = {
		{"eta", 0.3},
		{"eta_unc", 0.3},
		{"A", 15},
		{"A_unc", 15},
		{"b1", 4},
		{"b1_unc", 4},
		{"B", 8},
		{"B_unc", 8},
		{"p0", 0.1},
		{"p0_unc", 0.1},
		{"rho", 0.1},
		{"rho_unc", 0.1},
		{"si_tot", 3},
		{"si_tot_unc", 3},
	};

	// prepare data
	Stat st(parameters.size());

	vector<vector<TH2D*>> correlation_plots;
	for (unsigned int i = 0; i < parameters.size(); ++i)
	{
		vector<TH2D*> v;
		for (unsigned int j = 0; j < parameters.size(); ++j)
			v.push_back(new TH2D("", "", 25, -parameters[i].range, +parameters[i].range, 25, -parameters[j].range, +parameters[j].range));
		correlation_plots.push_back(move(v));
	}

	vector<ReducedResult> reducedResults;

	// load reference
	unique_ptr<TFile> f_ref(new TFile(ref_file.c_str()));
	TH1D *h_ref = (TH1D *) f_ref->Get(ref_object.c_str());
	Result r_ref(h_ref);

	// loop over seeds
	for (unsigned int seed = seed_min; seed <= seed_max; ++seed)
	{
		char buf[200];
		sprintf(buf, input_pattern.c_str(), seed);
		unique_ptr<TFile> f_in(new TFile(buf));

		TH1D *h_test = (TH1D *) f_in->Get("final/results");
		Result r_test(h_test);

		vector<double> diff(parameters.size());
		ReducedResult rr;

		for (unsigned int pi = 0; pi < parameters.size(); ++pi)
		{
			const auto &p = parameters[pi];

			bool p_is_unc = (p.tag.find("_unc") != string::npos);

			const double v_ref = r_ref.Get(p.tag);
			const double v_test = r_test.Get(p.tag);

			if (p_is_unc)
				diff[pi] = v_test;
			else
				diff[pi] = v_test - v_ref;

			if (p.tag == "B") rr.de_B = diff[pi];
			if (p.tag == "rho") rr.de_rho = diff[pi];
			if (p.tag == "si_tot") rr.de_si_tot = diff[pi];

			p.h_diff->Fill(diff[pi]);
		}

		reducedResults.push_back(rr);

		st.Fill(diff);

		for (unsigned int i = 0; i < parameters.size(); ++i)
		{
			for (unsigned int j = 0; j < parameters.size(); ++j)
				correlation_plots[i][j]->Fill(diff[i], diff[j]);
		}
	}

	// save results
	unique_ptr<TFile> f_out(new TFile(output.c_str(), "recreate"));

	for (unsigned int pi = 0; pi < parameters.size(); ++pi)
	{
		const auto &p = parameters[pi];

		gDirectory = f_out->mkdir(p.tag.c_str());

		p.h_diff->Write("h_diff");

		unique_ptr<TGraph> g_data(new TGraph());
		g_data->SetPoint(0, st.GetMean(pi), st.GetMeanUnc(pi));
		g_data->SetPoint(1, st.GetStdDev(pi), st.GetStdDevUncGauss(pi));
		g_data->Write("g_data");
	}

	gDirectory = f_out->mkdir("correlations");

	for (unsigned int i = 0; i < parameters.size(); ++i)
	{
		for (unsigned int j = 0; j < parameters.size(); ++j)
		{
			TH2D *h2 = correlation_plots[i][j];
			h2->SetName((parameters[j].tag + " vs " + parameters[i].tag).c_str());
			h2->SetTitle((";" + parameters[i].tag + ";" + parameters[j].tag).c_str());
			h2->Write();
		}
	}

	gDirectory = f_out->mkdir("fit");

	double r_de_B = 1., r_de_tot = 1.;
	for (const auto &p : parameters)
	{
		if (p.tag == "B") r_de_B = p.range;
		if (p.tag == "si_tot") r_de_tot = p.range;
	}

	MakeOneFit(reducedResults, r_de_B, r_de_tot);

	return 0;
}
