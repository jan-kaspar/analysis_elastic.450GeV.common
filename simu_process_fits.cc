#include "classes/command_line_tools.hh"
#include "classes/Result.hh"
#include "classes/Stat.hh"

#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TGraph.h"

#include <cstring>
#include <memory>

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

		shared_ptr<TH1D> h_diff;

		ParameterInfo(const string &_t) : tag(_t), h_diff(new TH1D("", "", 10, 0., 0.)) {}
	};

	vector<ParameterInfo> parameters = {
		{"eta"},
		{"eta_unc"},
		{"A"},
		{"A_unc"},
		{"b1"},
		{"b1_unc"},
		{"p0"},
		{"p0_unc"},
		{"rho"},
		{"rho_unc"},
		{"si_tot"},
		{"si_tot_unc"},
	};

	// prepare data
	Stat st(parameters.size());

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
		for (unsigned int pi = 0; pi < parameters.size(); ++pi)
		{
			const auto &p = parameters[pi];

			const double v_ref = r_ref.Get(p.tag);
			const double v_test = r_test.Get(p.tag);

			diff[pi] = v_test - v_ref;

			p.h_diff->Fill(diff[pi]);
		}

		st.Fill(diff);
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
		g_data->SetPoint(1, st.GetStdDev(pi), st.GetStdDevUnc(pi));
		g_data->Write("g_data");
	}

	return 0;
}
