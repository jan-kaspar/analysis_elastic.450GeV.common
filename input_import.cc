#include <TDirectory.h>
#include "classes/command_line_tools.hh"

#include "TFile.h"

#include <cstdio>
#include <cstring>
#include <sstream>
#include <vector>

using namespace std;

//----------------------------------------------------------------------------------------------------

void PrintUsage()
{
	printf("USAGE: program <option> <option>\n");
	printf("OPTIONS:\n");
	printf("    -dir-low-beta <dir>       input directory with low-beta data\n");
	printf("    -dir-high-beta <dir>      input directory with high-beta data\n");
	printf("    -binnings <string>        comma-separated list of binnings\n");
}

//----------------------------------------------------------------------------------------------------

int main(int argc, const char **argv)
{
	// defaults
	string dir_low_beta = "";
	string dir_high_beta = "";
	string binnings_list = "";

	// parse command line
	for (int argi = 1; (argi < argc) && (cl_error == 0); ++argi)
	{
		if (strcmp(argv[argi], "-h") == 0 || strcmp(argv[argi], "--help") == 0)
		{
			cl_error = 1;
			continue;
		}

		if (TestStringParameter(argc, argv, argi, "-dir-low-beta", dir_low_beta)) continue;
		if (TestStringParameter(argc, argv, argi, "-dir-high-beta", dir_high_beta)) continue;
		if (TestStringParameter(argc, argv, argi, "-binnings", binnings_list)) continue;

		printf("ERROR: unknown option '%s'.\n", argv[argi]);
		cl_error = 1;
	}

	if (cl_error)
	{
		PrintUsage();
		return 1;
	}

	// validate input
	if (dir_low_beta.empty() || dir_high_beta.empty() || binnings_list.empty())
	{
		printf("ERROR: insufficient input\n");
		return 2;
	}

	// parse list of binnings
	vector<string> binnings;

	{
		stringstream ss(binnings_list);
		string binning;
		while (getline(ss, binning, ','))
			binnings.push_back(binning);
	}

	// open files
	TFile *f_hist_low_beta = TFile::Open((dir_low_beta + "/data/merged.root").c_str());
	TFile *f_unc_low_beta = TFile::Open((dir_low_beta + "/studies/systematics/matrix.root").c_str());

	TFile *f_hist_high_beta = TFile::Open((dir_high_beta + "/data/merged.root").c_str());
	TFile *f_unc_high_beta = TFile::Open((dir_high_beta + "/studies/systematics/matrix.root").c_str());

	TFile *f_out = TFile::Open("import.root", "recreate");

	if (!f_hist_low_beta || !f_unc_low_beta || !f_hist_high_beta || !f_unc_high_beta || !f_out)
	{
		printf("ERROR: some files cannot be opened.\n");
		return 3;
	}

	for (const auto &binning : binnings)
	{
		printf("* %s\n", binning.c_str());

		TDirectory *d_binning = f_out->mkdir(binning.c_str());

		gDirectory = d_binning->mkdir("low_beta");
		f_hist_low_beta->Get((binning + "/merged/combined/h_dsdt").c_str())->Write("h_dsdt_cen_stat");
		f_unc_low_beta->Get(("matrices/all-but-norm/" + binning + "/cov_mat").c_str())->Write("m_dsdt_rel_syst");
		f_unc_low_beta->Get(("matrices/all-but-norm/" + binning + "/h_stddev").c_str())->Write("h_dsdt_rel_syst");

		gDirectory = d_binning->mkdir("high_beta");
		f_hist_high_beta->Get((binning + "/merged/combined/h_dsdt").c_str())->Write("h_dsdt_cen_stat");
		f_unc_high_beta->Get(("matrices/all-but-norm/" + binning + "/cov_mat").c_str())->Write("m_dsdt_rel_syst");
		f_unc_high_beta->Get(("matrices/all-but-norm/" + binning + "/h_stddev").c_str())->Write("h_dsdt_rel_syst");
	}

	// clean up
	delete f_hist_low_beta;
	delete f_unc_low_beta;

	delete f_hist_high_beta;
	delete f_unc_high_beta;

	delete f_out;

	return 0;
}
