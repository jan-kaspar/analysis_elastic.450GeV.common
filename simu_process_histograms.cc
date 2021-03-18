#include <TDirectory.h>
#include "classes/command_line_tools.hh"

#include "classes/Stat.hh"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

#include <cstdio>
#include <cstring>
#include <memory>
#include <map>
#include <sstream>

using namespace std;

//----------------------------------------------------------------------------------------------------

void PrintUsage()
{
	printf("USAGE: program <option> <option>\n");
	printf("OPTIONS:\n");
	printf("    -seed-min <integer>        first seed to process\n");
	printf("    -seed-max <integer>        last seed to process\n");
	printf("    -binnings <string>         comma-separated list of binnings");
}

//----------------------------------------------------------------------------------------------------

int main(int argc, const char **argv)
{
	// defaults
	unsigned int seed_min = 0;
	unsigned int seed_max = 0;

	string binnings_str;
	string datasets_str = "low_beta,high_beta";

	// parse command line
	for (int argi = 1; (argi < argc) && (cl_error == 0); ++argi)
	{
		if (strcmp(argv[argi], "-h") == 0 || strcmp(argv[argi], "--help") == 0)
		{
			cl_error = 1;
			continue;
		}

		if (TestUIntParameter(argc, argv, argi, "-seed-min", seed_min)) continue;
		if (TestUIntParameter(argc, argv, argi, "-seed-max", seed_max)) continue;

		if (TestStringParameter(argc, argv, argi, "-binnings", binnings_str)) continue;

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
	vector<string> datasets;
	{
		stringstream ss(binnings_str);
		string item;
		while (getline(ss, item, ','))
			binnings.push_back(item);

		stringstream ss2(datasets_str);
		while (getline(ss2, item, ','))
			datasets.push_back(item);
	}

	// prepare data structures
	struct Block
	{
		bool initialized = false;
		unique_ptr<Stat> st;
		unique_ptr<TH1D> h_ref;
	};

	map<string, map<string, Block>> data;

	// loop over seeds
	for (unsigned int seed = seed_min; seed <= seed_max; ++seed)
	{
		char buf[200];
		sprintf(buf, "seed_%i/histograms.root", seed);

		TFile *f_in = TFile::Open(buf);

		for (const auto &binning : binnings)
		{
			for (const auto &dataset : datasets)
			{
				TH1D *h = (TH1D *) f_in->Get((binning + "/" + dataset + "/h_dsdt_cen_stat").c_str());

				printf("    %p\n", h);

				auto &d = data[binning][dataset];

				if (!d.initialized)
				{
					d.initialized = true;

					d.st = make_unique<Stat>(h->GetNbinsX());
					d.h_ref = make_unique<TH1D>(*h);
					d.h_ref->SetDirectory(nullptr);
				}

				vector<double> bin_content(h->GetNbinsX());
				for (int bi = 1; bi <= h->GetNbinsX(); ++bi)
					bin_content[bi - 1] = h->GetBinContent(bi);

				d.st->Fill(bin_content);
			}
		}

		delete f_in;
	}

	// save output
	TFile *f_out = TFile::Open("process_histograms.root", "recreate");

	for (const auto &bit : data)
	{
		TDirectory *d_binning = f_out->mkdir(bit.first.c_str());

		for (const auto &dit : bit.second)
		{
			TDirectory *d_dataset = d_binning->mkdir(dit.first.c_str());
			gDirectory = d_dataset;

			auto &d = dit.second;

			d.h_ref->Reset();

			for (int bi = 1; bi <= d.h_ref->GetNbinsX(); ++bi)
			{
				const int i = bi - 1;
				const double mu = d.st->GetMean(i);
				const double si = d.st->GetStdDev(i);
				const double si_unc = d.st->GetStdDevUncGauss(i);

				d.h_ref->SetBinContent(bi, (mu > 0) ? si/mu : 0.);
				d.h_ref->SetBinError(bi, (mu > 0) ? si_unc/mu : 0.);
			}

			d.h_ref->Write("h_stddev_rel");

			int n_bins = d.h_ref->GetNbinsX();
			TH2D *h2_corr = new TH2D("h2_corr", ";bin idx;bin idx", n_bins, -0.5, double(n_bins) - 0.5, n_bins, -0.5, double(n_bins) - 0.5);

			for (int i = 0; i < n_bins; ++i)
			{
				for (int j = 0; j < n_bins; ++j)
					h2_corr->SetBinContent(i+1, j+1, d.st->GetCorrelation(i, j));
			}

			h2_corr->Write("h2_rel_corr");
		}
	}

	delete f_out;

	return 0;
}