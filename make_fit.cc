#include "classes/command_line_tools.hh"

#include "TFile.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "Fit/Fitter.h"
#include "TMinuitMinimizer.h"

#include <cstring>
#include <functional>
#include <vector>
#include <sstream>

using namespace std;

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

struct Model
{
	Model();

	void ResetCaches();
	void UpdateCaches();

	double Eval(double t);
};

Model::Model()
{
	//FIXME: implement
}

//----------------------------------------------------------------------------------------------------

void Model::ResetCaches()
{
	//FIXME: implement

}

//----------------------------------------------------------------------------------------------------

void Model::UpdateCaches()
{
	//FIXME: implement
}

//----------------------------------------------------------------------------------------------------

double Model::Eval(double t)
{
	//FIXME: implement

	return 0;
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

struct Bin
{
	double t_left, t_right, t_repr;
	double dsdt, dsdt_unc_stat;
};

//----------------------------------------------------------------------------------------------------

struct Dataset
{
	// vector of bin data
	vector<Bin> bins;

	// input matrix of relative systematic uncertainties
	TMatrixDSym m_unc_syst_rel;

	// inverted covariance matrix - for use in chi^2 calculation
	TMatrixD m_cov_inv;

	Dataset(const string &file, const string &directory);

	// use initial settings (bin representative points, covariance matrix)
	void UpdateInitial();

	// use settings (bin representative points, covariance matrix) based on model from previous iteration
	void Update(const Model &model);
};

//----------------------------------------------------------------------------------------------------

Dataset::Dataset(const string &file, const string &directory)
{
	//FIXME: implement

	// load h_dsdt_cen_stat, m_dsdt_rel_syst
}

//----------------------------------------------------------------------------------------------------

void Dataset::UpdateInitial()
{
	//FIXME: implement
}

//----------------------------------------------------------------------------------------------------

void Dataset::Update(const Model &model)
{
	//FIXME: implement
}

//----------------------------------------------------------------------------------------------------

struct Data
{
	vector<Dataset> datasets;

	// use initial settings (bin representative points, covariance matrix)
	void UpdateInitial()
	{
		for (auto &d : datasets)
			d.UpdateInitial();
	}

	// use settings (bin representative points, covariance matrix) based on model from previous iteration
	void Update(const Model &model)
	{
		for (auto &d : datasets)
			d.Update(model);
	}
};

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

struct Metric
{
	const Data &data;
	Model &model;

	Metric(const Data &d, Model &m) : data(d), model(m) {}

	void SetModelParameters(const double *par) const;

	double operator() (const double *par) const;
};

//----------------------------------------------------------------------------------------------------

void Metric::SetModelParameters(const double *par) const
{
	// FIXME: implement
}

//----------------------------------------------------------------------------------------------------

double Metric::operator() (const double *par) const
{
	// decode parameters
	SetModelParameters(par);

	const double eta = 1.;	// FIXME: assign real value, make it per dataset

	// loop over datasets
	double s2 = 0;
	for (const auto &d : data.datasets)
	{
		unsigned int dim = d.bins.size();

		// build vector of differences
		vector<double> diff(dim);
		for (int i = 0; i < dim; ++i)
			diff[i] = d.bins[i].dsdt - eta * model.Eval(d.bins[i].t_repr);

		// calculate sum of squares
		for (int i = 0; i < dim; ++i)
		{
			for (int j = 0; j < dim; ++j)
				s2 += diff[i] * d.m_cov_inv(i, j) * diff[j];
		}
	}

	return s2;
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

struct Minimization
{
	const Data &data;
	Model &model;
	const Metric &metric;

	Minimization(const Data &da, Model &mo, const Metric &me);

	void Minimize();

	void ResultToModel();
};

//----------------------------------------------------------------------------------------------------

Minimization::Minimization(const Data &da, Model &mo, const Metric &me) : data(da), model(mo), metric(me)
{
	// FIXME: implement
}

//----------------------------------------------------------------------------------------------------

void Minimization::Minimize()
{
	// FIXME: implement
}

//----------------------------------------------------------------------------------------------------

void Minimization::ResultToModel()
{
	// FIXME: implement
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

void PrintUsage()
{
	printf("USAGE: program <option> <option>\n");
	printf("OPTIONS:\n");
	//printf("    -cfg <file>       config file\n");
}

//----------------------------------------------------------------------------------------------------

int main(int argc, const char **argv)
{
	// defaults
	//string cfg_file = "config.py";

	string input_file = "";
	string input_datasets = "high_beta,low_beta";
	string input_binning = "sb1";

	unsigned int n_iterations = 3;

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

	// load data
	Data data;
	{
		stringstream ss(input_datasets);
		string item;
		while (getline(ss, item, ','))
			data.datasets.emplace_back(Dataset{input_file, input_binning+"/"+item});
	}

	// initializations
	Model model;

	Metric metric(data, model);

	Minimization minimization(data, model, metric);

	// run iterations
	for (unsigned int it = 0; it < n_iterations; ++it)
	{
		if (it == 0)
		{
			data.UpdateInitial();
			model.ResetCaches();
		} else {
			data.Update(model);
			model.UpdateCaches();
		}

		minimization.Minimize();

		minimization.ResultToModel();
	}

	// save results
	// TODO

	return 0;
}
