#include "classes/command_line_tools.hh"
#include "classes/HadronicFitModel.hh"
#include "classes/Result.hh"

#include "Elegent/Constants.h"
#include "Elegent/CoulombInterference.h"
#include "Elegent/Model.h"

#include "TFile.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "Fit/Fitter.h"
#include "TMinuitMinimizer.h"
#include "TSpline.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TMath.h"

#include <cstdio>
#include <cstring>
#include <functional>
#include <memory>
#include <vector>
#include <sstream>

using namespace std;

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

struct Model
{
	unsigned int n_b;

	unsigned int n_fit_parameters;

	unique_ptr<HadronicFitModel> hfm;

	bool caches_initialized;
	unique_ptr<TGraph> g_psi_re, g_psi_im;
	unique_ptr<TSpline3> s_psi_re, s_psi_im;

	enum Component { cC, cH, cCH };

	Model(unsigned int _n_b);

	void ResetCaches();
	void UpdateCaches();
	void WriteCaches() const;

	double Eval(double t, Component c = cCH) const;
};

//----------------------------------------------------------------------------------------------------

Model::Model(unsigned int _n_b) :
	n_b(_n_b),
	n_fit_parameters(n_b + 3),
	hfm(new HadronicFitModel()),
	caches_initialized(false)
{
	// initialise hadronic model
	hfm->modulusMode = HadronicFitModel::mmExp;
	hfm->phaseMode = HadronicFitModel::pmConstant;
	hfm->t1 = 0.5;
	hfm->t2 = 1.5;
	hfm->Print();

	// initialise Elegent
	using namespace Elegent;

	Constants::Init(2*450., cnts->mPP);
    cnts->Print();

	model = hfm.get();

	coulomb->mode = CoulombInterference::mKL;
	coulomb->ffType = coulomb->ffPuckett;
	coulomb->precision = 1E-3;
	coulomb->Print();
}

//----------------------------------------------------------------------------------------------------

void Model::ResetCaches()
{
	g_psi_re.reset();
	g_psi_im.reset();

	s_psi_re.reset();
	s_psi_im.reset();

	caches_initialized = false;
}

//----------------------------------------------------------------------------------------------------

void Model::UpdateCaches()
{
	g_psi_re = make_unique<TGraph>();
	g_psi_im = make_unique<TGraph>();

	double dmt = 1.;
	for (double mt = 1E-4; mt < 0.1; mt += dmt)
	{
		TComplex Psi = Elegent::coulomb->Psi_KL(-mt);

		int idx = g_psi_re->GetN();
		g_psi_re->SetPoint(idx, mt, Psi.Re());
		g_psi_im->SetPoint(idx, mt, Psi.Im());

		dmt = 0.0001;

		if (mt > 0.003)
			dmt = 0.001;

		if (mt > 0.02)
			dmt = 0.004;
	}

	s_psi_re = make_unique<TSpline3>("splinePsiRe", g_psi_re->GetX(), g_psi_re->GetY(), g_psi_re->GetN());
	s_psi_im = make_unique<TSpline3>("splinePsiRe", g_psi_im->GetX(), g_psi_im->GetY(), g_psi_im->GetN());

	caches_initialized = true;
}

//----------------------------------------------------------------------------------------------------

void Model::WriteCaches() const
{
	if (caches_initialized)
	{
		g_psi_re->Write("g_psi_re");
		g_psi_im->Write("g_psi_im");
	}
}

//----------------------------------------------------------------------------------------------------

double Model::Eval(double mt, Model::Component component) const
{
	using namespace Elegent;

	// amplitude components
	const TComplex F_C = coulomb->Amp_pure(-mt);
	const TComplex F_H = hfm->Amp(-mt);

	// amplitude choice
	TComplex F_T = 0.;

	if (component == cC)
		F_T = F_C;

	if (component == cH)
		F_T = F_H;

	if (component == cCH)
	{
		const TComplex Psi = (caches_initialized) ? TComplex(s_psi_re->Eval(mt), s_psi_im->Eval(mt)) : - coulomb->Phi_SWY(-mt);
		F_T = F_C + F_H * TComplex::Exp(i*Psi);
	}

	/*
	printf("mt=%.1E | FC: re=%+.1E, im=%+.1E, FH: re=%+.1E, im=%+.1E, Psi: re=%+.1E, im=%+.1E, FT: re=%+.1E, im=%+.1E | ds/dt = %.1E\n",
		mt, F_C.Re(), F_C.Im(), F_H.Re(), F_H.Im(), Psi.Re(), Psi.Im(), F_T.Re(), F_T.Im(), cnts->sig_fac * F_T.Rho2());
	*/

	return cnts->sig_fac * F_T.Rho2();

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
	// label
	string tag;

	// vector of bin data
	vector<Bin> bins;

	// input matrix of relative systematic uncertainties
	TMatrixDSym m_unc_syst_rel;

	// (inverted) covariance matrix - for use in chi^2 calculation
	TMatrixD m_cov, m_cov_inv;

	// normalisation factor
	double eta;

	Dataset(const string &_tag, const string &file, const string &directory, double t_min, double t_max);

	void Write() const;

	// update mutable elements (bin representative points, inverted covariance matrix)
	// if model != 0, then based on model from previous iteration
	void Update(const Model *model, bool use_stat_unc, bool use_syst_unc);
};

//----------------------------------------------------------------------------------------------------

Dataset::Dataset(const string &_tag, const string &file, const string &directory, double t_min, double t_max) :
	tag(_tag)
{
	unique_ptr<TFile> f_in(TFile::Open(file.c_str()));

	if (!f_in)
	{
		printf("ERROR: cannot open file '%s'\n", file.c_str());
		throw 1;
	}

	TH1D *h_dsdt_cen_stat = (TH1D *) f_in->Get((directory + "/h_dsdt_cen_stat").c_str());
	TMatrixD *m_dsdt_rel_syst_in = (TMatrixD *) f_in->Get((directory + "/m_dsdt_rel_syst").c_str());

	if (!h_dsdt_cen_stat || !m_dsdt_rel_syst_in)
	{
		printf("ERROR: cannot load input from file '%s', directory '%s'.\n", file.c_str(), directory.c_str());
		throw 1;
	}

	// define bin range
	const int bi_min = h_dsdt_cen_stat->FindBin(t_min);
	const int bi_max = h_dsdt_cen_stat->FindBin(t_max);
	const int dim = bi_max - bi_min + 1;

	// extract bin data
	bins.resize(dim);
	for (int bi = bi_min; bi <= bi_max; ++bi)
	{
		const double l = h_dsdt_cen_stat->GetBinLowEdge(bi);
		const double r = l + h_dsdt_cen_stat->GetBinWidth(bi);

		const double v = h_dsdt_cen_stat->GetBinContent(bi);
		const double u = h_dsdt_cen_stat->GetBinError(bi);

		int idx = bi - bi_min;
		bins[idx] = Bin{l, r, 0., v, u};
	}

	// extract systematics
	m_unc_syst_rel.ResizeTo(dim, dim);
	for (int i = 0; i < dim; ++i)
	{
		for (int j = 0; j < dim; ++j)
		{
			const int offset = bi_min - 1;
			m_unc_syst_rel(i, j) = (*m_dsdt_rel_syst_in)(i + offset, j + offset);
		}
	}
}

//----------------------------------------------------------------------------------------------------

void Dataset::Write() const
{
	unique_ptr<TGraphErrors> g_data_t_unc(new TGraphErrors), g_data_dsdt_unc_stat(new TGraphErrors),
		g_data_dsdt_unc_syst(new TGraphErrors), g_data_dsdt_unc_comb(new TGraphErrors);

	for (unsigned int i = 0; i < bins.size(); ++i)
	{
		const auto &b = bins[i];

		const double t_cen = (b.t_right + b.t_left) / 2.;
		const double t_unc = (b.t_right - b.t_left) / 2.;
		const double t_repr = b.t_repr;

		const double dsdt = b.dsdt / eta;
		const double dsdt_unc_stat = b.dsdt_unc_stat / eta;
		const double dsdt_unc_syst = sqrt(m_unc_syst_rel(i, i)) * dsdt / eta;
		const double dsdt_unc_comb = sqrt(m_cov(i, i)) / eta;

		int idx = g_data_t_unc->GetN();

		g_data_t_unc->SetPoint(idx, t_cen, dsdt);
		g_data_t_unc->SetPointError(idx, t_unc, 0.);

		g_data_dsdt_unc_stat->SetPoint(idx, t_repr, dsdt);
		g_data_dsdt_unc_stat->SetPointError(idx, 0., dsdt_unc_stat);

		g_data_dsdt_unc_syst->SetPoint(idx, t_repr, dsdt);
		g_data_dsdt_unc_syst->SetPointError(idx, 0., dsdt_unc_syst);

		g_data_dsdt_unc_comb->SetPoint(idx, t_repr, dsdt);
		g_data_dsdt_unc_comb->SetPointError(idx, 0., dsdt_unc_comb);
	}

	g_data_t_unc->Write("g_data_t_unc");
	g_data_dsdt_unc_stat->Write("g_data_dsdt_unc_stat");
	g_data_dsdt_unc_syst->Write("g_data_dsdt_unc_syst");
	g_data_dsdt_unc_comb->Write("g_data_dsdt_unc_comb");
}

//----------------------------------------------------------------------------------------------------

void Dataset::Update(const Model *model, bool use_stat_unc, bool use_syst_unc)
{
	// adjust bin representative points
	for (auto &bin : bins)
	{
		if (model == nullptr)
		{
			bin.t_repr = (bin.t_left + bin.t_right) / 2.;
		} else {
			double l = bin.t_left;
			double r = bin.t_right;
			const double w = r - l;

			// calculate integral
			unsigned int n_div = 100;
			const double wo = w / n_div;
			double I = 0.;
			for (unsigned int i = 0; i < n_div; ++i)
				I += model->Eval(l + (0.5 + double(i)) * wo);
			I *= wo;

			// find representative point
			double xr;
			while (r - l > w/10000)
			{
				xr = (r+l)/2.;
				
				if (model->Eval(xr) < I/w)
					r -= (r-l)/2.;
				else
					l += (r-l)/2.;
			}
			bin.t_repr = (r+l)/2.;
		}
	}

	// build uncertainty matrix
	const int dim = bins.size();

	m_cov.ResizeTo(dim, dim);

	for (int i = 0; i < dim; ++i)
	{
		for (int j = 0; j < dim; ++j)
		{
			// reference ds/dt to make the conversion relative -> absolute
			double y_ref_i = 0., y_ref_j = 0.;
			if (!model)
			{
				y_ref_i = bins[i].dsdt;
				y_ref_j = bins[j].dsdt;
			} else {
				y_ref_i = model->Eval(bins[i].t_repr);
				y_ref_j = model->Eval(bins[j].t_repr);
			}

			double &e = m_cov(i, j);

			e = 0;

			if (use_stat_unc && i == j)
				e += pow(bins[i].dsdt_unc_stat, 2);

			if (use_syst_unc)
				e += y_ref_i * m_unc_syst_rel(i, j) * y_ref_j;
		}
	}

	// invert covariance matrix
	m_cov_inv.ResizeTo(m_cov);
	m_cov_inv = m_cov;
	m_cov_inv.Invert();
}

//----------------------------------------------------------------------------------------------------

struct Data
{
	vector<Dataset> datasets;

	unsigned int NPoints() const
	{
		unsigned int n = 0;
		for (const auto &ds : datasets)
			n += ds.bins.size();
		return n;
	}

	// use initial settings (bin representative points, covariance matrix)
	void UpdateInitial(bool use_stat_unc, bool use_syst_unc)
	{
		for (auto &d : datasets)
			d.Update(nullptr, use_stat_unc, use_syst_unc);
	}

	// use settings (bin representative points, covariance matrix) based on model from previous iteration
	void Update(const Model &model, bool use_stat_unc, bool use_syst_unc)
	{
		for (auto &d : datasets)
			d.Update(&model, use_stat_unc, use_syst_unc);
	}
};

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

struct Metric
{
	Data &data;
	Model &model;

	Metric(Data &d, Model &m) : data(d), model(m) {}

	void SetModelParameters(const double *par);

	double CalculateChi2(bool debug) const;

	double operator() (const double *par);
};

//----------------------------------------------------------------------------------------------------

void Metric::SetModelParameters(const double *par)
{
	// order of parameters:
	//    eta
	//    a
	//    b_1, b_2, ...
	//    p0, ...

	for (auto &ds : data.datasets)
		ds.eta = par[0];

	model.hfm->a = par[1] * 1E6;

	for (unsigned int bi = 0; bi < model.n_b; ++bi)
	{
		if (bi == 0) model.hfm->b1 = par[2 + bi];
		if (bi == 1) model.hfm->b2 = par[2 + bi];
		if (bi == 2) model.hfm->b3 = par[2 + bi];
		if (bi == 3) model.hfm->b4 = par[2 + bi];
	}

	model.hfm->p0 = par[2 + model.n_b];
}

//----------------------------------------------------------------------------------------------------

double Metric::CalculateChi2(bool debug) const
{
	// loop over datasets
	double s2 = 0;
	for (const auto &d : data.datasets)
	{
		unsigned int dim = d.bins.size();

		// build vector of differences
		vector<double> diff(dim);
		for (unsigned int i = 0; i < dim; ++i)
			diff[i] = d.bins[i].dsdt - d.eta * model.Eval(d.bins[i].t_repr);

		// calculate sum of squares
		for (unsigned int i = 0; i < dim; ++i)
		{
			for (unsigned int j = 0; j < dim; ++j)
				s2 += diff[i] * d.m_cov_inv(i, j) * diff[j];
		}

		if (debug)
		{
			for (unsigned int i = 0; i < dim; ++i)
			{
				const double c = diff[i] * d.m_cov_inv(i, i) * diff[i];
				printf("    %s, bin %i: diff = %.3f, c = %.2E\n", d.tag.c_str(), i, diff[i], c);
			}
		}
	}

	return s2;
}

//----------------------------------------------------------------------------------------------------

double Metric::operator() (const double *par)
{
	SetModelParameters(par);

	return CalculateChi2(false);
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

struct InitialSettings
{
	bool eta_fixed;
	double eta;

	bool A_fixed;
	double A;

	bool b1_fixed;
	double b1;

	bool rho_fixed;
	double rho;
};

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

struct Minimization
{
	const Data &data;
	Model &model;
	Metric &metric;

	ROOT::Fit::Fitter fitter;

	Minimization(const Data &da, Model &mo, Metric &me, const InitialSettings &is);

	void Minimize();

	void PrintState() const;

	void PrintResults() const;

	void ResultToModel();

	void WriteGraphs() const;

	Result GetResults() const;
};

//----------------------------------------------------------------------------------------------------

Minimization::Minimization(const Data &da, Model &mo, Metric &me, const InitialSettings &is) : data(da), model(mo), metric(me)
{
	// initialize fitter
	double pStart[model.n_fit_parameters];
	fitter.SetFCN(model.n_fit_parameters, metric, pStart, 0, true);

	// initialize parameters
	{
		auto &cp = fitter.Config().ParSettings(0);
		cp.Set("eta", is.eta, 0.05, 0.70, 1.30);
		if (is.eta_fixed)
			cp.Fix();
	}
	// FIXME: handle fixed variables

	{
		const double a = sqrt(is.A / Elegent::cnts->sig_fac);
		auto &cp = fitter.Config().ParSettings(1);
		cp.Set("a", a / 1E6, a / 1E6 / 10);
		if (is.A_fixed)
			cp.Fix();
	}

	for (unsigned int i = 0; i < model.n_b; ++i)
	{
		char buf[100];
		sprintf(buf, "b%u", i+1);

		double v = 0., u = 0.;
		if (i == 0) { v = is.b1; u = v / 10.; }

		if (u <= 0.)
			u = 1.;

		auto &cp = fitter.Config().ParSettings(2 + i);
		cp.Set(buf, v, u);

		if (i == 0 && is.b1_fixed)
			cp.Fix();
	}

	{
		const double p0 = M_PI/2. - atan(is.rho);
		auto &cp = fitter.Config().ParSettings(model.n_b + 2);
		cp.Set("p0", p0, 0.01);
		if (is.rho_fixed)
			cp.Fix();
	}

	// set parameters to model
	double par[model.n_fit_parameters];
	for (unsigned int i = 0; i < model.n_fit_parameters; ++i)
		par[i] = fitter.Config().ParamsValues()[i];

	metric.SetModelParameters(par);
}

//----------------------------------------------------------------------------------------------------

Result Minimization::GetResults() const
{
	using namespace Elegent;

	Result r;
	
	const ROOT::Fit::FitResult &fr = fitter.Result();

	unsigned int ndf = data.NPoints() - model.n_fit_parameters;

	r.Set("chi2", fr.Chi2());
	r.Set("ndf", ndf);
	r.Set("chi2_norm", fr.Chi2() / ndf);
	r.Set("prob", TMath::Prob(fr.Chi2(), ndf));

	unsigned int idx;

	idx = 0;
	r.Set("eta", fr.Parameter(idx));
	r.Set("eta_unc", sqrt(fr.CovMatrix(idx, idx)));

	idx = 1;
	const double a = fr.Parameter(idx) * 1E6;
	const double a_unc = sqrt(fr.CovMatrix(idx, idx)) * 1E6;
	const double A = cnts->sig_fac * a*a;
	r.Set("A", A);
	r.Set("A_unc", cnts->sig_fac * 2*a * a_unc);

	idx = 2;
	r.Set("b1", fr.Parameter(idx));
	r.Set("b1_unc", sqrt(fr.CovMatrix(idx, idx)));
	r.Set("B", 2.*fr.Parameter(idx));
	r.Set("B_unc", 2.*sqrt(fr.CovMatrix(idx, idx)));

	idx = 2 + model.n_b;
	const double p0 = fr.Parameter(idx);
	const double p0_unc = sqrt(fr.CovMatrix(idx, idx));
	r.Set("p0", p0);
	r.Set("p0_unc", p0_unc);

	const double rho = cos(p0) / sin(p0);
	const double rho_unc = fabs(1. / sin(p0) / sin(p0)) * p0_unc;

	r.Set("rho", rho);
	r.Set("rho_unc", rho_unc);

	const double si_tot = sqrt( 16.*cnts->pi * cnts->sq_hbarc / (1. + rho * rho) * A );

	const double V_a_a = fr.CovMatrix(1, 1);
	const double V_a_p0 = fr.CovMatrix(1, idx);
	const double V_p0_p0 =  fr.CovMatrix(idx, idx);

	const double sc_A = 2. * cnts->sig_fac * fr.Parameter(1) * 1E6 * 1E6;
	const double sc_rho = 1. / sin(p0) / sin(p0);

	const double V_A_A = sc_A * V_a_a * sc_A;
	const double V_A_rho = sc_A * V_a_p0 * sc_rho;
	const double V_rho_rho = sc_rho * V_p0_p0 * sc_rho;

	const double der_A = si_tot/2. * 1./A;
	const double der_rho = si_tot/2. * 2.*rho / (1. + rho*rho);

	const double V_si_tot = der_A * V_A_A * der_A + 2.* der_A * V_A_rho * der_rho + der_rho * V_rho_rho * der_rho;
	const double si_tot_unc = sqrt(V_si_tot);

	r.Set("si_tot", si_tot);
	r.Set("si_tot_unc", si_tot_unc);

	return r;
}

//----------------------------------------------------------------------------------------------------

void Minimization::Minimize()
{
	fitter.FitFCN();
	fitter.FitFCN();
}

//----------------------------------------------------------------------------------------------------

void Minimization::PrintState() const
{
	for (unsigned int i = 0; i < fitter.Config().NPar(); ++i)
	{
		const auto &pd = fitter.Config().ParSettings(i);
		printf("idx %i [%3s]: val=%.3f, step=%.3f, min=%.3f, max=%.3f\n", i, pd.Name().c_str(),
			pd.Value(), pd.StepSize(), pd.LowerLimit(), pd.UpperLimit());
	}
}

//----------------------------------------------------------------------------------------------------

void Minimization::PrintResults() const
{
	const ROOT::Fit::FitResult &result = fitter.Result();

	unsigned int ndf = data.NPoints() - model.n_fit_parameters;

	printf("chi^2 = %.3f, ndf = %u, chi^2/ndf = %.3f\n", result.Chi2(), ndf, result.Chi2() / ndf);

	for (unsigned int i = 0; i < result.NPar(); ++i)
		printf("idx %u [%3s]: %+.3E +- %.3E\n", i, result.ParName(i).c_str(), result.Parameter(i), sqrt(result.CovMatrix(i, i)));
}

//----------------------------------------------------------------------------------------------------

void Minimization::ResultToModel()
{
	double par[model.n_fit_parameters];
	const ROOT::Fit::FitResult &result = fitter.Result();
	for (unsigned int i = 0; i < result.NPar(); ++i)
		par[i] = result.Parameter(i);

	metric.SetModelParameters(par);
}

//----------------------------------------------------------------------------------------------------

void Minimization::WriteGraphs() const
{
	TDirectory *d_top = gDirectory;

	for (const auto &ds : data.datasets)
	{
		gDirectory = d_top->mkdir(ds.tag.c_str());
		ds.Write();
	}

	gDirectory = d_top;

	unique_ptr<TGraph> g_fit_C(new TGraph), g_fit_H(new TGraph), g_fit_CH(new TGraph), g_fit_Z(new TGraph);

	g_fit_C->SetName("g_fit_C");
	g_fit_C->SetLineColor(4);

	g_fit_H->SetName("g_fit_H");
	g_fit_H->SetLineColor(8);

	g_fit_CH->SetName("g_fit_CH");
	g_fit_CH->SetLineColor(2);

	g_fit_Z->SetName("g_fit_Z");

	for (double t = 0.5E-4; t <= 0.1; )
	{
		const double v_C = model.Eval(t, Model::cC);
		const double v_H = model.Eval(t, Model::cH);
		const double v_CH = model.Eval(t, Model::cCH);

		int idx = g_fit_C->GetN();

		g_fit_C->SetPoint(idx, t, v_C);
		g_fit_H->SetPoint(idx, t, v_H);
		g_fit_CH->SetPoint(idx, t, v_CH);

		g_fit_Z->SetPoint(idx, t, (v_CH - v_C - v_H) / v_CH);

		double dt = 0.5E-5;
		if (t > 0.0008)
			dt = 0.5E-4;
		if (t > 0.006)
			dt = 0.002;
		if (t > 0.1)
			dt = 0.02;
		t += dt;
	}

	g_fit_C->Write();
	g_fit_H->Write();
	g_fit_CH->Write();
	g_fit_Z->Write();

	unique_ptr<TCanvas> c_fit_cmp(new TCanvas("c_fit_cmp"));
	g_fit_CH->Draw("al");
	g_fit_C->Draw("l");
	g_fit_H->Draw("l");
	c_fit_cmp->Write();
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

void PrintUsage()
{
	printf("USAGE: program <option> <option>\n");
	printf("OPTIONS:\n");

	printf("    -input-file <string>        input file\n");
	printf("    -input-datasets <string>    comma-separated list of datasets\n");
	printf("    -input-binning <string>     binning to use from input file\n");

	printf("    -t-min-low-beta <double>    lowest |t| to consider from low-beta dataset\n");
	printf("    -t-max-low-beta <double>    higher |t| to consider from low-beta dataset\n");
	printf("    -t-min-high-beta <double>   lowest |t| to consider from high-beta dataset\n");
	printf("    -t-max-high-beta <double>   higher |t| to consider from high-beta dataset\n");

	printf("    -use-stat-unc <bool>        whether to include statistical unc. to fit matrix\n");
	printf("    -use-syst-unc <bool>        whether to include statistical unc. to fit matrix\n");

	printf("    -n-b <int>                  number of b parameters\n");

	printf("    -use-eta-fixed <bool>       if eta shall be fixed in minimisation\n");
	printf("    -init-eta <double>          initial value of eta\n");

	printf("    -use-A-fixed <bool>         if A shall be fixed in minimisation\n");
	printf("    -init-A <double>            initial value of A\n");

	printf("    -use-b1-fixed <bool>        if b1 shall be fixed in minimisation\n");
	printf("    -init-b1 <double>           initial value of b1\n");

	printf("    -use-rho-fixed <bool>       if rho shall be fixed in minimisation\n");
	printf("    -init-rho <double>          initial value of rho\n");

	printf("    -n-iterations <integer>             number of fit iterations\n");
	printf("    -use-safe-first-iteration <bool>    don't use init model for adjustments\n");

	printf("    -output-root <string>               output ROOT file\n");
	printf("    -output-results <string>            output result file\n");
}

//----------------------------------------------------------------------------------------------------

int main(int argc, const char **argv)
{
	// defaults
	string input_file = "";
	string input_datasets = "high_beta,low_beta";
	string input_binning = "sb1";

	double t_min_low_beta = 0.017;
	double t_max_low_beta = 0.053;
	double t_min_high_beta = 3.1E-4;
	double t_max_high_beta = 0.023;

	bool use_stat_unc = true;
	bool use_syst_unc = true;

	unsigned int n_b = 1;
	bool use_safe_first_iteration = false;

	InitialSettings is;
	is.eta_fixed = false;
	is.eta = 1;
	is.A_fixed = false;
	is.A = 239.8;
	is.b1_fixed = false;
	is.b1 = 8.5;
	is.rho_fixed = false;
	is.rho = 0.10;

	unsigned int n_iterations = 3;

	string output_root = "make_fit.root";
	string output_results = "make_fit.results";

	// parse command line
	for (int argi = 1; (argi < argc) && (cl_error == 0); ++argi)
	{
		if (strcmp(argv[argi], "-h") == 0 || strcmp(argv[argi], "--help") == 0)
		{
			cl_error = 1;
			continue;
		}

		if (TestStringParameter(argc, argv, argi, "-input-file", input_file)) continue;
		if (TestStringParameter(argc, argv, argi, "-input-datasets", input_datasets)) continue;
		if (TestStringParameter(argc, argv, argi, "-input-binning", input_binning)) continue;

		if (TestDoubleParameter(argc, argv, argi, "-t-min-low-beta", t_min_low_beta)) continue;
		if (TestDoubleParameter(argc, argv, argi, "-t-max-low-beta", t_max_low_beta)) continue;
		if (TestDoubleParameter(argc, argv, argi, "-t-min-high-beta", t_min_high_beta)) continue;
		if (TestDoubleParameter(argc, argv, argi, "-t-max-high-beta", t_max_high_beta)) continue;

		if (TestBoolParameter(argc, argv, argi, "-use-stat-unc", use_stat_unc)) continue;
		if (TestBoolParameter(argc, argv, argi, "-use-syst-unc", use_syst_unc)) continue;

		if (TestUIntParameter(argc, argv, argi, "-n-b", n_b)) continue;

		if (TestBoolParameter(argc, argv, argi, "-use-eta-fixed", is.eta_fixed)) continue;
		if (TestDoubleParameter(argc, argv, argi, "-init-eta", is.eta)) continue;

		if (TestBoolParameter(argc, argv, argi, "-use-A-fixed", is.A_fixed)) continue;
		if (TestDoubleParameter(argc, argv, argi, "-init-A", is.A)) continue;

		if (TestBoolParameter(argc, argv, argi, "-use-b1-fixed", is.b1_fixed)) continue;
		if (TestDoubleParameter(argc, argv, argi, "-init-b1", is.b1)) continue;

		if (TestBoolParameter(argc, argv, argi, "-use-rho-fixed", is.rho_fixed)) continue;
		if (TestDoubleParameter(argc, argv, argi, "-init-rho", is.rho)) continue;

		if (TestUIntParameter(argc, argv, argi, "-n-iterations", n_iterations)) continue;
		if (TestBoolParameter(argc, argv, argi, "-use-safe-first-iteration", use_safe_first_iteration)) continue;

		if (TestStringParameter(argc, argv, argi, "-output-root", output_root)) continue;
		if (TestStringParameter(argc, argv, argi, "-output-results", output_results)) continue;

		printf("ERROR: unknown option '%s'.\n", argv[argi]);
		cl_error = 1;
	}

	if (cl_error)
	{
		PrintUsage();
		return 1;
	}

	// print settings
	printf("* input settings:\n");
	printf("    input_file=%s\n", input_file.c_str());
	printf("    input_datasets=%s\n", input_datasets.c_str());
	printf("    input_binning=%s\n", input_binning.c_str());

	printf("    t_min_low_beta=%.4f\n", t_min_low_beta);
	printf("    t_max_low_beta=%.4f\n", t_max_low_beta);
	printf("    t_min_high_beta=%.4f\n", t_min_high_beta);
	printf("    t_max_high_beta=%.4f\n", t_max_high_beta);

	printf("    use_stat_unc=%u\n", use_stat_unc);
	printf("    use_syst_unc=%u\n", use_syst_unc);

	printf("    n_b=%u\n", n_b);

	printf("    is.eta_fixed=%i\n", is.eta_fixed);
	printf("    is.eta=%.2f\n", is.eta);
	printf("    is.A_fixed=%i\n", is.A_fixed);
	printf("    is.A=%.2f\n", is.A);
	printf("    is.b1_fixed=%i\n", is.b1_fixed);
	printf("    is.b1=%.2f\n", is.b1);
	printf("    is.rho_fixed=%i\n", is.rho_fixed);
	printf("    is.rho=%.3f\n", is.rho);

	printf("    n_iterations=%u\n", n_iterations);
	printf("    use_safe_first_iteration=%u\n", use_safe_first_iteration);

	printf("    output_root=%s\n", output_root.c_str());
	printf("    output_results=%s\n", output_results.c_str());

	// open output file
	bool save_root = false;
	unique_ptr<TFile> f_root;
	if (!output_root.empty())
	{
		save_root = true;
		f_root = make_unique<TFile>(output_root.c_str(), "recreate");
	}

	// load data
	Data data;
	{
		stringstream ss(input_datasets);
		string item;
		while (getline(ss, item, ','))
		{
			double t_min = 0., t_max = 1.;
			if (item == "low_beta") { t_min = t_min_low_beta; t_max = t_max_low_beta; }
			if (item == "high_beta") { t_min = t_min_high_beta; t_max = t_max_high_beta; }
			data.datasets.emplace_back(Dataset{item, input_file, input_binning+"/"+item, t_min, t_max});
		}
	}

	if (data.datasets.empty())
	{
		printf("ERROR: no datasets given.\n");
		return 10;
	}
	
	printf("* datasets:\n");
	for (const auto &ds : data.datasets)
		printf("    %s, n_bins = %lu\n", ds.tag.c_str(), ds.bins.size());

	// initializations
	Model model(n_b);

	Metric metric(data, model);

	Minimization minimization(data, model, metric, is);

	// run iterations
	for (unsigned int it = 0; it < n_iterations; ++it)
	{
		printf("\n----- iteration %i -----\n", it);
	
		if (save_root)
		{
			char buf[100];
			sprintf(buf, "iteration %i", it);
			gDirectory = f_root->mkdir(buf);
		}

		if (it == 0 && use_safe_first_iteration)
		{
			data.UpdateInitial(use_stat_unc, use_syst_unc);
			model.ResetCaches();
		} else {
			data.Update(model, use_stat_unc, use_syst_unc);
			model.UpdateCaches();
		}

		printf("* before minimization:\n");
		minimization.PrintState();

		printf("check: chi^2 = %.3f\n", metric.CalculateChi2(false));

		minimization.Minimize();

		printf("* after minimization:\n");
		minimization.PrintResults();

		minimization.ResultToModel();

		printf("check: chi^2 = %.3f\n", metric.CalculateChi2(false));

		// save details/debug info
		if (save_root)
		{
			model.WriteCaches();
			minimization.WriteGraphs();
		}
	}
	
	printf("\n----- after iterations -----\n");
	
	// save results
	const auto &results = minimization.GetResults();

	if (save_root)
	{
		gDirectory = f_root->mkdir("final");
		results.ExportROOT()->Write("results");
	}

	if (!output_results.empty())
		results.Write(output_results);

	return 0;
}
