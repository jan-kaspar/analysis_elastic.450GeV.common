import root;
import pad_layout;
include "../common.asy";
//include "shared.asy";

string version = "default";

string types[] = {
	"both",
	"high",
};

string fit_iteration = "iteration 2";

string datasets[];
pen d_pens[];
datasets.push("high_beta"); d_pens.push(red + squarecap);
datasets.push("low_beta"); d_pens.push(blue + squarecap);

TGraph_errorBar = None;

mark mHL = mark((-1,0)--(+1,0)--cycle);
mark mVL = mark((0,-1)--(0,+1)--cycle);
mark mSY = mSq+false;

//----------------------------------------------------------------------------------------------------

struct Results
{
	real chi2, chi2_norm;
	real ndf;
	real prob;

	real eta, eta_unc;
	real A, A_unc;
	real b1, b1_unc;
	real p0, p0_unc;
	real rho, rho_unc;
	real si_tot, si_tot_unc;
};

//----------------------------------------------------------------------------------------------------

Results LoadResults(RootObject h)
{
	Results r;

	RootObject axis = h.oExec("GetXaxis");

	for (int bi = 1; bi <= h.iExec("GetNbinsX"); ++bi)
	{
		string l = axis.sExec("GetBinLabel", bi);
		
		real v = h.rExec("GetBinContent", bi);

		if (l == "chi2") r.chi2 = v;
		if (l == "ndf") r.ndf = v;
		if (l == "chi2_norm") r.chi2_norm = v;
		if (l == "prob") r.prob = v;

		if (l == "eta") r.eta = v;
		if (l == "eta_unc") r.eta_unc = v;

		if (l == "A") r.A = v;
		if (l == "A_unc") r.A_unc = v;

		if (l == "b1") r.b1 = v;
		if (l == "b1_unc") r.b1_unc = v;

		if (l == "p0") r.p0 = v;
		if (l == "p0_unc") r.p0_unc = v;

		if (l == "rho") r.rho = v;
		if (l == "rho_unc") r.rho_unc = v;

		if (l == "si_tot") r.si_tot = v;
		if (l == "si_tot_unc") r.si_tot_unc = v;
	}

	return r;
}

//----------------------------------------------------------------------------------------------------

void DrawSystematicMark(real l, real r, real y_cen, real y_unc, pen p)
{
	real e = 0.003;

	pair L = Scale((l, y_cen - y_unc));
	pair R = Scale((r, y_cen - y_unc));
	draw((L+(0, +e))--L--R--(R+(0, +e)), p);

	pair L = Scale((l, y_cen + y_unc));
	pair R = Scale((r, y_cen + y_unc));
	draw((L+(0, -e))--L--R--(R+(0, -e)), p);
}

//----------------------------------------------------------------------------------------------------

void DrawSystematics(RootObject g1, RootObject g2, pen p)
{
	for (int i = 0; i < g2.iExec("GetN"); ++i)
	{
		real ax[] = {0};
		real ay[] = {0};

		g2.vExec("GetPoint", i, ax, ay);
		real t = ax[0];
		real dsdt = ay[0];
		real dsdt_unc = g2.rExec("GetErrorY", i);

		real t_unc = g1.rExec("GetErrorX", i);

		draw(Scale((t-t_unc, dsdt-dsdt_unc))--Scale((t+t_unc, dsdt-dsdt_unc)), p);
		draw(Scale((t-t_unc, dsdt+dsdt_unc))--Scale((t+t_unc, dsdt+dsdt_unc)), p);
	}
}

//----------------------------------------------------------------------------------------------------

void DrawOne(string fit_file, real t_min, real t_max, real x_size)
{
	NewPad("$|t|\ung{GeV^2}$", "$\d\si/\d t\ung{mb/GeV}$", xSize = x_size);

	scale(Linear, Log);

	for (int dsi : datasets.keys)
	{
		string ds = datasets[dsi]; 
		pen p = d_pens[dsi];

		RootObject g_data_t_unc = RootGetObject(fit_file, fit_iteration + "/" + ds + "/g_data_t_unc", error=false);
		RootObject g_data_dsdt_unc_stat = RootGetObject(fit_file, fit_iteration + "/" + ds + "/g_data_dsdt_unc_stat", error=false);
		RootObject g_data_dsdt_unc_syst = RootGetObject(fit_file, fit_iteration + "/" + ds + "/g_data_dsdt_unc_syst", error=false);

		if (g_data_t_unc.valid)
		{
			draw(g_data_t_unc, "p", p);
			draw(g_data_dsdt_unc_stat, "p", p);
			DrawSystematics(g_data_t_unc, g_data_dsdt_unc_syst, p);

			AddToLegend(replace(datasets[dsi], "_", " "), p);
		}
	}

	RootObject g_fit_CH = RootGetObject(fit_file, fit_iteration + "/g_fit_CH");

	draw(g_fit_CH, "l", heavygreen, "fit C+H");

	limits((t_min, 5e1), (t_max, 3e3), Crop);
}

//----------------------------------------------------------------------------------------------------

void DrawOneRelative(string fit_file, real t_min, real t_max, real x_size)
{
	NewPad("$|t|\ung{GeV^2}$", "$(\d\si/\d t - \hbox{fit}) / \hbox{fit}$", xSize = x_size);

	RootObject g_fit_CH = RootGetObject(fit_file, fit_iteration + "/g_fit_CH");

	for (int dsi : datasets.keys)
	{
		string ds = datasets[dsi]; 
		pen p = d_pens[dsi];

		RootObject g_data_t_unc = RootGetObject(fit_file, fit_iteration + "/" + ds + "/g_data_t_unc", error=false);
		RootObject g_data_dsdt_unc_stat = RootGetObject(fit_file, fit_iteration + "/" + ds + "/g_data_dsdt_unc_stat", error=false);
		RootObject g_data_dsdt_unc_syst = RootGetObject(fit_file, fit_iteration + "/" + ds + "/g_data_dsdt_unc_syst", error=false);

		if (g_data_t_unc.valid)
		{
			for (int i = 0; i < g_data_t_unc.iExec("GetN"); ++i)
			{
				real ax[] = {0};
				real ay[] = {0};

				g_data_t_unc.vExec("GetPoint", i, ax, ay);
				real t_cen = ax[0];
				real dsdt = ay[0];
				real t_unc = g_data_t_unc.rExec("GetErrorX", i);

				g_data_dsdt_unc_stat.vExec("GetPoint", i, ax, ay);
				real t_repr = ax[0];
				real dsdt_unc_stat = g_data_dsdt_unc_stat.rExec("GetErrorY", i);

				real dsdt_unc_syst = g_data_dsdt_unc_syst.rExec("GetErrorY", i);

				real dsdt_ref = g_fit_CH.rExec("Eval", t_repr);

				real dsdt_rel = dsdt / dsdt_ref - 1.;
				real dsdt_rel_unc_stat = dsdt_unc_stat / dsdt_ref;
				real dsdt_rel_unc_syst = dsdt_unc_syst / dsdt_ref;

				draw((t_cen - t_unc, dsdt_rel)--(t_cen + t_unc, dsdt_rel), p);

				draw((t_repr, dsdt_rel - dsdt_rel_unc_stat)--(t_repr, dsdt_rel + dsdt_rel_unc_stat), p);

				DrawSystematicMark(t_cen - t_unc, t_cen + t_unc, dsdt_rel, dsdt_rel_unc_syst, p);
			}

			//AddToLegend(replace(datasets[dsi], "_", " "), p);
		}
	}

	AddToLegend("bin width", mHL+4pt);
	AddToLegend("statistical unc.", mVL+4pt);
	AddToLegend("systematic unc.", mSY+4pt);

	draw((t_min, 0)--(t_max, 0), heavygreen);

	limits((t_min, -0.20), (t_max, +0.20), Crop);
}

//----------------------------------------------------------------------------------------------------

void DrawOneFit(string type)
{
	string fit_file = topDir + "data/fits/" + version + "/fit_" + type + ".root";

	Results r = LoadResults(RootGetObject(fit_file, "final/results"));

	NewRow();

	NewPad(false);
	string l = "\vbox{\hbox{\bf version: " + version + "}\hbox{\bf type: " + type +  "}\hbox{iteration: " + fit_iteration + "}";
	l += format("\hbox{$\ch^2/ndf = %#.2f", r.chi2) + format("/%.0f", r.ndf) + format("= %#.2f$}", r.chi2_norm);
	l += format("\hbox{$\hbox{p-value} = %#.2f$}", r.prob);
	l += format("\hbox{$\et = (%#.2f", r.eta) + format("\pm %#.2f)$}", r.eta_unc);
	l += format("\hbox{$A = (%#.2f", r.A) + format("\pm %#.2f)\un{mb/GeV^2}$}", r.A_unc);
	l += format("\hbox{$b_1 = (%#.2f", r.b1) + format("\pm %#.2f)\un{GeV^{-2}}$}", r.b1_unc);
	l += format("\hbox{$\rh = (%#.2f", r.rho) + format("\pm %#.2f)$}", r.rho_unc);
	l += format("\hbox{$\si_{\rm tot} = (%#.2f", r.si_tot) + format("\pm %#.2f)\un{mb}$}", r.si_tot_unc);
	l += "}";
	label(l);

	DrawOne(fit_file, 0, 0.005, 6cm);
	DrawOne(fit_file, 0, 0.060, 10cm);
	AttachLegend();

	NewRow();
	NewPad(false);
	DrawOneRelative(fit_file, 0, 0.005, 6cm);
	DrawOneRelative(fit_file, 0, 0.060, 10cm);
	AttachLegend();
}

//----------------------------------------------------------------------------------------------------

for (int ti : types.keys)
{
	DrawOneFit(types[ti]);
}

GShipout(hSkip=1mm, vSkip=0mm);
