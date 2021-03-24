import root;
import pad_layout;
include "../common.asy";
//include "shared.asy";

string fit_file = topDir + "data/simu/bla/exp1_con_rho0.10/errors_none/seed_1/fit_both.root";

string fit_iteration = "iteration 2";

string datasets[];
pen d_pens[];
datasets.push("high_beta"); d_pens.push(red + squarecap);
datasets.push("low_beta"); d_pens.push(blue + squarecap);

TGraph_errorBar = None;

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

void DrawOne(int h_idx, real t_min, real t_max, real x_size)
{
	NewPad("$|t|\ung{GeV^2}$", "$\d\si/\d t\ung{mb/GeV}$", xSize = x_size);
	SetGridHint(h_idx, 0);

	scale(Linear, Log);

	for (int dsi : datasets.keys)
	{
		string ds = datasets[dsi]; 
		pen p = d_pens[dsi];

		RootObject g_data_t_unc = RootGetObject(fit_file, fit_iteration + "/" + ds + "/g_data_t_unc");
		RootObject g_data_dsdt_unc_stat = RootGetObject(fit_file, fit_iteration + "/" + ds + "/g_data_dsdt_unc_stat");
		RootObject g_data_dsdt_unc_syst = RootGetObject(fit_file, fit_iteration + "/" + ds + "/g_data_dsdt_unc_syst");

		draw(g_data_t_unc, "p", p);
		draw(g_data_dsdt_unc_stat, "p", p);
		DrawSystematics(g_data_t_unc, g_data_dsdt_unc_syst, p);
	}

	RootObject g_fit_CH = RootGetObject(fit_file, fit_iteration + "/g_fit_CH");

	draw(g_fit_CH, "l", heavygreen);

	xlimits(t_min, t_max, Crop);
}

//----------------------------------------------------------------------------------------------------

DrawOne(0, 0, 0.005, 6cm);

DrawOne(1, 0, 0.100, 10cm);

GShipout(vSkip=0mm);
