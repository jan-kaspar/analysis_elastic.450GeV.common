import root;
import pad_layout;
include "../common.asy";
include "shared.asy";

version = "iv2_TEST";

string f = topDir + "data/simu/" + version + "/models.root";

string f_fit = "/afs/cern.ch/work/j/jkaspar/work/analyses/elastic/450GeV/common/master/benchmarks/test2.root";

string models[] = {
	"exp1_con_rho0.05",
	"exp1_con_rho0.10",
	"exp1_con_rho0.15"
};

string model_ref = "exp1_con_rho0.10";

TGraph_x_min = 3e-4;
TGraph_x_max = 0.05;

//----------------------------------------------------------------------------------------------------

void DrawDiff(RootObject g, RootObject g_ref, bool rel, pen p, string l="")
{
	guide gd;

	for (int i = 0; i < g.iExec("GetN"); ++i)
	{
		real ax[] = {0};
		real ay[] = {0};
		g.vExec("GetPoint", i, ax, ay);

		if (ax[0] < TGraph_x_min || ax[0] > TGraph_x_max)
			continue;

		real y_ref = g_ref.rExec("Eval", ax[0]);
		
		real v = ay[0] - y_ref;
		if (rel)
			v /= y_ref;

		gd = gd -- (ax[0], v);
	}

	draw(gd, p, l);
}

//----------------------------------------------------------------------------------------------------

RootObject g_ref = RootGetObject(f, model_ref + "/g_dsdt_CH");

NewPad("$|t|\ung{GeV^2}$", "$\d\si^{\rm model}/\d t - \d\si^{\rm reference}/\d t \ung{mb/GeV^2}$");

for (int mi : models.keys)
	DrawDiff(RootGetObject(f, models[mi] + "/g_dsdt_CH"), g_ref, false, StdPen(mi), replace(models[mi], "_", "\_"));

RootObject fit = RootGetObject(f_fit, "fit_0.005_0.020|ff");
TF1_x_min = -inf; TF1_x_max = +inf; draw(fit, green);
TF1_x_min = 0; TF1_x_max = 0.005; draw(fit, green+dotted);

RootObject fit = RootGetObject(f_fit, "fit_0.015_0.025|ff");
TF1_x_min = -inf; TF1_x_max = +inf; draw(fit, magenta);
TF1_x_min = 0; TF1_x_max = 0.015; draw(fit, magenta+dotted);

limits((0, -10), (0.025, +10), Crop);


NewPad("$|t|\ung{GeV^2}$", "${\d\si^{\rm model}/\d t - \d\si^{\rm reference}\d t\over \d\si^{\rm reference}/\d t} \ung{mb/GeV^2}$");

for (int mi : models.keys)
	DrawDiff(RootGetObject(f, models[mi] + "/g_dsdt_CH"), g_ref, true, StdPen(mi), replace(models[mi], "_", "\_"));
limits((0, -0.06), (0.025, +0.06), Crop);



//----------------------------------------------------------------------------------------------------

frame f_legend = BuildLegend();
NewPad(false);
add(f_legend);


GShipout(hSkip=5mm, vSkip=0mm);
