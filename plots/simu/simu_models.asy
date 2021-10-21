import root;
import pad_layout;
include "../common.asy";
include "shared.asy";

version = "iv2_v3";

string models[] = {
	"exp1_con_rho0.05",
	"exp1_con_rho0.10",
	"exp1_con_rho0.15"
};

TGraph_x_min = 3e-4;
TGraph_x_max = 0.1;

//----------------------------------------------------------------------------------------------------

void DrawOne(string object, string y_label, bool y_log, real y_min, real y_max)
{
	NewPad("$|t|\ung{GeV^2}$", y_label);
	
	if (y_log)
		scale(Linear, Log);
	
	for (int mi : models.keys)
	{
		string f = topDir + "data/simu/" + version + "/models.root";
		pen p = StdPen(mi);
		draw(RootGetObject(f, models[mi] + "/" + object), p, replace(models[mi], "_", "\_"));
	}

	limits((0, y_min), (0.03, y_max), Crop);

	yaxis(XEquals(0.005, false), green);
	yaxis(XEquals(0.020, false), magenta);
}

//----------------------------------------------------------------------------------------------------

DrawOne("g_dsdt_H", "$\d\si^{\rm H}/\d t\ung{mb/GeV^2}$", false, 120, 250);

DrawOne("g_dsdt_CH", "$\d\si^{\rm C+H}/\d t\ung{mb/GeV^2}$", false, 120, 250);

DrawOne("g_dsdt_CH", "$\d\si^{\rm C+H}/\d t\ung{mb/GeV^2}$", true, 100, 1000);

frame f_legend = BuildLegend();
NewPad(false);
add(f_legend);


GShipout(hSkip=1mm, vSkip=0mm);
