import root;
import pad_layout;
include "../common.asy";
//include "shared.asy";

//string version = "iv4-2_v1";
string version = "TEST";

string types[], t_steps[][];
types.push("HighBeta-single-nuis-0.020"); t_steps.push(new string[] {""});
types.push("HighBeta-single-nuis-0.020-fix-b1"); t_steps.push(new string[] {""});

string fit_iteration = "iteration 2";

string datasets[];
pen d_pens[];
datasets.push("high_beta"); d_pens.push(red + squarecap);
//datasets.push("low_beta"); d_pens.push(blue + squarecap);

TGraph_errorBar = None;

xSizeDef = 12cm;
xTicksDef = LeftTicks(0.005, 0.001);

//----------------------------------------------------------------------------------------------------

void DrawOneFit(string type)
{
	string fit_file = topDir + "data/fits/" + version + "/fit_" + type + "/fit.root";

	NewPad("$|t|\ung{GeV^2}$", "systematic effect of $\d\si/\d t\ung{mb/GeV^2}$");
	
	for (int dsi : datasets.keys)
	{
		draw(RootGetObject(fit_file, fit_iteration + "/" + datasets[dsi] + "/g_syst_eff"), "l,p", red, mCi+1pt+d_pens[dsi]);
	}

	//limits((0, -3), (0.02, +3), Crop);
}

//----------------------------------------------------------------------------------------------------

for (int ti : types.keys)
{
	for (int si : t_steps[ti].keys)
		DrawOneFit(types[ti] + t_steps[ti][si]);

	GShipout("syst_effect_" + types[ti] + ".pdf", hSkip=1mm, vSkip=0mm);
}

