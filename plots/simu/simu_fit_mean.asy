import root;
import pad_layout;
include "../common.asy";
include "shared.asy";

string fit_types[], f_labels[];
fit_types.push("high"); f_labels.push("high beta");
fit_types.push("both"); f_labels.push("high and low beta");

string errors[], e_labels[];
pen e_pens[];
errors.push("none"); e_labels.push("none"); e_pens.push(black);
errors.push("stat"); e_labels.push("stat"); e_pens.push(blue);
errors.push("syst"); e_labels.push("syst"); e_pens.push(red);
errors.push("norm"); e_labels.push("norm"); e_pens.push(magenta);
errors.push("stat,syst,norm"); e_labels.push("stat,syst,norm"); e_pens.push(heavygreen);

string parameters[], p_captions[], p_labels[];
parameters.push("eta"); p_captions.push("$\eta$"); p_labels.push("$\De\eta$");
parameters.push("A"); p_captions.push("$A$"); p_labels.push("$\De A\ung{mb/GeV^2}$");
parameters.push("b1"); p_captions.push("$b_1$"); p_labels.push("$\De b_1 \ung{GeV^{-2}}$");
//parameters.push("p0"); p_captions.push("$p_0$"); p_labels.push("$\De p_0$");
parameters.push("rho"); p_captions.push("$\rho$"); p_labels.push("$\De\rho$");
parameters.push("si_tot"); p_captions.push("$\si_{\rm tot}$"); p_labels.push("$\De\si_{\rm tot}$");

//----------------------------------------------------------------------------------------------------

string FillTickLabels(real x)
{
	if (x >=0 && x < errors.length)
	{
		return e_labels[(int) x];
	} else {
		return "";
	}
}

yTicksDef = RightTicks(rotate(0)*Label(""), FillTickLabels, Step=1, step=0);

//----------------------------------------------------------------------------------------------------

NewPad(false);
for (int pari : parameters.keys)
	NewPadLabel(p_captions[pari]);

for (int fi : fit_types.keys)
{
	NewRow();

	NewPadLabel(f_labels[fi]);

	for (int pari : parameters.keys)
	{
		NewPad(p_labels[pari]);

		for (int ei : errors.keys)
		{
			string f = topDir + "data/simu/" + version + "/" + model + "/errors_" + errors[ei] + "/process_fits_" + fit_types[fi] + ".root";
			RootObject g_data = RootGetObject(f, parameters[pari] + "/g_data");

			real ax[] = {0.};
			real ay[] = {0.};

			g_data.vExec("GetPoint", 0, ax, ay); real mean = ax[0], mean_unc = ay[0];
			g_data.vExec("GetPoint", 1, ax, ay); real rms = ax[0], rms_unc = ay[0];

			real y = ei;
			real dy = 0.1;

			draw((mean, y-dy)--(mean, y+dy), red);
			draw((mean - mean_unc, y)--(mean + mean_unc, y), red);
		}

		draw((0, 0)--(0, errors.length-1), heavygreen + dashed);

		//AttachLegend(BuildLegend(NW), NE);
	}
}

GShipout(vSkip=0mm);
