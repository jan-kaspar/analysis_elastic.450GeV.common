import root;
import pad_layout;
include "../common.asy";
include "shared.asy";

version = "iv2_v3";

string models[];
models.push("exp1_con_rho0.05");
models.push("exp1_con_rho0.10");
models.push("exp1_con_rho0.15");

string fit_type = "HighBeta-sequence-0.005-0.020";

string errors[], e_labels[];
pen e_pens[];
errors.push("none"); e_labels.push("none"); e_pens.push(black);
//errors.push("stat"); e_labels.push("stat"); e_pens.push(blue);
//errors.push("syst"); e_labels.push("syst"); e_pens.push(red);
//errors.push("norm"); e_labels.push("norm"); e_pens.push(magenta);
//errors.push("stat,syst,norm"); e_labels.push("stat,syst,norm"); e_pens.push(heavygreen);

string parameters[], p_captions[], p_labels[];
parameters.push("eta"); p_captions.push("$\eta$"); p_labels.push("$\De\eta$");
parameters.push("A"); p_captions.push("$A$"); p_labels.push("$\De A\ung{mb/GeV^2}$");
parameters.push("b1"); p_captions.push("$b_1$"); p_labels.push("$\De b_1 \ung{GeV^{-2}}$");
//parameters.push("p0"); p_captions.push("$p_0$"); p_labels.push("$\De p_0$");
parameters.push("rho"); p_captions.push("$\rho$"); p_labels.push("$\De\rho$");
parameters.push("si_tot"); p_captions.push("$\si_{\rm tot}$"); p_labels.push("$\De\si_{\rm tot}$");

//----------------------------------------------------------------------------------------------------

NewPad(false);
AddToLegend("input version: " + replace(version_input, "_", "\_"));
AddToLegend("binning: " + replace(binning, "_", "\_"));
AddToLegend("fit type: " + replace(fit_type, "_", "\_"));
AttachLegend();

for (int pari : parameters.keys)
	NewPadLabel(p_captions[pari]);

for (int mi : models.keys)
{
	string model = models[mi];

	NewRow();

	NewPadLabel(replace(models[mi], "_", "\_"));

	for (int pari : parameters.keys)
	{
		NewPad(p_labels[pari]);

		for (int ei : errors.keys)
		{
			string f = topDir + "data/simu/" + version + "/" + model + "/errors_" + errors[ei] + "/process_fits_" + fit_type + ".root";
			string o = parameters[pari] + "/h_diff";
			draw(RootGetObject(f, o), "vl,lR", e_pens[ei], errors[ei]);
		}

		AttachLegend(BuildLegend(NW), NE);
	}
}

GShipout(vSkip=0mm);
