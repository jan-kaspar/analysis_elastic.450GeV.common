import root;
import pad_layout;

string f = "../coulomb_effect_on_B.root";

string rhos[], rho_labels[];
rhos.push("0.10"); rho_labels.push("$\rh = 0.10$");
rhos.push("0.15"); rho_labels.push("$\rh = 0.15$");

string cnis[], cni_labels[];
cnis.push("KL"); cni_labels.push("KL");
cnis.push("AmpSum"); cni_labels.push("amplitude sum");

string ranges[];
pen r_pens[];
ranges.push("0.005_0.020"); r_pens.push(red);
ranges.push("0.005_0.010"); r_pens.push(blue);
ranges.push("0.010_0.015"); r_pens.push(magenta);
ranges.push("0.015_0.020"); r_pens.push(heavygreen);
ranges.push("0.020_0.025"); r_pens.push(cyan);

xSizeDef = 8cm;

//----------------------------------------------------------------------------------------------------

NewPad(false);
for (int ii : cnis.keys)
	NewPadLabel(cnis[ii]);

for (int rhi : rhos.keys)
{
	NewRow();

	NewPadLabel(rho_labels[rhi]);

	for (int ii : cnis.keys)
	{
		NewPad("$|t|\ung{GeV^2}$", "$\d\si^{C+H}/\d t\ung{mb/GeV^2}$");
		//scale(Linear, Log);

		string dir = "exp1_con_rho" + rhos[rhi] + "_" + cnis[ii];
		draw(RootGetObject(f, dir + "/g_dsdt_CH"), black);

		for (int ri : ranges.keys)
		{
			RootObject fit = RootGetObject(f, dir + "/fit_" + ranges[ri]);
			string lab = replace(ranges[ri], "_", "--") + format(": $B = %#.2f$", fit.rExec("GetParameter", 1));
			draw(fit, "l", r_pens[ri], lab);
		}

		limits((0, 150), (0.03, 250), Crop);
		AttachLegend();
	}
}

GShipout(vSkip=0mm);
