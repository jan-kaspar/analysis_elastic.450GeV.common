import root;
import pad_layout;
include "shared.asy";

string topDir = "../../";

string datasets[], d_labels[];
datasets.push("low_beta"); d_labels.push("low beta");
datasets.push("high_beta"); d_labels.push("high beta");

string errors[], e_labels[];
errors.push("stat"); e_labels.push("stat");
errors.push("syst"); e_labels.push("syst");
errors.push("norm"); e_labels.push("norm");
errors.push("full"); e_labels.push("full");

//----------------------------------------------------------------------------------------------------

void DrawRelUnc(RootObject h, pen p)
{
	guide g;

	for (int bi = 1; bi <= h.iExec("GetNbinsX"); ++bi)
	{
		real l = h.rExec("GetBinLowEdge", bi);
		real r = l + h.rExec("GetBinWidth", bi);

		real v = h.rExec("GetBinContent", bi);
		real u = h.rExec("GetBinError", bi);

		real u_rel = (v > 0) ? u/v : 0.;

		g = g--(l, u_rel)--(r, u_rel);
	}

	draw(g, p);
}

//----------------------------------------------------------------------------------------------------

NewPad(false);
for (int ei : errors.keys)
	NewPadLabel(e_labels[ei]);

for (int dsi : datasets.keys)
{
	NewRow();

	NewPadLabel(d_labels[dsi]);

	for (int ei : errors.keys)
	{
		NewPad("$|t|\ung{GeV^2}$", "$\d\si/\d t$ relative error");

		string f = topDir + "data/simu/" + version + "/" + model + "/errors_" + errors[ei] + "/process_histograms.root";
		string o = binning + "/" + datasets[dsi] + "/h_stddev_rel";
		draw(RootGetObject(f, o), "vl,eb", red);

		pen p_ref = blue;

		if (errors[ei] == "stat")
		{
			string f = topDir + "data/input/" + version_input + "/import.root";
			string o = binning + "/" + datasets[dsi] + "/h_dsdt_cen_stat";
			DrawRelUnc(RootGetObject(f, o), blue);
		}

		if (errors[ei] == "syst")
		{
			string f = topDir + "data/input/" + version_input + "/import.root";
			string o = binning + "/" + datasets[dsi] + "/h_dsdt_rel_syst";
			draw(RootGetObject(f, o), "vl", blue);
		}

		if (errors[ei] == "norm")
		{
			real unc = 0.10;

			draw((0, unc)--(0.1, unc), blue);
		}
	}
}
