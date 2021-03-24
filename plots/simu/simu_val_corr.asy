import root;
import pad_layout;
include "../common.asy";
include "shared.asy";

string datasets[], d_labels[];
datasets.push("low_beta"); d_labels.push("low beta");
datasets.push("high_beta"); d_labels.push("high beta");

string errors[], e_labels[];
errors.push("stat"); e_labels.push("stat");
errors.push("syst"); e_labels.push("syst");
errors.push("norm"); e_labels.push("norm");
errors.push("stat,syst,norm"); e_labels.push("stat,syst,norm");

TH2_palette = Gradient(blue, yellow, red);

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
		NewPad("bin index", "bin index");

		string f = topDir + "data/simu/" + version + "/" + model + "/errors_" + errors[ei] + "/process_histograms.root";
		string o = binning + "/" + datasets[dsi] + "/h2_rel_corr";
		draw(RootGetObject(f, o));
	}
}
