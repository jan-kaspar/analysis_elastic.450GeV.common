import root;
import pad_layout;
include "../common.asy";
include "shared.asy";

string errors[], e_labels[];
pen e_pens[];
//errors.push("none"); e_labels.push("none"); e_pens.push(black);
//errors.push("stat"); e_labels.push("stat"); e_pens.push(blue);
//errors.push("syst"); e_labels.push("syst"); e_pens.push(red);
//errors.push("norm"); e_labels.push("norm"); e_pens.push(magenta);
errors.push("stat,syst,norm"); e_labels.push("stat,syst,norm"); e_pens.push(heavygreen);

string fit_types[], f_labels[];
fit_types.push("HighBeta-single-0.020"); f_labels.push("HighBeta-single-0.020");

//----------------------------------------------------------------------------------------------------

for (int ei : errors.keys)
{
	for (int fi : fit_types.keys)
	{
		NewRow();

		NewPadLabel("\vbox{\hbox{" + e_labels[ei] + "}\hbox{" + f_labels[fi] + "}}");
		
		string f = topDir + "data/simu/" + version + "/" + model + "/errors_" + errors[ei] + "/process_fits_" + fit_types[fi] + ".root";

		NewPad("$\De B\ung{GeV^2}$", "$\De\rh$");
		RootObject h2 = RootGetObject(f, "correlations/rho vs B");
		draw(h2);
		real r = 8, sl = 8.28e-3;
		draw((-r, -r * sl)--(+r, +r * sl), dashed);
		limits((-8, -0.1), (+8, +0.1), Crop);

		NewPad("$\De\si_{\rm tot}\ung{mb}$", "$\De\rh$");
		RootObject h2 = RootGetObject(f, "correlations/rho vs si_tot");
		draw(h2);
		real r = 3, sl = -8.48e-3;
		draw((-r, -r * sl)--(+r, +r * sl), dashed);
		limits((-3, -0.1), (+3, +0.1), Crop);
	}
}
