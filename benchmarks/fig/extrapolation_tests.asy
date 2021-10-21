import root;
import pad_layout;

string topDir = "../";

string f = topDir + "test1.root";

string fits[];
pen f_pens[];
fits.push("from_0.005_to_0.015"); f_pens.push(blue);
fits.push("from_0.005_to_0.020"); f_pens.push(red);

TH1_x_max = 0.026;

//----------------------------------------------------------------------------------------------------

real GetParameter(RootObject hist, string parName)
{
	for (int bi = 1; bi <= hist.iExec("GetNbinsX"); ++bi)
	{
		string binName = hist.oExec("GetXaxis").sExec("GetBinLabel", bi);
		if (binName == parName)
			return hist.rExec("GetBinContent", bi);
	}

	return -1234;
}

//----------------------------------------------------------------------------------------------------

NewPad("$|t|\ung{GeV}$", "$\d\si/\d t\ung{mb/GeV^2}$");

draw(RootGetObject(f, "h_dsdt_high"), "eb", black);

for (int fi : fits.keys)
{
	RootObject fit = RootGetObject(f, "extrapolation checks/" + fits[fi] + "/h_dsdt_high_with_fit|ff");
	RootObject results = RootGetObject(f, "extrapolation checks/" + fits[fi] + "/h_results");

	real B = fit.rExec("GetParameter", 1);
	real B_unc = fit.rExec("GetParError", 1);

	string l = replace(fits[fi], "_", " ");

	TF1_x_min = -inf;
	TF1_x_max = +inf;
	draw(fit, f_pens[fi], l);

	TF1_x_min = 0;
	TF1_x_max = 0.025;
	draw(fit, f_pens[fi] + dashed);

	AddToLegend(format("$A = (%.1f", GetParameter(results, "A")) + format(" \pm %.1f)\un{mv/GeV^2}$", GetParameter(results, "A_unc")));
	AddToLegend(format("$B = (%.2f", GetParameter(results, "B")) + format(" \pm %.2f)\un{GeV^{-2}}$", GetParameter(results, "B_unc")));
	AddToLegend(format("$\rho = (%.2f", GetParameter(results, "rho")) + format(" \pm %.2f)$", GetParameter(results, "rho_unc")));
	AddToLegend(format("$\si_{\rm tot} = (%.2f", GetParameter(results, "si_tot")) + format(" \pm %.2f)\un{mb}$", GetParameter(results, "si_tot_unc")));
}

limits((0, 100), (0.025, 300), Crop);
AttachLegend(NW, NE);
