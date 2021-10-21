import root;
import pad_layout;

string topDir = "../";

string f = topDir + "test1.root";

string fits[];
pen f_pens[];
fits.push("from_5.0E-03_to_1.0E-02"); f_pens.push(blue);
fits.push("from_1.0E-02_to_1.5E-02"); f_pens.push(red);
fits.push("from_1.5E-02_to_2.0E-02"); f_pens.push(cyan);
fits.push("from_2.0E-02_to_2.5E-02"); f_pens.push(magenta);

TH1_x_max = 0.026;

//----------------------------------------------------------------------------------------------------

NewPad("$|t|\ung{GeV}$", "$\d\si/\d t\ung{mb/GeV^2}$");

draw(RootGetObject(f, "h_dsdt_high"), "eb", black);

for (int fi : fits.keys)
{
	RootObject fit = RootGetObject(f, "slope checks/" + fits[fi] + "|ff");

	real B = fit.rExec("GetParameter", 1);
	real B_unc = fit.rExec("GetParError", 1);

	string l = replace(fits[fi], "_", " ") + format(": $B = (%.2f", B) + format(" \pm %.2f)\un{GeV^{-2}}$", B_unc);

	draw(fit, f_pens[fi], l);
}

limits((0, 100), (0.025, 300), Crop);
AttachLegend(NW, NE);
