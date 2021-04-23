import root;
import pad_layout;
include "../common.asy";
//include "shared.asy";

string input_file = topDir + "data/input/version_0/import.root";

string binning = "sb1";

string datasets[], d_labels[];
pen d_pens[];
real d_scale[];
datasets.push("high_beta"); d_labels.push("high $\be^*$"); d_pens.push(red + squarecap); d_scale.push(1);
datasets.push("low_beta"); d_labels.push("low $\be^*$"); d_pens.push(blue + squarecap); d_scale.push(1);

TGraph_errorBar = None;

//----------------------------------------------------------------------------------------------------

void DrawSystematics(transform t, RootObject h_cen, RootObject h_rel, pen p)
{
	for (int bi = 1; bi <= h_cen.iExec("GetNbinsX"); ++bi)
	{
		real t_left = h_cen.rExec("GetBinLowEdge", bi);
		real t_right = t_left + h_cen.rExec("GetBinWidth", bi);

		real dsdt = h_cen.rExec("GetBinContent", bi);

		real dsdt_unc_rel = h_rel.rExec("GetBinContent", bi);
		real dsdt_unc = dsdt * dsdt_unc_rel;

		if (dsdt <= 0)
			continue;

		draw(t * (Scale((t_left, dsdt-dsdt_unc))--Scale((t_right, dsdt-dsdt_unc))), p);
		draw(t * (Scale((t_left, dsdt+dsdt_unc))--Scale((t_right, dsdt+dsdt_unc))), p);
	}
}

//----------------------------------------------------------------------------------------------------

void DrawOne(int h_idx, real t_min, real t_max, real x_size)
{
	NewPad("$|t|\ung{GeV^2}$", "$\d\si/\d t\ung{mb/GeV}$", xSize = x_size);
	SetGridHint(h_idx, 0);

	scale(Linear, Log);

	for (int dsi : datasets.keys)
	{
		string ds = datasets[dsi]; 
		pen p = d_pens[dsi];

		RootObject h_dsdt_cen_stat = RootGetObject(input_file, binning + "/" + ds + "/h_dsdt_cen_stat");
		RootObject h_dsdt_rel_syst = RootGetObject(input_file, binning + "/" + ds + "/h_dsdt_rel_syst");

		real s = d_scale[dsi];
		transform t = shift(0, log10(s));

		draw(t, h_dsdt_cen_stat, "vl,eb", p, d_labels[dsi]);
		DrawSystematics(t, h_dsdt_cen_stat, h_dsdt_rel_syst, p);
	}

	xlimits(t_min, t_max, Crop);
}

//----------------------------------------------------------------------------------------------------

DrawOne(0, 0, 0.005, 6cm);

DrawOne(1, 0, 0.100, 10cm);

AttachLegend();

GShipout(hSkip=1mm, vSkip=0mm);
