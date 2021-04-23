import root;
import pad_layout;
include "../common.asy";

string elegentDir = "/afs/cern.ch/work/j/jkaspar/work/software/Elegent/production/CMSSW_11_0_0/";

//----------------------------------------------------------------------------------------------------

xSizeDef = 10cm;
ySizeDef = 8cm;

//----------------------------------------------------------------------------------------------------

string model_tags[];
pen model_pens[];

void AddModel(string tag, pen p)
{
	model_tags.push(tag);
	model_pens.push(gray);
}

AddModel("block [06]", heavygreen);
AddModel("bourrely [03]", blue);
AddModel("dl [13]", heavygreen + dashed);
AddModel("ferreira [14]", blue + dashed);
AddModel("godizov [14]", red + dashed);
AddModel("islam (hp) [06,09]", black);
AddModel("islam (lxg) [06,09]", black+dashed);
AddModel("jenkovszky [11]", magenta);
//AddModel("petrov (2p) [02]", red+dashed);
AddModel("petrov (3p) [02]", red);

//----------------------------------------------------------------------------------------------------

string modes[] = { "PH" };

//----------------------------------------------------------------------------------------------------

string f_elegent = elegentDir + "scripts/distributions/t-distributions,pp,900GeV.root";

string f_fit = topDir + "data/simu/neco/models.root";

TGraph_x_max = 10.0;

//xTicksDef = LeftTicks(0.005, 0.001);

for (int mdi : modes.keys)
{
	NewPad("$|t|\ung{GeV^2}$", "$\d\si/\d t\ung{mb/GeV^2}$");
	scale(Linear, Log);

	for (int mi : model_tags.keys)
	{
		TGraph_reducePoints = 10;
		draw(RootGetObject(f_elegent, "full range/" + model_tags[mi] + "/" + modes[mdi] + "/differential cross-section", search=false), "l", model_pens[mi], model_tags[mi]);	
	}

	TGraph_reducePoints = 1;
	draw(RootGetObject(f_fit, "exp1_con_rho0.10/g_dsdt_H"), "l", red + 2pt, "fit model");

	//limits((0, 4e1), (0.03, 1e4), Crop);

	if (mdi == 0)
		AttachLegend(BuildLegend(NW), NE);
}
