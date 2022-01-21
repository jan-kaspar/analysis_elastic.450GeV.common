import root;
import pad_layout;
include "../common.asy";
//include "shared.asy";

//string version = "iv4-2_v1";
string version = "TEST";

string types[], t_steps[][];
types.push("HighBeta-single-nuis-0.020"); t_steps.push(new string[] {""});
types.push("HighBeta-single-nuis-0.020-fix-b1"); t_steps.push(new string[] {""});

//----------------------------------------------------------------------------------------------------

struct NuisancePar
{
	string tag;
	real val, unc;
};

NuisancePar nps[];

//----------------------------------------------------------------------------------------------------

void Update(string tag, bool set_val, real val, bool set_unc, real unc)
{
	bool found;
	for (NuisancePar np : nps)
	{
		if (np.tag == tag)
		{
			if (set_val) np.val = val;
			if (set_unc) np.unc = unc;

			found = true;
			break;
		}
	}

	if (!found)
	{
		NuisancePar np;
		np.tag = tag;
		np.val = val;
		np.unc = unc;

		nps.push(np);
	}
}

//----------------------------------------------------------------------------------------------------

string NPTickLabels(real y)
{
	if (y >=0 && y < nps.length)
	{
		return nps[(int) y].tag;
	} else {
		return "";
	}
}

yTicksDef = RightTicks(Label(""), NPTickLabels, Step=1, step=0);

//----------------------------------------------------------------------------------------------------

void DrawOneFit(string type)
{
	// get results
	string fit_file = topDir + "data/fits/" + version + "/fit_" + type + "/fit.root";
	RootObject h_data = RootGetObject(fit_file, "final/results");

	// extract list of nuisance parameters
	nps.delete();

	for (int bi = 1; bi <= h_data.iExec("GetNbinsX"); ++bi)
	{
		RootObject axis = h_data.oExec("GetXaxis");
		string l = axis.sExec("GetBinLabel", bi);

		if (find(l, "nuis:") != 0)
			continue;

		real bc = h_data.rExec("GetBinContent", bi);

		l = substr(l, 5);

		if (find(l, "_unc") > 0)
		{
			l = substr(l, 0, length(l) - 4);
			Update(l, false, 0, true, bc);
		} else {
			Update(l, true, bc, false, 0);
		}
	}

	// make plot
	NewPad("nuisance parameter", "", ySize = 12cm);
	for (int npi : nps.keys)
	{
		//write(np.tag + ": ", np.val, np.unc);

		NuisancePar np = nps[npi];
		real y = npi;

		draw((np.val, y), mCi+2pt+blue);
		draw((np.val - np.unc, y)--(np.val + np.unc, y), blue);
	}

	limits((-2, -1), (+2, nps.length), Crop);
}

//----------------------------------------------------------------------------------------------------

for (int ti : types.keys)
{
	for (int si : t_steps[ti].keys)
		DrawOneFit(types[ti] + t_steps[ti][si]);

	GShipout("nuisance_parameters_" + types[ti] + ".pdf", hSkip=1mm, vSkip=0mm);
}

