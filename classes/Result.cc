#include "Result.hh"

#include <cstdio>
#include <fstream>
#include <sstream>

using namespace std;

//----------------------------------------------------------------------------------------------------

double Result::Get(const string &name, double default_value) const
{
	auto it = data.find(name);
	if (it == data.end())
		return default_value;
	else
		return it->second;
}

//----------------------------------------------------------------------------------------------------

void Result::Set(const string &name, double value)
{
	data[name] = value;
}

//----------------------------------------------------------------------------------------------------

void Result::Print() const
{
	for (const auto &entry : data)
		printf("%s-->%.5E\n", entry.first.c_str(), entry.second);
}

//----------------------------------------------------------------------------------------------------

void Result::Load(const std::string &file)
{
	ifstream f_in(file);

	string line;
	while (getline(f_in, line))
	{
		stringstream ss(line);

		string key;
		getline(ss, key, '=');

		string value;
		getline(ss, value, '=');

		data[key] = atof(value.c_str());
	}
}

//----------------------------------------------------------------------------------------------------

void Result::Write(const std::string &file) const
{
	ofstream f_out;
	f_out.open(file);

	for (const auto &entry : data)
		f_out << entry.first << "=" << entry.second << endl;

	f_out.close();
}

//----------------------------------------------------------------------------------------------------

std::unique_ptr<TH1D> Result::ExportROOT() const
{
	unique_ptr<TH1D> h(new TH1D("", "", data.size(), 0., data.size()));

	int bi = 0;
	for (const auto &it : data)
	{
		bi++;

		h->GetXaxis()->SetBinLabel(bi, it.first.c_str());
		h->SetBinContent(bi, it.second);
	}

	return h;
}

//----------------------------------------------------------------------------------------------------

void Result::ImportROOT(const TH1D *h)
{
	for (int bi = 1; bi <= h->GetNbinsX(); ++bi)
		data[h->GetXaxis()->GetBinLabel(bi)] = h->GetBinContent(bi);
}