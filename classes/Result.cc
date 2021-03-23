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
