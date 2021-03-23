#ifndef _Result_hh_
#define _Result_hh_

#include <map>
#include <string>
#include <cmath>

//----------------------------------------------------------------------------------------------------

struct Result
{
	std::map<std::string, double> data;

	double Get(const std::string &name, double default_value = NAN) const;
	void Set(const std::string &name, double value);

	void Load(const std::string &file);
	void Write(const std::string &file) const;
};

#endif
