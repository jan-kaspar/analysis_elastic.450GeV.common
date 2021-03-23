#ifndef _Result_hh_
#define _Result_hh_

#include "TH1D.h"

#include <map>
#include <memory>
#include <string>

//----------------------------------------------------------------------------------------------------

struct Result
{
	std::map<std::string, double> data;

	Result() {}
	Result(const TH1D *h) { ImportROOT(h); }

	double Get(const std::string &name, double default_value = -1) const;
	void Set(const std::string &name, double value);

	void Load(const std::string &file);
	void Write(const std::string &file) const;

	std::unique_ptr<TH1D> ExportROOT() const;
	void ImportROOT(const TH1D *h);
};

#endif
