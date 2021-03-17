#ifndef _Stat_hh_
#define _Stat_hh_

#include "TMatrixDSym.h"

#include <vector>
#include <string>
#include <cmath>

using namespace std;

//----------------------------------------------------------------------------------------------------

struct Stat
{
	unsigned int dim;
	double S1;
	vector<double> Sv, Svv, Svvv, Svvvv;

	vector< vector<double> > Sxy, Sxxy, Sxyy, Sxxyy;

	vector<string> labels;

	void Init(unsigned int _dim = 1);

	Stat() {}

	Stat(unsigned int _dim)
	{
		Init(_dim);
	}

	void SetLabels(const vector<string> &_l);

	template <class T>
	void Fill(const T &v)
	{
		S1 += 1.;
		for (unsigned int i = 0; i < dim; i++)
		{
			Sv[i] += v[i];
			Svv[i] += v[i]*v[i];
			Svvv[i] += v[i]*v[i]*v[i];
			Svvvv[i] += v[i]*v[i]*v[i]*v[i];

			for (unsigned int j = 0; j < dim; j++)
			{
				Sxy[i][j] += v[i] * v[j];
				Sxxy[i][j] += v[i]*v[i] * v[j];
				Sxyy[i][j] += v[i] * v[j]*v[j];
				Sxxyy[i][j] += v[i]*v[i] * v[j]*v[j];
			}
		}
	}

	void Fill(double v1, double v2 = 0., double v3 = 0., double v4 = 0., double v5 = 0.);

	string QLabel(unsigned int i) const;
	
	//--------------------
	// 1D getters
	//--------------------

	double GetEntries() const;

	double GetMean(unsigned int i) const;

	double GetStdDev(unsigned int i) const;

	double GetMeanUnc(unsigned int i) const;

	double GetStdDevUnc(unsigned int i) const;

	// approximation of GetStdDevUnc valid for Gaussian distributions
	double GetStdDevUncGauss(unsigned int i) const;
	
	//--------------------
	// 2D getters
	//--------------------

	double GetCovariance(unsigned int i, unsigned int j) const;

	double GetCorrelation(unsigned int i, unsigned int j) const;

	double GetCovarianceUnc(unsigned int i, unsigned int j) const;

	double GetCorrelationUnc(unsigned int i, unsigned int j) const;

	TMatrixDSym GetCovarianceMatrix() const;

	//--------------------
	// print methods
	//--------------------

	void PrintStat() const;

	void PrintMeanAndStdDev() const;
	
	void PrintCovariance() const;

	void PrintCorrelation() const;
};

#endif
