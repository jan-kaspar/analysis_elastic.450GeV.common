#include "Stat.hh"

//----------------------------------------------------------------------------------------------------

void Stat::Init(unsigned int _dim)
{
	dim = _dim;

	S1 = 0.;
	for (unsigned int i = 0; i < dim; i++)
	{
		Sv.push_back(0);
		Svv.push_back(0);
		Svvv.push_back(0);
		Svvvv.push_back(0);

		vector<double> temp;
		for (unsigned int j = 0; j < dim; j++)
		{
			temp.push_back(0);
		}
		Sxy.push_back(temp);
		Sxxy.push_back(temp);
		Sxyy.push_back(temp);
		Sxxyy.push_back(temp);
	}
}

//----------------------------------------------------------------------------------------------------

void Stat::SetLabels(const vector<string> &_l)
{
	labels.resize(dim);
	for (unsigned int i = 0; i < dim; i++)
		labels[i] = _l[i];
}

//----------------------------------------------------------------------------------------------------

void Stat::Fill(double v1, double v2, double v3, double v4, double v5)
{
	vector<double> v(5);
	v[0] = v1;
	v[1] = v2;
	v[2] = v3;
	v[3] = v4;
	v[4] = v5;

	Fill(v);
}

//----------------------------------------------------------------------------------------------------

string Stat::QLabel(unsigned int i) const
{
	if (labels.empty())
	{
		char buf[10];
		sprintf(buf, "qu.%3i", i);
		return buf;
	} else
		return labels[i];
}

//----------------------------------------------------------------------------------------------------

double Stat::GetEntries() const
{
	return S1;
}

//----------------------------------------------------------------------------------------------------

double Stat::GetMean(unsigned int i) const
{
	double mu = Sv[i] / S1;
	return mu;
}

//----------------------------------------------------------------------------------------------------

double Stat::GetStdDev(unsigned int i) const
{
	double v = (Svv[i] - Sv[i]*Sv[i] / S1) / (S1 - 1.);
	double s = (v > 0.) ? sqrt(v) : 0.;
	return s;
}

//----------------------------------------------------------------------------------------------------

double Stat::GetMeanUnc(unsigned int i) const
{
	double mu_unc = GetStdDev(i) / sqrt(S1);
	return mu_unc;
}

//----------------------------------------------------------------------------------------------------

double Stat::GetStdDevUnc(unsigned int i) const
{
	double mu = GetMean(i);
	double s = GetStdDev(i);
	double v = s*s;

	double sum = Svvvv[i] - 4.*mu*Svvv[i] + 6.*mu*mu*Svv[i] - 4.*mu*mu*mu*Sv[i] + mu*mu*mu*mu*S1;
	double E4 = (S1 > 1.) ? sum / (S1 - 1.) : 0.;

	double v_var = (S1 > 3.) ? (E4 - (S1 - 3.)/(S1 - 1.)*v*v) / S1 : 0.;
	double s_var = v_var / 4. / v;
	double s_s = (s_var > 0.) ? sqrt(s_var) : 0.;
	return s_s;
}

//----------------------------------------------------------------------------------------------------

double Stat::GetStdDevUncGauss(unsigned int i) const
{
	double s = GetStdDev(i);
	double s_s = (S1 > 0.) ? s / sqrt(2. * S1) : 0.;
	return s_s;
}

//----------------------------------------------------------------------------------------------------

double Stat::GetCovariance(unsigned int i, unsigned int j) const
{
	double C = (S1 > 1.) ? (Sxy[i][j] - Sv[i]*Sv[j] / S1) / (S1 - 1.) : 0.;
	return C;
}

//----------------------------------------------------------------------------------------------------

double Stat::GetCorrelation(unsigned int i, unsigned int j) const
{
	double C = GetCovariance(i, j);
	double den = GetStdDev(i) * GetStdDev(j);
	double rho = (den > 0.) ? C / den : 0.;
	return rho;
}

//----------------------------------------------------------------------------------------------------

double Stat::GetCovarianceUnc(unsigned int i, unsigned int j) const
{
	double mx = GetMean(i);
	double my = GetMean(j);
	double sx = GetStdDev(i);
	double sy = GetStdDev(j);
	double C = GetCovariance(i, j);

	double sum =
		Sxxyy[i][j] 
		-2.*Sxyy[i][j]*mx - 2.*Sxxy[i][j]*my
		+ 4.*Sxy[i][j]*mx*my
		+ Svv[i]*my*my + Svv[j]*mx*mx
		- 2.*Sv[i]*mx*my*my - 2.*Sv[j]*mx*mx*my
		+ mx*mx*my*my;
	double D = (S1 > 1.) ? sum / (S1 - 1.) : 0.;

	double C_var = (S1 > 2.) ? (D + sx*sx*sy*sy/(S1 - 1.) - (S1-2.)/(S1-1.)*C*C) / S1 : 0.;
	double C_s = (C_var > 0.) ? sqrt(C_var) : 0.;

	return C_s;
}

//----------------------------------------------------------------------------------------------------

double Stat::GetCorrelationUnc(unsigned int i, unsigned int j) const
{
	// WARNING: the calculation below assumes no correlation between C, si_i and si_j, which
	// might not be correct - in that case it gives an upper bound for the uncertainty

	double C = GetCovariance(i, j), C_unc = GetCovarianceUnc(i, j);
	double si_i = GetStdDev(i), si_i_unc = GetStdDevUnc(i);
	double si_j = GetStdDev(j), si_j_unc = GetStdDevUnc(j);
	double rho = C / (si_i * si_j);
	double sum =
		(C != 0. && si_i != 0. && si_j != 0.) ? pow(C_unc / C, 2.) + pow(si_i_unc / si_i, 2.) + pow(si_j_unc / si_j, 2.) : 0.;
	double rho_unc = fabs(rho) * sqrt(sum);
	return rho_unc;
}

//----------------------------------------------------------------------------------------------------

TMatrixDSym Stat::GetCovarianceMatrix() const
{
	TMatrixDSym m(dim);

	for (unsigned int i = 0; i < dim; i++)
		for (unsigned int j = 0; j < dim; j++)
			m(i, j) = GetCovariance(i, j);

	return m;
}

//----------------------------------------------------------------------------------------------------

void Stat::PrintStat() const
{
	printf("entries: %.3E\n", S1);
}

//----------------------------------------------------------------------------------------------------

void Stat::PrintMeanAndStdDev() const
{
	for (unsigned int i = 0; i < dim; i++)
	{
		double mu = GetMean(i);
		double mu_unc = GetMeanUnc(i);
		double s = GetStdDev(i);
		double s_unc = GetStdDevUnc(i);
		printf("%s: mean %+.3E +- %.3E, std. dev. = %.3E +- %.3E\n", QLabel(i).c_str(), mu, mu_unc, s, s_unc);
	}
}

//----------------------------------------------------------------------------------------------------

void Stat::PrintCovariance() const
{
	printf("      ");
	for (unsigned int i = 0; i < dim; i++)
		printf("   %6s", QLabel(i).c_str());
	printf("\n");

	for (unsigned int i = 0; i < dim; i++)
	{
		printf("%6s", QLabel(i).c_str());
		for (unsigned int j = 0; j < dim; j++)
			printf("   %+.3f", GetCovariance(i, j));
		printf("\n");
	}
}

//----------------------------------------------------------------------------------------------------

void Stat::PrintCorrelation() const
{
	printf("      ");
	for (unsigned int i = 0; i < dim; i++)
		printf("   %6s", QLabel(i).c_str());
	printf("\n");

	for (unsigned int i = 0; i < dim; i++)
	{
		printf("%6s", QLabel(i).c_str());
		for (unsigned int j = 0; j < dim; j++)
			printf("   %+.3f", GetCorrelation(i, j));
		printf("\n");
	}
}
