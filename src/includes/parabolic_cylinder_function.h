#pragma once
#define _USE_MATH_DEFINES
#include <cmath> // gamma function std::tgamma(arg)
#include <array>
#include <complex>
#include <omp.h>
namespace PCF
{
	double pu(double a, double x);   // idx = 0
	double dpu(double a, double x);
	double pulx(double a, double x);
	double dpulx(double a, double x);
	double D(double a, double x);
	double dD_dx(double a, double x);

	double pv(double a, double x);
	double dpv(double a, double x);
	double pvlx(double a, double x);
	double dpvlx(double a, double x);   // idx = 9
	std::complex<double> cgamma(double x, double y, int kf);
	double pw(double x, double y);
	double dpw(double a, double x);
	double pwlx(double a, double x);
	double dpwlx(double a, double x);
	double call_D(double x, double a);  
	double call_dD(double x, double a); 
	double call_pcf_name(double a, double x, int idx);
	void call_pcf_name_array(double a, double* x, int idx, int x_size, double* val);
}