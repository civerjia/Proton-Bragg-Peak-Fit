#pragma once
namespace BP
{
	double IDD(double depth, double R0, double sigma, double epsilon, double Phi);
	double IDDV2(double depth, double R0, double sigma, double epsilon, double Phi, double a, double b, double c, double d);
	double Bragg_hat(double z, double r, double s, double e, double p);
	double Bragg(double z, double r, double s, double e, double p);
	void IDD_array(double* depth, double* dose_o, int size, double R0, double sigma, double epsilon, double Phi);
	void IDD_array_v2(double* depth, double* dose_o, int size, double R0, double sigma, double epsilon, double Phi, double a, double b, double c, double d);
	void grad(double z, double* grad4_o, double r, double s, double e, double p);
	void grad_array(double* z, double* grad_o, int size, double* para);
	void fitBragg(double* z, double* dose, int size, double* para, double* lambda, int niter, double* idd_o, double* loss_o, double* para_o, double& Lmin);
	void IDD_array_new(double* depth, double* dose_o, double* para_i, int size, int para_size);
	void get_mean_grad(double* depth, double* grad_o, double* para_i, int size, int para_size);
	void get_jacobian(double* depth, double* grad_o, double* para_i, int size, int para_size);
}