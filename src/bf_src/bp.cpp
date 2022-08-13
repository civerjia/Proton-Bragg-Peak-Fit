#define _USE_MATH_DEFINES
#include "bp.h"
#include "parabolic_cylinder_function.h"
#include <omp.h>
#include <cmath>
#include <vector>
double BP::IDD(double depth, double R0, double sigma, double epsilon, double Phi)
{
    // reference : Bortfeld T. An analytical approximation of the Bragg curve for therapeutic proton beams. Med Phys. 1997; 24(12) : 2024 - 2033. doi : 10.1118 / 1.598116
    double dose{};
    // analytical bragg curve
    if (depth < R0 - 10.0 * sigma)
    {
        dose = Bragg_hat(depth, R0, sigma, epsilon, Phi);
    }
    else if (R0 - 10.0 * sigma <= depth && depth <= R0 + 5.0 * sigma)
    {
        dose = Bragg(depth, R0, sigma, epsilon, Phi);
    }
    else
    {
        dose = 0.0;
    }
    return dose;
}
double BP::IDDV2(double depth, double R0, double sigma, double epsilon, double Phi, double a, double b, double c, double d)
{
    // reference : Zhang, Xiaodong, et al. "Parameterization of multiple Bragg curves for scanning proton beams using simultaneous fitting of multiple curves." Physics in Medicine & Biology 56.24 (2011): 7725.
    double dose{};
    // analytical bragg curve
    if (depth < R0 - 10.0 * sigma)
    {
		dose = Bragg_hat(depth, R0, sigma, epsilon, Phi) + a * pow(R0 - depth, 3.0) + b * pow(R0 - depth, 2.0) + c * (R0 - depth) + d;
    }
    else if (R0 - 10.0 * sigma <= depth && depth <= R0 + 5.0 * sigma)
    {
        dose = Bragg(depth, R0, sigma, epsilon, Phi);
    }
    else
    {
        dose = 0.0;
    }
    return dose;
}
double BP::Bragg_hat(double z, double r, double s, double e, double p)
{
    double val{};
    val = (p / (1 + 0.012 * r)) * (17.93 * pow(r - z, -0.435) + (0.444 + 31.7 * e / r) * pow(r - z, 0.565));
    return val;
}
double BP::Bragg(double z, double r, double s, double e, double p)
{
    double val{};
    val = (p / (1 + 0.012 * r)) * std::exp(-(r - z) * (r - z) / (4.0 * s * s)) * pow(s, 0.565) * ((11.26 / s) * PCF::call_D(-(r - z) / s, -0.565) + (0.157 + 11.26 * e / r) * PCF::call_D(-(r - z) / s, -1.565));
    return val;
}
void BP::IDD_array(double* depth, double* dose_o, int size, double R0, double sigma, double epsilon, double Phi)
{
#pragma omp parallel for
    for (int i = 0; i < size; ++i)
    {
        *(dose_o + i) = BP::IDD(depth[i], R0, sigma, epsilon, Phi);
    }
}
void BP::IDD_array_v2(double* depth, double* dose_o, int size, double R0, double sigma, double epsilon, double Phi, double a, double b, double c, double d)
{
#pragma omp parallel for
    for (int i = 0; i < size; ++i)
    {
        *(dose_o + i) = BP::IDDV2(depth[i], R0, sigma, epsilon, Phi, a, b, c, d);
    }
}
void BP::grad(double z, double* grad4_o, double r, double s, double e, double p)
{
    double dBhat_dr{}, dBhat_de{}, dBhat_dp{}, dB_dr{}, dB_ds{}, dB_de{}, dB_dp{};
    dBhat_dr = (-0.012 * p * (17.93 / pow((r - z), 0.435) + (0.444 + 86.1695 / r) * pow((r - z), 0.565))) / pow((1 + 0.012 * r), 2.0) +
        (p * (-7.79955 / pow((r - z), 1.435) + (0.565 * (0.444 + 86.1695 / r)) / pow((r - z), 0.435) - (86.1695 * pow((r - z), 0.565)) / (r * r))) / (1 + 0.012 * r);
    dBhat_de = (31.7 * p * pow((r - z), 0.565)) / (r + 0.012 * r * r);
    dBhat_dp = (17.93 / pow((r - z), 0.435) + (0.444 + 86.1695 / r) * pow((r - z), 0.565)) / (1 + 0.012 * r);

    dB_dr = (-0.012 * p * pow(s, 0.565) * ((11.26 * PCF::call_D(-((r - z) / s), -0.565)) / s + (0.157 + (11.26 * e) / r) * PCF::call_D(-((r - z) / s), -1.565))) / (std::exp(((r - z) * (r - z)) / (4.0 * s * s)) * pow((1 + 0.012 * r), 2.0)) -
        (p * (r - z) * ((11.26 * PCF::call_D(-((r - z) / s), -0.565)) / s + (0.157 + (11.26 * e) / r) * PCF::call_D(-((r - z) / s), -1.565))) / (2 * std::exp(((r - z) * (r - z)) / (4.0 * s * s)) * (1 + 0.012 * r) * pow(s, 1.435)) +
        (p * pow(s, 0.565) * ((-11.26 * e * PCF::call_D(-((r - z) / s), -1.565)) / (r * r) - (11.26 * PCF::call_dD(-((r - z) / s), -0.565)) / (s * s) - ((0.157 + (11.26 * e) / r) * PCF::call_dD(-((r - z) / s), -1.565)) / s)) / (exp(((r - z) * (r - z)) / (4.0 * (s * s))) * (1 + 0.012 * r));
    dB_ds = (0.565 * p * ((11.26 * PCF::call_D(-(r - z) / s, -0.565)) / s + (0.157 + (11.26 * e) / r) * PCF::call_D(-(r - z) / s, -1.565))) / (std::exp(((r - z) * (r - z)) / (4.0 * (s * s))) * (1 + 0.012 * r) * pow(s, 0.435)) +
        (p * ((r - z) * (r - z)) * ((11.26 * PCF::call_D(-(r - z) / s, -0.565)) / s + (0.157 + (11.26 * e) / r) * PCF::call_D(-(r - z) / s, -1.565))) / (2 * std::exp(((r - z) * (r - z)) / (4.0 * (s * s))) * (1 + 0.012 * r) * pow(s, 2.435)) +
        (p * pow(s, 0.565) * ((-11.26 * PCF::call_D(-(r - z) / s, -0.565)) / (s * s) + (11.26 * (r - z) * PCF::call_dD(-(r - z) / s, -0.565)) / (s * s * s) + ((0.157 + (11.26 * e) / r) * (r - z) * PCF::call_dD(-(r - z) / s, -1.565)) / (s * s))) / (exp(((r - z) * (r - z)) / (4.0 * (s * s))) * (1 + 0.012 * r));

    dB_de = (11.26 * p * (pow(s, 0.565)) * PCF::call_D(-(r - z) / s, -1.565)) / (exp(((r - z) * (r - z)) / (4.0 * (s * s))) * (1 + 0.012 * r) * r);
    dB_dp = (pow(s, 0.565) * ((11.26 * PCF::call_D(-((r - z) / s), -0.565)) / s + (0.157 + (11.26 * e) / r) * PCF::call_D(-((r - z) / s), -1.565))) / (exp(((r - z) * (r - z)) / (4.0 * (s * s))) * (1 + 0.012 * r));
    if (z < r - 10.0 * s)
    {
        grad4_o[0] = dBhat_dr;
        grad4_o[1] = 0.0;
        grad4_o[2] = dBhat_de;
        grad4_o[3] = dBhat_dp;
    }
    else if ((r - 10.0 * s <= z) && (z <= r + 5.0 * s))
    {
        grad4_o[0] = dB_dr;
        grad4_o[1] = dB_ds;
        grad4_o[2] = dB_de;
        grad4_o[3] = dB_dp;
    }
    else
    {
        grad4_o[0] = 0.0;
        grad4_o[1] = 0.0;
        grad4_o[2] = 0.0;
        grad4_o[3] = 0.0;
    }
}
void BP::grad_array(double* z, double* grad_o, int size, double* para)
{

    double grad4[4]{};

    for (int i = 0; i < size; ++i)
    {
        BP::grad(z[i], grad4, para[0], para[1], para[2], para[3]);
        for (int j = 0; j < 4; ++j)
        {
            grad_o[j + i * 4] = grad4[j];
        }
    }
}
double norm_l2(double* data1, double* data2, int size)
{
    double val{ 0.0 };
#pragma omp parallel for reduction(+:val)
    for (int i = 0; i < size; ++i)
    {
        val += (data1[i] - data2[i]) * (data1[i] - data2[i]);
    }
    //return std::sqrt(val);
    return (val);
}
inline void update_gradient(int size, double* z, double* dose, double* idd_dose, double* grad4_o, double* para)
{
    double dose_diff{};
    std::vector<double> grad4(4);
    for (std::size_t i = 0; i < size; ++i)
    {
        dose_diff = idd_dose[i] - dose[i];
        // get gradient 
        BP::grad(z[i], grad4.data(), para[0], para[1], para[2], para[3]);
        // calculate the batch mean of gradient
        for (std::size_t j = 0; j < 4; ++j)
        {
            *(grad4_o + j) += grad4[j] * dose_diff / double(size);
        }
    }
}
void BP::fitBragg(double* z, double* dose, int size, double* para, double* lambda, int niter, double* idd_o, double* loss_o, double* para_o, double& Lmin)
{
    // z : depth in cm
    // dose : measured depth dose
    std::vector<double> para0(4), para1(4), para2(4), para_best(4), s(4);
    double t0{ 1.0 }, t1{ 1.0 };
    double L{};
    // predicted dose of analytical model
    std::vector<double> idd_dose(size);
    // get initial predict

    // record the best loss value
    BP::IDD_array(z, idd_dose.data(), size, para[0], para[1], para[2], para[3]);
    Lmin = 0.5 * norm_l2(dose, idd_dose.data(), size);

    for (int j = 0; j < 4; ++j)
    {
        para_best[j] = para[j];
    }
    para0 = para_best;
    para1 = para_best;
    para2 = para_best;
    for (int i = 0; i < niter; ++i)
    {
        // accelerated gradient method
        t1 = (1.0 + std::sqrt(1.0 + 4.0 * t0 * t0)) / 2.0;
        for (int j = 0; j < 4; ++j)
        {
            s[j] = para1[j] + (para1[j] - para0[j]) * (t0 - 1.0) / t1;
        }
        // update gradient
        BP::IDD_array(z, idd_dose.data(), size, s[0], s[1], s[2], s[3]);
        std::vector<double> grad4(4);
        update_gradient(size, z, dose, idd_dose.data(), grad4.data(), s.data());
        // gradient descent
        for (int j = 0; j < 4; ++j)
        {
            double temp{ s[j] - lambda[j] * grad4[j] };
            // positive proximal projection
            if (temp < 0)
            {
                para2[j] = 0;
            }
            else
            {
                para2[j] = temp;
            }
        }
        t0 = t1;
        para0 = para1;
        para1 = para2;
        BP::IDD_array(z, idd_dose.data(), size, para1[0], para1[1], para1[2], para1[3]);
        L = 0.5 * norm_l2(dose, idd_dose.data(), size);
        loss_o[i] = L;
        if (L < Lmin)
        {
            para_best = para1;
            Lmin = L;
        }
    }
    for (int j = 0; j < 4; ++j)
    {
        para_o[j] = para_best[j];
    }
    for (int j = 0; j < size; ++j)
    {
        idd_o[j] = idd_dose[j];
    }
}


void BP::IDD_array_new(double* depth, double* dose_o, double* para_i, int size, int para_size)
{
#pragma omp parallel for
    for (int i = 0; i < size; ++i)
    {
        double temp{};
        for (int k = 0; k < para_size; k = k + 4)
        {
            temp += BP::IDD(depth[i], para_i[k], para_i[k+1], para_i[k+2], para_i[k+3]);
        }
        dose_o[i] = temp;
    }
}
void BP::get_mean_grad(double* depth, double* grad_o, double* para_i, int size, int para_size)
{
#pragma omp parallel for 
    for (int k = 0; k < para_size; k = k + 4)
    {
        double r = para_i[k];
        double s = para_i[k+1];
        double e = para_i[k+2];
        double p = para_i[k+3];
        double grad4_o[4]{};
        for (int i = 0; i < size; ++i)
        {
            double grad4_temp[4]{};
			BP::grad(depth[i], grad4_temp, r, s, e, p);
            grad4_o[0] += grad4_temp[0];
            grad4_o[1] += grad4_temp[1];
            grad4_o[2] += grad4_temp[2];
            grad4_o[3] += grad4_temp[3];
        }
        grad_o[k] = grad4_o[0] / double(size);
        grad_o[k+1] = grad4_o[1] / double(size);
        grad_o[k+2] = grad4_o[2] / double(size);
        grad_o[k+3] = grad4_o[3] / double(size);
    }
}

void BP::get_jacobian(double* depth, double* grad_o, double* para_i, int size, int para_size)
{
#pragma omp parallel for 
    for (int k = 0; k < para_size; k = k + 4)
    {
        double r = para_i[k];
        double s = para_i[k + 1];
        double e = para_i[k + 2];
        double p = para_i[k + 3];
        for (int i = 0; i < size; ++i)
        {
            double grad4_temp[4]{};
            BP::grad(depth[i], grad4_temp, r, s, e, p);
            grad_o[i+ k*size] = grad4_temp[0];
			grad_o[i + (k + 1) * size] = grad4_temp[1];
            grad_o[i + (k + 2) * size] = grad4_temp[2];
            grad_o[i + (k + 3) * size ] = grad4_temp[3];
        }
    }
}
