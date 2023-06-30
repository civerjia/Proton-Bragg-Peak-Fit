#include "bf.h"
#include <iostream>

// see: Bortfeld, T. (1999). An analytical approximation of the Bragg curve for therapeutic proton beams. Med. Phys. 26, 1847-1851.
// bortfeld funtion contains parabolic cylinder function, which is slow to compute
// this code is based on the approximation of the parabolic cylinder function
// author: Shuang Zhou 
// date: 2023-06-29
// email: 740110712@qq.com
template <class T>
T Dn0565(T z)
{
    // this is the approximation of the parabolic cylinder function when z->-inf with fixed a=-1.565
    T scale = 1.591188981970586;                   // sqrt(2pi)/gamma(0.565)
    T term1 = pow(-z, T(-0.435)) * exp(z * z / T(4.0)); // (-z)^(-a-1)e^(-z^2/4)
    // 1+(a+1)(a+2)/(2*z^2)
    T term2 = T(1.0) + T(0.624225) / (T(2.0) * z * z);
    return scale * term1 * term2;
}
template <class T>
T Dn1565(T z)
{
    // this is the approximation of the parabolic cylinder function when z->inf with fixed a=-1.565
    T scale = 2.816263684903692;                  // sqrt(2pi)/gamma(1.565)
    T term1 = pow(-z, T(0.565)) * exp(z * z / T(4.0)); // (-z)^(-a-1)e^(-z^2/4)
    // 1+(a+1)(a+2)/(2*z^2)
    T term2 = T(1.0) - T(0.245775) / (T(2.0) * z * z);
    return scale * term1 * term2;
}
template <class T>
T Dp0565(T z)
{
    // this is the approximation of the parabolic cylinder function when z->inf with fixed a=-1.565
    T term1 = pow(z, T(-0.565)) * exp(-z * z / T(4.0)); // (z)^(a)e^(-z^2/4)
    // 1-a(a-1)/(4*z^2)
    T term2 = T(1.0) - T(0.884225) / (T(4.0) * z * z);
    return term1 * term2;
}
template <class T>
T Dp1565(T z)
{
    // this is the approximation of the parabolic cylinder function when z->inf with fixed a=-1.565
    T term1 = pow(z, T(-1.565)) * exp(-z * z / T(4.0)); // (z)^(a)e^(-z^2/4)
    // 1-a(a-1)/(4*z^2)
    T term2 = T(1.0) - T(4.014225) / (T(4.0) * z * z);
    return term1 * term2;
}

template <class T>
T F1(T a, T b, T x)
{
    // this is a approximation of the hypergeometric function 1F1(a,b,x)
    T b_x = b - x;
    T y = T(2.0) * a / (b_x + sqrt(b_x * b_x + T(4.0) * a * x));
    T term1 = pow(b, b - T(0.5)) / sqrt(y * y / a + (T(1.0) - y) * (T(1.0) - y) / (b - a));
    T term2 = pow(y / a, a) * pow((T(1.0) - y) / (b - a), b - a);
    T term3 = exp(x * y);
    return term1 * term2 * term3;
}
template <class T>
T Dmid0565(T z)
{
    // this is the approximation of the parabolic cylinder function when |z|<1.6 with fixed a=-0.565
    // 2^(a/2)e^(-z^2/4)
    T scale = T(0.822165078627196) * exp(-z * z / T(4.0));
    // sqrt(pi)/gamma((1-a)/2)
    T term1 = T(1.496398901525717) * F1(T(0.2825), T(0.5), z * z / T(2.0));
    // sqrt(2pi)/gamma(-a/2)
    T term2 = T(0.786549997591402) * z * F1(T(0.7825), T(1.5), z * z / T(2.0));
    return scale * (term1 - term2);
}
template <class T>
T Dmid1565(T z)
{
    // this is the approximation of the parabolic cylinder function when |z|<1.6 with fixed a=-1.565
    // 2^(a/2)e^(-z^2/4)
    T scale = T(0.581358502352061) * exp(-z * z / T(4.0));
    // sqrt(pi)/gamma((1-a)/2)
    T term1 = T(1.968760485094311) * F1(T(0.7825), T(0.5), z * z / T(2.0));
    // sqrt(2pi)/gamma(-a/2)
    T term2 = T(2.116227621257871) * z * F1(T(1.2825), T(1.5), z * z / T(2.0));
    return scale * (term1 - term2);
}
template <class T>
T D0565(T z)
{
    if (z < -1.9)
    {
        return Dn0565(z);
    }
    else if (z > 1.8)
    {
        return Dp0565(z);
    }
    else
    {
        return Dmid0565(z);
    }
}
template <class T>
T D1565(T z)
{
    if (z < -1.4)
    {
        return Dn1565(z);
    }
    else if (z > 1.6)
    {
        return Dp1565(z);
    }
    else
    {
        return Dmid1565(z);
    }
}

template <class T>
T Bragg_hat(T z, T r, T s, T e, T p)
{
    // z : depth
    // r : range R80 (80% of the maximum dose)
    // s : sigma (energy spread)
    // e : epsilon
    // p : phi (fluence)
    T r1012 = T(1 + 0.012 * r);
    T r_z = r - z;
    return (p / r1012) * ((T)17.93 * pow(r_z, -(T)0.435) + ((T)0.444 + (T)31.7 * e / r) * pow(r_z, (T)0.565));
}

template<class T>
T Bragg(T z, T r, T s, T e, T p)
{
    // z : depth
    // r : range R80 (80% of the maximum dose)
    // s : sigma (energy spread)
    // e : epsilon Fraction of primary fluence contribut-ing to the ‘‘tail’’ of the energy spectrum
    // p : phi (fluence)
    T r1012 = T(1 + 0.012 * r);
    T r_z = r - z;
    T r_z_sqr = r_z * r_z;
    T s_sqr = s * s;
    T exp_val = exp(-r_z_sqr / (T(4.0) * s_sqr));
    return (p / r1012) * exp_val * pow(s, (T)0.565) * (((T)11.26 / s) * D0565(-r_z / s) + ((T)0.157 + (T)11.26 * e / r) * D1565(-r_z / s));
}

template<class T>
T IDD(T depth, T R0, T sigma, T epsilon, T Phi)
{
    // reference : Bortfeld T. An analytical approximation of the Bragg curve for therapeutic proton beams. Med Phys. 1997; 24(12) : 2024 - 2033. doi : 10.1118 / 1.598116
    // epsilon : affect the slope of the entrance dose
    T dose{};
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
        dose = (T)0.0;
    }
    return dose;
}

template<class T>
void grad(T z, T* grad4_o, T r, T s, T e, T p)
{
    // analytical gradient of Bragg_hat
    T dBhat_dr{}, dBhat_de{}, dBhat_dp{};
    // define some commonly used values
    T r_sqr = r * r;
    T r_z = r - z;

    T pow_r_z_0565 = pow(r_z, T(0.565));
    T pow_r_z_0435 = pow(r_z, T(0.435));
    T r_hat = T(1 + 0.012 * r);
    T pow_1012r_2 = pow(r_hat, 2);

    T Bhat2_scale = ((T)0.444 + (T)31.7*e / r);
    dBhat_dr = (-(T)0.012 * p * ((T)17.93 / pow_r_z_0435 + Bhat2_scale * pow_r_z_0565)) / pow_1012r_2 +
        (p * (-(T)7.79955 / pow(r_z, (T)1.435) + ((T)0.565 * Bhat2_scale) / pow_r_z_0435 - ((T)31.7*e * pow_r_z_0565) / (r_sqr))) / r_hat;
    dBhat_de = ((T)31.7 * p * pow_r_z_0565) / (r + (T)0.012 * r_sqr);
    dBhat_dp = ((T)17.93 / pow_r_z_0435 + Bhat2_scale * pow_r_z_0565) / r_hat;

    // use numerical differentiation for the gradient of Bragg function
    T dB_dr{}, dB_ds{}, dB_de{}, dB_dp{};
    T h = 1e-6;
    T B0, B1;
    B0 = IDD(z, r, s, e, p);
    B1 = IDD(z, r + h, s, e, p);
    dB_dr = (B1 - B0) / h;
    B1 = IDD(z, r, s + h, e, p);
    dB_ds = (B1 - B0) / h;
    B1 = IDD(z, r, s, e + h, p);
    dB_de = (B1 - B0) / h;
    B1 = IDD(z, r, s, e, p + h);
    dB_dp = (B1 - B0) / h;

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

template<class T>
void BP_fast::IDD_array_N(T* depth, T* dose_o, T* para_i, int size, int para_size)
{
#pragma omp parallel for
    for (int i = 0; i < size; ++i)
    {
        T temp{};
        for (int k = 0; k < para_size; k = k + 4)// multi-bf mixture
        {
            temp += IDD(depth[i], para_i[k], para_i[k + 1], para_i[k + 2], para_i[k + 3]);
        }
        dose_o[i] = temp;
    }
}

template<class T>
void BP_fast::get_jacobian(T* depth, T* grad_o, T* para_i, int size, int para_size)
{
#pragma omp parallel for 
    for (int k = 0; k < para_size; k = k + 4)
    {
        T r = para_i[k];
        T s = para_i[k + 1];
        T e = para_i[k + 2];
        T p = para_i[k + 3];
        for (int i = 0; i < size; ++i)
        {
            T grad4_temp[4]{};
            grad(depth[i], grad4_temp, r, s, e, p);
            grad_o[i + k * size] = grad4_temp[0];
            grad_o[i + (k + 1) * size] = grad4_temp[1];
            grad_o[i + (k + 2) * size] = grad4_temp[2];
            grad_o[i + (k + 3) * size] = grad4_temp[3];
        }
    }
}
template void BP_fast::IDD_array_N(double* depth, double* dose_o, double* para_i, int size, int para_size);
template void BP_fast::get_jacobian(double* depth, double* grad_o, double* para_i, int size, int para_size);
template void BP_fast::IDD_array_N(float* depth, float* dose_o, float* para_i, int size, int para_size);
template void BP_fast::get_jacobian(float* depth, float* grad_o, float* para_i, int size, int para_size);


int main(){
    std::cout << "Hello World!\n";
    std::cout << Dmid0565(-1.0) << std::endl;
}