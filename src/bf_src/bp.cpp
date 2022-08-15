#include "bp.h"

template<class T>
T BP::IDD(T depth, T R0, T sigma, T epsilon, T Phi)
{
    // reference : Bortfeld T. An analytical approximation of the Bragg curve for therapeutic proton beams. Med Phys. 1997; 24(12) : 2024 - 2033. doi : 10.1118 / 1.598116
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
T BP::IDDV2(T depth, T R0, T sigma, T epsilon, T Phi, T a, T b, T c, T d)
{
    // reference : Zhang, Xiaodong, et al. "Parameterization of multiple Bragg curves for scanning proton beams using simultaneous fitting of multiple curves." Physics in Medicine & Biology 56.24 (2011): 7725.
    T dose{};
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
        dose = (T)0.0;
    }
    return dose;
}
template<class T>
T BP::Bragg_hat(T z, T r, T s, T e, T p)
{
    T val{};
    T r1012 = T(1 + 0.012 * r);
    T r_z = r - z;
    val = (p / r1012) * ((T)17.93 * pow(r_z, -(T)0.435) + ((T)0.444 + (T)31.7 * e / r) * pow(r_z, (T)0.565));
    return val;
}
template<class T>
T BP::Bragg(T z, T r, T s, T e, T p)
{
    T val{};
    T r1012 = T(1 + 0.012 * r);
    T r_z = r - z;
    T r_z_sqr = r_z * r_z;
    T s_sqr = s * s;
    T exp_val = exp(-r_z_sqr / (T(4.0) * s_sqr));
    T D0565 = static_cast<T>(PCF::call_D(-r_z / s, -0.565));
    T D1565 = static_cast<T>(PCF::call_D(-r_z / s, -1.565));
    val = (p / r1012) * exp_val * pow(s, (T)0.565) * (((T)11.26 / s) * D0565 + ((T)0.157 + (T)11.26 * e / r) * D1565);
    return val;
}
template<class T>
void BP::IDD_array(T* depth, T* dose_o, int size, T R0, T sigma, T epsilon, T Phi)
{
#pragma omp parallel for
    for (int i = 0; i < size; ++i)
    {
        *(dose_o + i) = BP::IDD(depth[i], R0, sigma, epsilon, Phi);
    }
}
template<class T>
void BP::IDD_array_v2(T* depth, T* dose_o, int size, T R0, T sigma, T epsilon, T Phi, T a, T b, T c, T d)
{
#pragma omp parallel for
    for (int i = 0; i < size; ++i)
    {
        *(dose_o + i) = BP::IDDV2(depth[i], R0, sigma, epsilon, Phi, a, b, c, d);
    }
}
template<class T>
void BP::grad(T z, T* grad4_o, T r, T s, T e, T p)
{
    T dBhat_dr{}, dBhat_de{}, dBhat_dp{}, dB_dr{}, dB_ds{}, dB_de{}, dB_dp{};
    // define some commonly used values
    T r_sqr = r * r;
    T r_z = r - z;
    T r_z_sqr = r_z * r_z;
    T s_sqr = s * s;

    T pow_s_0565 = pow(s, T(0.565));
    T pow_r_z_0565 = pow(r_z, T(0.565));
    T pow_r_z_0435 = pow(r_z, T(0.435));
    T r_hat = T(1 + 0.012 * r);
    T pow_1012r_2 = pow(r_hat, 2);

    T exp_val = exp(r_z_sqr / ((T)4.0 * s_sqr));

    T D_0565 = static_cast<T>(PCF::call_D(-(r_z / s), -0.565));
    T D_1565 = static_cast<T>(PCF::call_D(-(r_z / s), -1.565));

    T dD_0565 = static_cast<T>(PCF::call_dD(-r_z / s, -0.565));
    T dD_1565 = static_cast<T>(PCF::call_dD(-r_z / s, -1.565));

    T B2_scale = ((T)0.157 + ((T)11.26 * e) / r);
    T Bhat2_scale = ((T)0.444 + (T)31.7*e / r);
    dBhat_dr = (-(T)0.012 * p * ((T)17.93 / pow_r_z_0435 + Bhat2_scale * pow_r_z_0565)) / pow_1012r_2 +
        (p * (-(T)7.79955 / pow(r_z, (T)1.435) + ((T)0.565 * Bhat2_scale) / pow_r_z_0435 - ((T)31.7*e * pow_r_z_0565) / (r_sqr))) / r_hat;
    dBhat_de = ((T)31.7 * p * pow_r_z_0565) / (r + (T)0.012 * r_sqr);
    dBhat_dp = ((T)17.93 / pow_r_z_0435 + Bhat2_scale * pow_r_z_0565) / r_hat;

    dB_dr = (-(T)0.012 * p * pow_s_0565 * (((T)11.26 * D_0565) / s + B2_scale * D_1565)) / (exp_val * pow_1012r_2) -
        (p * r_z * (((T)11.26 * D_0565) / s + B2_scale * D_1565)) / ((T)2.0 * exp_val * r_hat * pow(s, (T)1.435)) +
        (p * pow_s_0565 * ((-(T)11.26 * e * D_1565) / (r_sqr) - ((T)11.26 * dD_0565) / s_sqr - (B2_scale * dD_1565) / s)) / (exp_val * r_hat);
    dB_ds = ((T)0.565 * p * (((T)11.26 * D_0565) / s + B2_scale * D_1565)) / (exp_val * r_hat * pow(s, (T)0.435)) +
        (p * r_z_sqr * (((T)11.26 * D_0565) / s + B2_scale * D_1565)) / ((T)2.0 * exp_val * r_hat * pow(s, (T)2.435)) +
        (p * pow_s_0565 * ((-(T)11.26 * D_0565) / s_sqr + ((T)11.26 * r_z * dD_0565) / (s_sqr * s) + (B2_scale * r_z * dD_1565) / s_sqr)) / (exp_val * r_hat);

    dB_de = ((T)11.26 * p * (pow_s_0565) * D_1565) / (exp_val * r_hat * r);
    dB_dp = (pow_s_0565 * (((T)11.26 * D_0565) / s + B2_scale * D_1565)) / (exp_val * r_hat);
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
void BP::grad_array(T* z, T* grad_o, int size, T* para)
{

    T grad4[4]{};

    for (int i = 0; i < size; ++i)
    {
        BP::grad(z[i], grad4, para[0], para[1], para[2], para[3]);
        for (int j = 0; j < 4; ++j)
        {
            grad_o[j + i * 4] = grad4[j];
        }
    }
}
template<class T>
T norm_l2(T* data1, T* data2, int size)
{
    T val{ 0.0 };
#pragma omp parallel for reduction(+:val)
    for (int i = 0; i < size; ++i)
    {
        val += (data1[i] - data2[i]) * (data1[i] - data2[i]);
    }
    //return std::sqrt(val);
    return (val);
}
template<class T>
inline void update_gradient(int size, T* z, T* dose, T* idd_dose, T* grad4_o, T* para)
{
    T dose_diff{};
    std::vector<T> grad4(4);
    for (std::size_t i = 0; i < size; ++i)
    {
        dose_diff = idd_dose[i] - dose[i];
        // get gradient 
        BP::grad(z[i], grad4.data(), para[0], para[1], para[2], para[3]);
        // calculate the batch mean of gradient
        for (std::size_t j = 0; j < 4; ++j)
        {
            *(grad4_o + j) += grad4[j] * dose_diff / T(size);
        }
    }
}
template<class T>
void BP::fitBragg(T* z, T* dose, int size, T* para, T* lambda, int niter, T* idd_o, T* loss_o, T* para_o, T& Lmin)
{
    // z : depth in cm
    // dose : measured depth dose
    std::vector<T> para0(4), para1(4), para2(4), para_best(4), s(4);
    T t0{ 1.0 }, t1{ 1.0 };
    T L{};
    // predicted dose of analytical model
    std::vector<T> idd_dose(size);
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
        std::vector<T> grad4(4);
        update_gradient(size, z, dose, idd_dose.data(), grad4.data(), s.data());
        // gradient descent
        for (int j = 0; j < 4; ++j)
        {
            T temp{ s[j] - lambda[j] * grad4[j] };
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


template<class T>
void BP::IDD_array_N(T* depth, T* dose_o, T* para_i, int size, int para_size)
{
#pragma omp parallel for
    for (int i = 0; i < size; ++i)
    {
        T temp{};
        for (int k = 0; k < para_size; k = k + 4)// multi-bf mixture
        {
            temp += BP::IDD(depth[i], para_i[k], para_i[k + 1], para_i[k + 2], para_i[k + 3]);
        }
        dose_o[i] = temp;
    }
}
template<class T>
void BP::get_mean_grad(T* depth, T* grad_o, T* para_i, int size, int para_size)
{
#pragma omp parallel for 
    for (int k = 0; k < para_size; k = k + 4)
    {
        T r = para_i[k];
        T s = para_i[k + 1];
        T e = para_i[k + 2];
        T p = para_i[k + 3];
        T grad4_o[4]{};
        for (int i = 0; i < size; ++i)
        {
            T grad4_temp[4]{};
            BP::grad(depth[i], grad4_temp, r, s, e, p);
            grad4_o[0] += grad4_temp[0];
            grad4_o[1] += grad4_temp[1];
            grad4_o[2] += grad4_temp[2];
            grad4_o[3] += grad4_temp[3];
        }
        grad_o[k] = grad4_o[0] / T(size);
        grad_o[k + 1] = grad4_o[1] / T(size);
        grad_o[k + 2] = grad4_o[2] / T(size);
        grad_o[k + 3] = grad4_o[3] / T(size);
    }
}

template<class T>
void BP::get_jacobian(T* depth, T* grad_o, T* para_i, int size, int para_size)
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
            BP::grad(depth[i], grad4_temp, r, s, e, p);
            grad_o[i + k * size] = grad4_temp[0];
            grad_o[i + (k + 1) * size] = grad4_temp[1];
            grad_o[i + (k + 2) * size] = grad4_temp[2];
            grad_o[i + (k + 3) * size] = grad4_temp[3];
        }
    }
}


template double BP::IDD(double depth, double R0, double sigma, double epsilon, double Phi);
template double BP::IDDV2(double depth, double R0, double sigma, double epsilon, double Phi, double a, double b, double c, double d);
template double BP::Bragg_hat(double z, double r, double s, double e, double p);
template double BP::Bragg(double z, double r, double s, double e, double p);
template void BP::IDD_array(double* depth, double* dose_o, int size, double R0, double sigma, double epsilon, double Phi);
template void BP::IDD_array_v2(double* depth, double* dose_o, int size, double R0, double sigma, double epsilon, double Phi, double a, double b, double c, double d);
template void BP::grad(double z, double* grad4_o, double r, double s, double e, double p);
template void BP::grad_array(double* z, double* grad_o, int size, double* para);
template void BP::fitBragg(double* z, double* dose, int size, double* para, double* lambda, int niter, double* idd_o, double* loss_o, double* para_o, double& Lmin);
template void BP::IDD_array_N(double* depth, double* dose_o, double* para_i, int size, int para_size);
template void BP::get_mean_grad(double* depth, double* grad_o, double* para_i, int size, int para_size);
template void BP::get_jacobian(double* depth, double* grad_o, double* para_i, int size, int para_size);

template float BP::IDD(float depth, float R0, float sigma, float epsilon, float Phi);
template float BP::IDDV2(float depth, float R0, float sigma, float epsilon, float Phi, float a, float b, float c, float d);
template float BP::Bragg_hat(float z, float r, float s, float e, float p);
template float BP::Bragg(float z, float r, float s, float e, float p);
template void BP::IDD_array(float* depth, float* dose_o, int size, float R0, float sigma, float epsilon, float Phi);
template void BP::IDD_array_v2(float* depth, float* dose_o, int size, float R0, float sigma, float epsilon, float Phi, float a, float b, float c, float d);
template void BP::grad(float z, float* grad4_o, float r, float s, float e, float p);
template void BP::grad_array(float* z, float* grad_o, int size, float* para);
template void BP::fitBragg(float* z, float* dose, int size, float* para, float* lambda, int niter, float* idd_o, float* loss_o, float* para_o, float& Lmin);
template void BP::IDD_array_N(float* depth, float* dose_o, float* para_i, int size, int para_size);
template void BP::get_mean_grad(float* depth, float* grad_o, float* para_i, int size, int para_size);
template void BP::get_jacobian(float* depth, float* grad_o, float* para_i, int size, int para_size);