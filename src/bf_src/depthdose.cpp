#include <cmath>
#include <omp.h>
#include <vector>
#include "mex.h"
#include "matrix.h"

float p_pos565_f[5] = {0.019680,-5.017049,-147.880863,82.969327,-181.139311};
float q_pos565_f[4] = {92.365631,230.948158,-58.187022,343.580628};
float p_neg565_f[5] = {0.016028,-0.890829,-5.293047,-16.249304,-19.868614};
float q_neg565_f[4] = {1.678568,6.734086,6.837704,38.044967};

float p_pos1565_f[5] = {-0.037653,-5.735521,-31.968762,-130.893484,7.228251};
float q_pos1565_f[4] = {11.740156,52.921439,120.487416,-6.724969};
float p_neg1565_f[5] = {-0.017402,2.611226,1.037757,18.393094,0.268162};
float q_neg1565_f[4] = {-2.449875,4.034374,-16.971481,-0.246841};

double p_pos565[5] = {0.019680,-5.017049,-147.880863,82.969327,-181.139311};
double q_pos565[4] = {92.365631,230.948158,-58.187022,343.580628};
double p_neg565[5] = {0.016028,-0.890829,-5.293047,-16.249304,-19.868614};
double q_neg565[4] = {1.678568,6.734086,6.837704,38.044967};

double p_pos1565[5] = {-0.037653,-5.735521,-31.968762,-130.893484,7.228251};
double q_pos1565[4] = {11.740156,52.921439,120.487416,-6.724969};
double p_neg1565[5] = {-0.017402,2.611226,1.037757,18.393094,0.268162};
double q_neg1565[4] = {-2.449875,4.034374,-16.971481,-0.246841};

// function f(x) = (p0*x^5 + p1*x^4 + p2*x^3 + p3*x^2 + p4*x)/(x^4+q0*x^3+q1*x^2+q2*x+q3)
template<class T>
T f(T *p, T *q, T x){
    T numerator = x * (p[4] + x * (p[3] + x * (p[2] + x * (p[1] + x * p[0]))));
    T denominator = q[3] + x * (q[2] + x * (q[1] + x * (q[0] + x)));
    return numerator / denominator;
}

float D0565(float x){
    // parabolic cylinder function D(-0.565, x)
    if (x < 0){
        return 1.230286920530541f*std::exp(f(p_neg565_f, q_neg565_f, x) + x*x/4.0f);
    }
    else{
        return 1.230286920530541f*std::exp(f(p_pos565_f, q_pos565_f, x) - x*x/4.0f);
    }
}
double D0565(double x){
    // parabolic cylinder function D(-0.565, x)
    if (x < 0){
        return 1.230286920530541*std::exp(f(p_neg565, q_neg565, x) + x*x/4.0);
    }
    else{
        return 1.230286920530541*std::exp(f(p_pos565, q_pos565, x) - x*x/4.0);
    }
}
float f0565(float x){
    // parabolic cylinder function D(-0.565, x)*exp(-x^2/4)
    if (x < 0){
        return 1.230286920530541f*std::exp(f(p_neg565_f, q_neg565_f, x));
    }
    else{
        return 1.230286920530541f*std::exp(f(p_pos565_f, q_pos565_f, x) - x*x/2.0f);
    }
}
double f0565(double x){
    // parabolic cylinder function D(-0.565, x)*exp(-x^2/4)
    if (x < 0){
        return 1.230286920530541*std::exp(f(p_neg565, q_neg565, x));
    }
    else{
        return 1.230286920530541*std::exp(f(p_pos565, q_pos565, x) - x*x/2.0);
    }
}

float D1565(float x){
    // parabolic cylinder function D(-1.565, x)
    if (x < 0){
        return 1.144555647104346f*std::exp(f(p_neg1565_f, q_neg1565_f, x) + x*x/4.0f);
    }
    else{
        return 1.144555647104346f*std::exp(f(p_pos1565_f, q_pos1565_f, x) - x*x/4.0f);
    }
}
double D1565(double x){
    // parabolic cylinder function D(-1.565, x)
    if (x < 0){
        return 1.144555647104346*std::exp(f(p_neg1565, q_neg1565, x) + x*x/4.0);
    }
    else{
        return 1.144555647104346*std::exp(f(p_pos1565, q_pos1565, x) - x*x/4.0);
    }
}

float f1565(float x){
    // parabolic cylinder function D(-1.565, x)*exp(-x^2/4)
    if (x < 0){
        return 1.144555647104346f*std::exp(f(p_neg1565_f, q_neg1565_f, x));
    }
    else{
        return 1.144555647104346f*std::exp(f(p_pos1565_f, q_pos1565_f, x) - x*x/2.0f);
    }
}
double f1565(double x){
    // parabolic cylinder function D(-1.565, x)*exp(-x^2/4)
    if (x < 0){
        return 1.144555647104346*std::exp(f(p_neg1565, q_neg1565, x));
    }
    else{
        return 1.144555647104346*std::exp(f(p_pos1565, q_pos1565, x) - x*x/2.0);
    }
}

float Bragg(float z, float r, float s, float e, float p){
    // z : depth
    // r : range R80 (80% of the maximum dose)
    // s : sigma (energy spread)
    // e : epsilon Fraction of primary fluence contribut-ing to the ‘‘tail’’ of the energy spectrum
    // p : phi (fluence)
    float r1012 = (1.0f + 0.012f * r);
    float x = (r - z) / s;
    // return (p / r1012) * std::pow(s, 0.565f) * ((11.26f / s) * f0565(-x) + (0.157f + 11.26f * e / r) * f1565(-x));
    return (p / r1012) * std::pow(s, 0.565f) * ((11.26 / s) * f0565(-x) + (0.157f + 11.26f * e / r) * f1565(-x));
}
double Bragg(double z, double r, double s, double e, double p){
    double r1012 = (1.0 + 0.012 * r);
    double x = (r - z) / s;
    // return (p / r1012) * std::pow(s, 0.565) * ((11.26 / s) * f0565(-x) + (0.157 + 11.26 * e / r) * f1565(-x));
    return (p / r1012) * std::pow(s, 0.565) * ((11.26 / s) * f0565(-x) + (0.157 + 11.26 * e / r) * f1565(-x));
}

float Bragg_hat(float z, float r, float s, float e, float p)
{
    // z : depth
    // r : range R80 (80% of the maximum dose)
    // s : sigma (energy spread)
    // e : epsilon
    // p : phi (fluence)
    float r1012 = (1.0f + 0.012f * r);
    float r_z = r - z;
    // return (p / r1012) * (17.93f * std::pow(r_z, -0.435f) + (0.444f + 31.7f * e / r) * std::pow(r_z, 0.565f));
    return (p / r1012) * (17.974454939899370f * std::pow(r_z, -0.435f) + (0.441595774885904f+31.671136466339373f * e / r) * std::pow(r_z, 0.565f));
}
double Bragg_hat(double z, double r, double s, double e, double p)
{
    double r1012 = (1.0 + 0.012 * r);
    double r_z = r - z;
    // return (p / r1012) * (17.93 * std::pow(r_z, -0.435) + (0.444 + 31.7 * e / r) * std::pow(r_z, 0.565));
    return (p / r1012) * (17.974454939899370 * std::pow(r_z, -0.435) + (0.441595774885904+31.671136466339373 * e / r) * std::pow(r_z, 0.565));
}

template<class T>
T IDD(T depth, T R0, T sigma, T epsilon, T Phi)
{
    // reference : Bortfeld float. An analytical approximation of the Bragg curve for therapeutic proton beams. Med Phys. 1997; 24(12) : 2024 - 2033. doi : 10.1118 / 1.598116
    T dose{};
    T edge = R0 - T(10.0) * sigma;
    // analytical bragg curve
    if (depth < edge)
    {
        dose = Bragg_hat(depth, R0, sigma, epsilon, Phi);
    }
    else if (edge <= depth && depth <= R0 + T(5.0) * sigma)
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
T IDDV2(T depth, T R0, T sigma, T epsilon, T Phi, T a, T b, T c)
{
    // reference : Zhang, Xiaodong, et al. "Parameterization of multiple Bragg curves for scanning proton beams using simultaneous fitting of multiple curves." Physics in Medicine & Biology 56.24 (2011): 7725.
    // fixed discontinuity of Zhang's formula
    T dose{};
    T edge = R0 - T(10.0) * sigma;
    // analytical bragg curve
    if (depth < edge)
    {
        dose = Bragg_hat(depth, R0, sigma, epsilon, Phi) + a * std::pow(edge - depth, 3.0) + b * std::pow(edge - depth, 2.0) + c * (edge - depth);
    }
    else if (edge <= depth && depth <= R0 + T(5.0) * sigma)
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
void IDD_array_v2(T* depth, T* dose_o, T* para_i, int N_depth, int para_size)
{
#pragma omp parallel for
    for (int i = 0; i < N_depth; ++i)
    {
        T temp{};
        for (int k = 0; k < para_size; k = k + 7)// multi-bf mixture
        {
            temp += IDDV2(depth[i], para_i[k], para_i[k + 1], para_i[k + 2], para_i[k + 3], para_i[k + 4], para_i[k + 5], para_i[k + 6]);
        }
        dose_o[i] = temp;
    }
}
template<class T>
void IDD_array_v1(T* depth, T* dose_o, T* para_i, int N_depth, int para_size)
{
#pragma omp parallel for
    for (int i = 0; i < N_depth; ++i)
    {
        T temp{};
        for (int k = 0; k < para_size; k = k + 4)// multi-bf mixture
        {
            temp += IDD(depth[i], para_i[k], para_i[k + 1], para_i[k + 2], para_i[k + 3]);
        }
        dose_o[i] = temp;
    }
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mxClassID id0 = mxGetClassID(prhs[0]);
    mxClassID id1 = mxGetClassID(prhs[1]);
    // common parts
    const mwSize *dim_Z = mxGetDimensions(prhs[0]);
    const mwSize *dim_para = mxGetDimensions(prhs[1]);
    int Nz = static_cast<int>(dim_Z[0] * dim_Z[1]);
    int N_para = static_cast<int>(dim_para[0] * dim_para[1]);
    if (nrhs == 2)
    {
        if (mxIsSingle(prhs[0]) && mxIsSingle(prhs[1])){
            // single precision
            if (dim_para[0] % 7 == 0)
            {
                float * Z = (float *)mxGetPr(prhs[0]);
                float * para = (float *)mxGetPr(prhs[1]);
                {
                    // get dose
                    const mwSize size[2]{dim_Z[0], dim_Z[1]};
                    plhs[0] = mxCreateNumericArray(2, size, mxSINGLE_CLASS, mxREAL);
                    float *idd_o_ptr;
                    idd_o_ptr = (float *)mxGetPr(plhs[0]);
                    IDD_array_v2(Z, idd_o_ptr, para, Nz, N_para);
                }
            }
            else if (dim_para[0] % 4 == 0)
            {
                float * Z = (float *)mxGetPr(prhs[0]);
                float * para = (float *)mxGetPr(prhs[1]);
                {
                    // get dose
                    const mwSize size[2]{dim_Z[0], dim_Z[1]};
                    plhs[0] = mxCreateNumericArray(2, size, mxSINGLE_CLASS, mxREAL);
                    float *idd_o_ptr;
                    idd_o_ptr = (float *)mxGetPr(plhs[0]);
                    IDD_array_v1(Z, idd_o_ptr, para, Nz, N_para);
                }
            }
            else{
                mexPrintf("mod(N_para(%d),7) != 0 or mod(N_para(%d),4) != 0!\n", N_para);
                mexErrMsgTxt("Parameter Size not match!\n");
            }
        }
        else if (mxIsDouble(prhs[0]) && mxIsDouble(prhs[1]))
        {
            // double precision
            mexPrintf("double pricision\n");
            if (dim_para[0] % 7 == 0)
            {
                double * Z = (double *)mxGetPr(prhs[0]);
                double * para = (double *)mxGetPr(prhs[1]);
                {
                    // get dose
                    const mwSize size[2]{dim_Z[0], dim_Z[1]};
                    plhs[0] = mxCreateNumericArray(2, size, mxDOUBLE_CLASS, mxREAL);
                    double *idd_o_ptr;
                    idd_o_ptr = (double *)mxGetPr(plhs[0]);
                    IDD_array_v2(Z, idd_o_ptr, para, Nz, N_para);
                }
            }
            else if (dim_para[0] % 4 == 0)
            {
                double * Z = (double *)mxGetPr(prhs[0]);
                double * para = (double *)mxGetPr(prhs[1]);
                {
                    // get dose
                    const mwSize size[2]{dim_Z[0], dim_Z[1]};
                    plhs[0] = mxCreateNumericArray(2, size, mxDOUBLE_CLASS, mxREAL);
                    double *idd_o_ptr;
                    idd_o_ptr = (double *)mxGetPr(plhs[0]);
                    IDD_array_v1(Z, idd_o_ptr, para, Nz, N_para);
                }
            }
            else{
                mexPrintf("mod(N_para(%d),7) != 0 or mod(N_para(%d),4) != 0!\n", N_para);
                mexErrMsgTxt("Parameter Size not match!\n");
            }

        }
        else
        {
            mexErrMsgTxt("Inconsistency float point input data type!, Should be all double or single\n");
        }
    }
    else{
        mexErrMsgTxt("Wrong number of input arguments!\n");
    }
}