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

template<class T>
T dfdx(T *p, T *q, T x){
    // derivative of f
    // part1 = (5p0*x^4 + 4p1*x^3 + 3p2*x^2 + 2p3*x^1 + p4)/(x^4+q0*x^3+q1*x^2+q2*x+q3)
    T part1 = (T(5) * p[0] + x * (T(4) * p[1] + x * (T(3) * p[2] + x * (T(2) * p[3] + x * T(1) * p[4])))) / (x * (x * (x * (x + q[0]) + q[1]) + q[2]) + q[3]);
    // part2 = (p0*x^5 + p1*x^4 + p2*x^3 + p3*x^2 + p4*x)*(4x^3+3q0*x^2+2q1*x+q2)/(x^4+q0*x^3+q1*x^2+q2*x+q3)^2
    T part2 = (p[0] + x * (p[1] + x * (p[2] + x * (p[3] + x * p[4])))) * (x * (x * (T(4) * x + T(3) * q[0]) + T(2) * q[1]) + q[2]) / ((x * (x * (x * (x + q[0]) + q[1]) + q[2]) + q[3]) * (x * (x * (x * (x + q[0]) + q[1]) + q[2]) + q[3]));
    return part1 - part2;
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
float df0565(float x){
    // derivative of parabolic cylinder function D(-0.565, x)*exp(-x^2/4)
    if (x < 0){
        return 1.230286920530541f*std::exp(f(p_neg565_f, q_neg565_f, x)) * dfdx(p_neg565_f, q_neg565_f, x);
    }
    else{
        return 1.230286920530541f*std::exp(f(p_pos565_f, q_pos565_f, x) - x*x/2.0f) * dfdx(p_pos565_f, q_pos565_f, x);
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
double df0565(double x){
    // derivative of parabolic cylinder function D(-0.565, x)*exp(-x^2/4)
    if (x < 0){
        return 1.230286920530541*std::exp(f(p_neg565, q_neg565, x)) * dfdx(p_neg565, q_neg565, x);
    }
    else{
        return 1.230286920530541*std::exp(f(p_pos565, q_pos565, x) - x*x/2.0) * dfdx(p_pos565, q_pos565, x);
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
float df1565(float x){
    // derivative of parabolic cylinder function D(-1.565, x)*exp(-x^2/4)
    if (x < 0){
        return 1.144555647104346f*std::exp(f(p_neg1565_f, q_neg1565_f, x)) * dfdx(p_neg1565_f, q_neg1565_f, x);
    }
    else{
        return 1.144555647104346f*std::exp(f(p_pos1565_f, q_pos1565_f, x) - x*x/2.0f) * dfdx(p_pos1565_f, q_pos1565_f, x);
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
double df1565(double x){
    // derivative of parabolic cylinder function D(-1.565, x)*exp(-x^2/4)
    if (x < 0){
        return 1.144555647104346*std::exp(f(p_neg1565, q_neg1565, x)) * dfdx(p_neg1565, q_neg1565, x);
    }
    else{
        return 1.144555647104346*std::exp(f(p_pos1565, q_pos1565, x) - x*x/2.0) * dfdx(p_pos1565, q_pos1565, x);
    }
}

float B2(float z, float r, float s, float e, float p){
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
void dB2(float z, float r, float s, float e, float p, float *d){
    // derivative of B2 respect to r, s, e, p
    // B2 = (p / r1012) * s^0.565 * ((11.26 / s) * f0565(-x) + (0.157 + 11.26 * e / r) * f1565(-x));
    // B2 = p * A * (part1 + part2);
    // dB2_dr 
    float A = 1.0f/(1.0f + 0.012f * r);
    float dA_dr = -0.012f * A * A;
    float x = (r - z) / s;
    float dx_dr = 1.0f / s;
    float dx_ds = (z - r) / (s * s);
    float df0565_dx = -df0565(-x);
    float df1565_dx = -df1565(-x);
    float df0565_dr = df0565_dx * dx_dr;
    float df1565_dr = df1565_dx * dx_dr;
    float df0565_ds = df0565_dx * dx_ds;
    float df1565_ds = df1565_dx * dx_ds;
    // part1  = s^-0.435 * 11.26 * f0565(-x);
    float part1 = std::pow(s, -0.435f) * 11.26f * f0565(-x);
    float dpart1_dr = std::pow(s, -0.435f) * 11.26f * df0565_dr;
    float dpart1_ds = -0.435f * std::pow(s, -1.435f) * 11.26f * f0565(-x) - std::pow(s, -0.435f) * 11.26f * df0565_ds;
    // part2 = s^0.565 * (0.157 + 11.26 * e / r) * f1565(-x);
    float part2 = std::pow(s, 0.565f) * (0.157f + 11.26f * e / r) * f1565(-x);
    float dpart2_dr = std::pow(s, 0.565f) * 0.157f * df1565_dr - std::pow(s, 0.565f) * 11.26f * e / (r * r) * f1565(-x) + std::pow(s, 0.565f) * 11.26f * e / r * df1565_dr;
    float dpart2_ds = 0.565f * std::pow(s, -0.435f) * (0.157f + 11.26f * e / r) * f1565(-x) + std::pow(s, 0.565f) * (0.157f + 11.26f * e / r) * df1565_ds;
    float dpart2_de = std::pow(s, 0.565f) * 11.26f / r * f1565(-x);

    // dB2_dr = p * A * (dpart1_dr + dpart2_dr) + dA_dr * p * (part1 + part2);
    d[0] = p * A * (dpart1_dr + dpart2_dr) + dA_dr * p * (part1 + part2);
    // dB2_ds = p * A * (dpart1_ds + dpart2_ds);
    d[1] = p * A * (dpart1_ds + dpart2_ds);
    // dB2_de = p * A * dpart2_de;
    d[2] = p * A * dpart2_de;
    // dB2_dp = A * (part1 + part2);
    d[3] = A * (part1 + part2);

}
double B2(double z, double r, double s, double e, double p){
    double r1012 = (1.0 + 0.012 * r);
    double x = (r - z) / s;
    // return (p / r1012) * std::pow(s, 0.565) * ((11.26 / s) * f0565(-x) + (0.157 + 11.26 * e / r) * f1565(-x));
    return (p / r1012) * std::pow(s, 0.565) * ((11.26 / s) * f0565(-x) + (0.157 + 11.26 * e / r) * f1565(-x));
}
void dB2(double z, double r, double s, double e, double p, double *d){
    // derivative of B2 respect to r, s, e, p
    // B2 = (p / r1012) * s^0.565 * ((11.26 / s) * f0565(-x) + (0.157 + 11.26 * e / r) * f1565(-x));
    // B2 = p * A * (part1 + part2);
    // dB2_dr 
    double A = 1.0/(1.0 + 0.012 * r);
    double dA_dr = -0.012 * A * A;
    double x = (r - z) / s;
    double dx_dr = 1.0 / s;
    double dx_ds = (z - r) / (s * s);
    double df0565_dx = -df0565(-x);
    double df1565_dx = -df1565(-x);
    double df0565_dr = df0565_dx * dx_dr;
    double df1565_dr = df1565_dx * dx_dr;
    double df0565_ds = df0565_dx * dx_ds;
    double df1565_ds = df1565_dx * dx_ds;
    // part1  = s^-0.435 * 11.26 * f0565(-x);
    double part1 = std::pow(s, -0.435) * 11.26 * f0565(-x);
    double dpart1_dr = std::pow(s, -0.435) * 11.26 * df0565_dr;
    double dpart1_ds = -0.435 * std::pow(s, -1.435) * 11.26 * f0565(-x) - std::pow(s, -0.435) * 11.26 * df0565_ds;
    // part2 = s^0.565 * (0.157 + 11.26 * e / r) * f1565(-x);
    double part2 = std::pow(s, 0.565) * (0.157 + 11.26 * e / r) * f1565(-x);
    double dpart2_dr = std::pow(s, 0.565) * 0.157 * df1565_dr - std::pow(s, 0.565) * 11.26 * e / (r * r) * f1565(-x) + std::pow(s, 0.565) * 11.26 * e / r * df1565_dr;
    double dpart2_ds = 0.565 * std::pow(s, -0.435) * (0.157 + 11.26 * e / r) * f1565(-x) + std::pow(s, 0.565) * (0.157 + 11.26 * e / r) * df1565_ds;
    double dpart2_de = std::pow(s, 0.565) * 11.26 / r * f1565(-x);

    // dB2_dr = p * A * (dpart1_dr + dpart2_dr) + dA_dr * p * (part1 + part2);
    d[0] = p * A * (dpart1_dr + dpart2_dr) + dA_dr * p * (part1 + part2);
    // dB2_ds = p * A * (dpart1_ds + dpart2_ds);
    d[1] = p * A * (dpart1_ds + dpart2_ds);
    // dB2_de = p * A * dpart2_de;
    d[2] = p * A * dpart2_de;
    // dB2_dp = A * (part1 + part2);
    d[3] = A * (part1 + part2);

}

float B1(float z, float r, float s, float e, float p)
{
    // z : depth
    // r : range R80 (80% of the maximum dose)
    // s : sigma (energy spread)
    // e : epsilon
    // p : phi (fluence)
    float r1012 = (1.0f + 0.012f * r);
    float r_z = r - z;
    return (p / r1012) * (17.93f * std::pow(r_z, -0.435f) + (0.444f + 31.7f * e / r) * std::pow(r_z, 0.565f));
    // return (p / r1012) * (17.974454939899370f * std::pow(r_z, -0.435f) + (0.441595774885904f+31.671136466339373f * e / r) * std::pow(r_z, 0.565f));
}
void dB1(float z, float r, float s, float e, float p, float *d){
    // derivative of B1 respect to r, s, e, p
    // B1 = p * (1/(1+0.012*r)) * (17.93 * (r-z)^(-0.435) + (0.444 + 31.7*e/r) * (r-z)^0.565);
    // B1 = p * A * (part1 + part2);
    float A = 1.0f/(1.0f + 0.012f * r);
    float dA_dr = -0.012f * A * A;
    float r_z = r - z;
    float dr_z_dr = 1.0f;

    float part1 = 17.93f * std::pow(r_z, -0.435f);
    float dpart1_dr = -0.435f * 17.93f * std::pow(r_z, -1.435f) * dr_z_dr;
    float part2 = (0.444f + 31.7f * e / r) * std::pow(r_z, 0.565f);
    float dpart2_dr = 0.565f * 0.444f * std::pow(r_z, -0.435f) * dr_z_dr - 31.7f * (e / (r * r)) * std::pow(r_z, 0.565f) + 31.7f * e / r * 0.565f * std::pow(r_z, -0.435f) * dr_z_dr;
    float dpart2_de = 31.7f / r * std::pow(r_z, 0.565f);

    // dB1_dr = p * A * (dpart1_dr + dpart2_dr) + dA_dr * p * (part1 + part2);
    d[0] = p * A * (dpart1_dr + dpart2_dr) + dA_dr * p * (part1 + part2);
    // dB1_ds = 0;
    d[1] = 0;
    // dB1_de = p * A * dpart2_de;
    d[2] = p * A * dpart2_de;
    // dB1_dp = A * (part1 + part2);
    d[3] = A * (part1 + part2);
}
double B1(double z, double r, double s, double e, double p)
{
    double r1012 = (1.0 + 0.012 * r);
    double r_z = r - z;
    return (p / r1012) * (17.93 * std::pow(r_z, -0.435) + (0.444 + 31.7 * e / r) * std::pow(r_z, 0.565));
    // return (p / r1012) * (17.974454939899370 * std::pow(r_z, -0.435) + (0.441595774885904+31.671136466339373 * e / r) * std::pow(r_z, 0.565));
}
void dB1(double z, double r, double s, double e, double p, double *d){
    // derivative of B1 respect to r, s, e, p
    // B1 = p * (1/(1+0.012*r)) * (17.93 * (r-z)^(-0.435) + (0.444 + 31.7*e/r) * (r-z)^0.565);
    // B1 = p * A * (part1 + part2);
    double A = 1.0/(1.0 + 0.012 * r);
    double dA_dr = -0.012 * A * A;
    double r_z = r - z;
    double dr_z_dr = 1.0;

    double part1 = 17.93 * std::pow(r_z, -0.435);
    double dpart1_dr = -0.435 * 17.93 * std::pow(r_z, -1.435) * dr_z_dr;
    double part2 = (0.444 + 31.7 * e / r) * std::pow(r_z, 0.565);
    double dpart2_dr = 0.565 * 0.444 * std::pow(r_z, -0.435) * dr_z_dr - 31.7 * (e / (r * r)) * std::pow(r_z, 0.565) + 31.7 * e / r * 0.565 * std::pow(r_z, -0.435) * dr_z_dr;
    double dpart2_de = 31.7 / r * std::pow(r_z, 0.565);

    // dB1_dr = p * A * (dpart1_dr + dpart2_dr) + dA_dr * p * (part1 + part2);
    d[0] = p * A * (dpart1_dr + dpart2_dr) + dA_dr * p * (part1 + part2);
    // dB1_ds = 0;
    d[1] = 0;
    // dB1_de = p * A * dpart2_de;
    d[2] = p * A * dpart2_de;
    // dB1_dp = A * (part1 + part2);
    d[3] = A * (part1 + part2);
}

template<class T>
T D4(T depth, T R0, T sigma, T epsilon, T Phi)
{
    // reference : Bortfeld float. An analytical approximation of the Bragg curve for therapeutic proton beams. Med Phys. 1997; 24(12) : 2024 - 2033. doi : 10.1118 / 1.598116
    T dose{};
    // analytical bragg curve
    if (depth < R0 - T(10.0) * sigma)
    {
        dose = B1(depth, R0, sigma, epsilon, Phi);
    }
    else if (R0 - T(10.0) * sigma <= depth && depth <= R0 + T(5.0) * sigma)
    {
        dose = B2(depth, R0, sigma, epsilon, Phi);
    }
    else
    {
        dose = (T)0.0;
    }
    return dose;
}
template<class T>
void dD4(T depth, T R0, T sigma, T epsilon, T Phi, T *d){
    // derivative of D4 respect to R0, sigma, epsilon, Phi
    if (depth < R0 - T(10.0) * sigma){
        dB1(depth, R0, sigma, epsilon, Phi, d);
    }
    else if (R0 - T(10.0) * sigma <= depth && depth <= R0 + T(5.0) * sigma){
        dB2(depth, R0, sigma, epsilon, Phi, d);
    }
    else{
        d[0] = 0;
        d[1] = 0;
        d[2] = 0;
        d[3] = 0;
    }
}

template<class T>
T D7(T depth, T R0, T sigma, T epsilon, T Phi, T a, T b, T c)
{
    // reference : Zhang, Xiaodong, et al. "Parameterization of multiple Bragg curves for scanning proton beams using simultaneous fitting of multiple curves." Physics in Medicine & Biology 56.24 (2011): 7725.
    // fixed discontinuity of Zhang's formula
    T dose{};
    T edge = R0 - T(10.0) * sigma;
    // analytical bragg curve
    if (depth < edge)
    {
        dose = B1(depth, R0, sigma, epsilon, Phi) + a * std::pow(edge - depth, 3.0) + b * std::pow(edge - depth, 2.0) + c * (edge - depth);
    }
    else if (edge <= depth && depth <= R0 + T(5.0) * sigma)
    {
        dose = B2(depth, R0, sigma, epsilon, Phi);
    }
    else
    {
        dose = (T)0.0;
    }
    return dose;
}
template<class T>
void dD7(T depth, T R0, T sigma, T epsilon, T Phi, T a, T b, T c, T *derivative){
    T edge = R0 - T(10.0) * sigma;
    T dedge_dr = 1.0;
    T dedge_ds = -T(10.0);
    if (depth < R0 - T(10.0) * sigma){
        dB1(depth, R0, sigma, epsilon, Phi, derivative);
        derivative[0] += T(3.0) * a * std::pow(edge - depth, 2.0) * dedge_dr + T(2.0) * b * (edge - depth) * dedge_dr + c * dedge_dr;
        derivative[1] += T(3.0) * a * std::pow(edge - depth, 2.0) * dedge_ds + T(2.0) * b * (edge - depth) * dedge_ds + c * dedge_ds;
        derivative[4] = std::pow(R0 - depth, 3.0);
        derivative[5] = std::pow(R0 - depth, 2.0);
        derivative[6] = R0 - depth;
    }
    else if (edge <= depth && depth <= R0 + T(5.0) * sigma){
        dB2(depth, R0, sigma, epsilon, Phi, derivative);
        derivative[4] = 0;
        derivative[5] = 0;
        derivative[6] = 0;
    }
    else{
        derivative[0] = 0;
        derivative[1] = 0;
        derivative[2] = 0;
        derivative[3] = 0;
        derivative[4] = 0;
        derivative[5] = 0;
        derivative[6] = 0;
    }
}

template<class T>
void D7_multi(T* depth, T* dose_o, T* para_i, int N_depth, int para_size)
{
    // use multiple-bf mixture
#pragma omp parallel for
    for (int i = 0; i < N_depth; ++i)
    {
        T temp{};
        for (int k = 0; k < para_size; k = k + 7)// multi-bf mixture
        {
            temp += D7(depth[i], para_i[k], para_i[k + 1], para_i[k + 2], para_i[k + 3], para_i[k + 4], para_i[k + 5], para_i[k + 6]);
        }
        dose_o[i] = temp;
    }
}
template<class T>
void dD7_multi(T* depth, T* derivative, T* para_i, int N_depth, int para_size)
{
    // get derivative of multiple-bf mixture
    // derivative size = para_size * N_depth
#pragma omp parallel for
    for (int i = 0; i < N_depth; ++i)
    {
        for (int k = 0; k < para_size; k = k + 7)// multi-bf mixture
        {
            T derivative_bf[7]{};
            dD7(depth[i], para_i[k], para_i[k + 1], para_i[k + 2], para_i[k + 3], para_i[k + 4], para_i[k + 5], para_i[k + 6], derivative_bf);
            // copy to derivative with memory copy
            memcpy(derivative + k + i * para_size, derivative_bf, 7 * sizeof(T));
        }
    }
}
template<class T>
void D4_multi(T* depth, T* dose_o, T* para_i, int N_depth, int para_size)
{
    // use multiple-bf mixture
#pragma omp parallel for
    for (int i = 0; i < N_depth; ++i)
    {
        T temp{};
        for (int k = 0; k < para_size; k = k + 4)// multi-bf mixture
        {
            temp += D4(depth[i], para_i[k], para_i[k + 1], para_i[k + 2], para_i[k + 3]);
        }
        dose_o[i] = temp;
    }
}
template<class T>
void dD4_multi(T* depth, T* derivative, T* para_i, int N_depth, int para_size)
{
    // get derivative of multiple-bf mixture
    // derivative size = para_size * N_depth
#pragma omp parallel for
    for (int i = 0; i < N_depth; ++i)
    {
        for (int k = 0; k < para_size; k = k + 4)// multi-bf mixture
        {
            T derivative_bf[4]{};
            dD4(depth[i], para_i[k], para_i[k + 1], para_i[k + 2], para_i[k + 3], derivative_bf);
            // copy to derivative with memory copy
            memcpy(derivative + k + i * para_size, derivative_bf, 4 * sizeof(T));
        }
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
    int N_bp = static_cast<int>(dim_para[1]);// number of Bragg peak (bp). para size (4,N_bp) or (7,N_bp)
    int N_para = static_cast<int>(dim_para[0] * dim_para[1]);
    if (nrhs == 2)
    { // input 2 arguments, depth and parameter
        if (mxIsSingle(prhs[0]) && mxIsSingle(prhs[1])){
            // single precision
            if (dim_para[0] == 7)
            {
                float * Z = (float *)mxGetPr(prhs[0]);
                float * para = (float *)mxGetPr(prhs[1]);
                {
                    // get dose
                    const mwSize size[2]{dim_Z[0], dim_Z[1]};
                    plhs[0] = mxCreateNumericArray(2, size, mxSINGLE_CLASS, mxREAL);
                    float *idd_o_ptr;
                    idd_o_ptr = (float *)mxGetPr(plhs[0]);
                    D7_multi(Z, idd_o_ptr, para, Nz, N_para);
                }
            }
            else if (dim_para[0] == 4)
            {
                float * Z = (float *)mxGetPr(prhs[0]);
                float * para = (float *)mxGetPr(prhs[1]);
                {
                    // get dose
                    const mwSize size[2]{dim_Z[0], dim_Z[1]};
                    plhs[0] = mxCreateNumericArray(2, size, mxSINGLE_CLASS, mxREAL);
                    float *idd_o_ptr;
                    idd_o_ptr = (float *)mxGetPr(plhs[0]);
                    D4_multi(Z, idd_o_ptr, para, Nz, N_para);
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
            // mexPrintf("double pricision\n");
            if (dim_para[0] == 7)
            {
                double * Z = (double *)mxGetPr(prhs[0]);
                double * para = (double *)mxGetPr(prhs[1]);
                {
                    // get dose
                    const mwSize size[2]{dim_Z[0], dim_Z[1]};
                    plhs[0] = mxCreateNumericArray(2, size, mxDOUBLE_CLASS, mxREAL);
                    double *idd_o_ptr;
                    idd_o_ptr = (double *)mxGetPr(plhs[0]);
                    D7_multi(Z, idd_o_ptr, para, Nz, N_para);
                }
            }
            else if (dim_para[0] == 4)
            {
                double * Z = (double *)mxGetPr(prhs[0]);
                double * para = (double *)mxGetPr(prhs[1]);
                {
                    // get dose
                    const mwSize size[2]{dim_Z[0], dim_Z[1]};
                    plhs[0] = mxCreateNumericArray(2, size, mxDOUBLE_CLASS, mxREAL);
                    double *idd_o_ptr;
                    idd_o_ptr = (double *)mxGetPr(plhs[0]);
                    D4_multi(Z, idd_o_ptr, para, Nz, N_para);
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
    else if(nrhs ==3){
        // input 3 arguments, depth, parameter, and derivative
        if (mxIsSingle(prhs[0]) && mxIsSingle(prhs[1])){
            // single precision
            if (dim_para[0] == 7)
            {
                float * Z = (float *)mxGetPr(prhs[0]);
                float * para = (float *)mxGetPr(prhs[1]);
                {
                    // get dose
                    const mwSize size[2]{mwSize(N_para), mwSize(Nz)};
                    plhs[0] = mxCreateNumericArray(2, size, mxSINGLE_CLASS, mxREAL);
                    float *derivative_o_ptr;
                    derivative_o_ptr = (float *)mxGetPr(plhs[0]);
                    dD7_multi(Z, derivative_o_ptr, para, Nz, N_para);
                }
            }
            else if (dim_para[0] == 4)
            {
                float * Z = (float *)mxGetPr(prhs[0]);
                float * para = (float *)mxGetPr(prhs[1]);
                {
                    // get dose
                    const mwSize size[2]{mwSize(N_para), mwSize(Nz)};
                    plhs[0] = mxCreateNumericArray(2, size, mxSINGLE_CLASS, mxREAL);
                    float *derivative_o_ptr;
                    derivative_o_ptr = (float *)mxGetPr(plhs[0]);
                    dD4_multi(Z, derivative_o_ptr, para, Nz, N_para);
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
            // mexPrintf("double pricision\n");
            if (dim_para[0] == 7)
            {
                double * Z = (double *)mxGetPr(prhs[0]);
                double * para = (double *)mxGetPr(prhs[1]);
                {
                    // get dose
                    const mwSize size[2]{mwSize(N_para), mwSize(Nz)};
                    plhs[0] = mxCreateNumericArray(2, size, mxDOUBLE_CLASS, mxREAL);
                    double *derivative_o_ptr;
                    derivative_o_ptr = (double *)mxGetPr(plhs[0]);
                    dD7_multi(Z, derivative_o_ptr, para, Nz, N_para);
                }
            }
            else if (dim_para[0] == 4)
            {
                double * Z = (double *)mxGetPr(prhs[0]);
                double * para = (double *)mxGetPr(prhs[1]);
                {
                    // get dose
                    const mwSize size[2]{mwSize(N_para), mwSize(Nz)};
                    plhs[0] = mxCreateNumericArray(2, size, mxDOUBLE_CLASS, mxREAL);
                    double *derivative_o_ptr;
                    derivative_o_ptr = (double *)mxGetPr(plhs[0]);
                    dD4_multi(Z, derivative_o_ptr, para, Nz, N_para);
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