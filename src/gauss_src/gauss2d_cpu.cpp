// // dose3d.cpp : This file contains the 'mexFunction' function. Program execution begins and ends there.
// //
 #define _USE_MATH_DEFINES
// #include <array>
// #include <vector>
 #include <cmath>
// #include <iostream>
// #include <omp.h>
// // #include "mex.h"
#include "gauss2d.h"


template<class T>
inline T gauss1d(T x, T A, T mu, T sigma)
{
    T c{ M_SQRT1_2 * sqrt(M_1_PI) / sigma };
    return A * c * exp(-0.5 * pow((x - mu) / sigma, 2));
}
template<class T>
void gauss1DGradient(std::vector<T> X, std::vector<T>  para, T* ouput, int Nx, int N_gaussian)
{
#pragma omp parallel for firstprivate(X,para,Nx,N_gaussian)
    for (int nx = 0; nx < Nx; ++nx)
    {
        T x{ X[nx] };
        for (int ng = 0; ng < N_gaussian; ++ng)
        {
            T A = para[3 * ng];
            T mu = para[3 * ng + 1];
            T sigma = para[3 * ng + 2];
            T w1 = A * ((x - mu) / (sigma * sigma));
            T w2 = A * (mu * mu - sigma * sigma + x * x - 2 * mu * x) / (sigma * sigma * sigma);
            T G = gauss1d(x, 1.0, mu, sigma);
            ouput[3 * N_gaussian * nx + 3 * ng] = G;// dGdA
            ouput[3 * N_gaussian * nx + 3 * ng + 1] = w1 * G;// dGdmu
            ouput[3 * N_gaussian * nx + 3 * ng + 2] = w2 * G;// dGdsigma
        }
    }
}
template<class T>
inline T gauss2d(T x, T y, T A, T mux, T muy, T sigma)
{
    if ((sigma < 1e-7) | (A < 1e-7))
    {
        return 0.0;
    }
    else
    {
        T one_sigma_sqr = T(1.0) / (sigma * sigma);
        T half_1_sigma2 = T(0.5) * one_sigma_sqr;
        return A * M_1_PI * half_1_sigma2 * exp(-half_1_sigma2 * (pow((x - mux), 2) + pow((y - muy), 2)));
    }
}


template<class T>
inline T dot2(T v1, T v2, T u1, T u2)
{
    return v1 * u1 + v2 * u2;
}

template<class T>
inline T mvn2d(T x, T y, T A, T mux, T muy, T sigma1, T sigma2, T beta)
{
    if ((sigma1 < 1e-7) | (sigma2 < 1e-7) | (A < 1e-7))
    {
        // float number effective decimal digits is 7
        return 0.0;
    }
    else
    {
        // Rotation matrix R, coordinate X = [x,y]-[mux,muy], Y = RX
        T sinb = std::sin(beta);
        T cosb = std::cos(beta);
        T Y1 = cosb * (x - mux) - sinb * (y - muy);// beta is angle in rad
        T Y2 = sinb * (x - mux) + cosb * (y - muy);
        // S = Y' * M^-1 * Y
		T S1 = (Y1 / sigma1) * (Y1 / sigma1);
        T S2 = (Y2 / sigma2) * (Y2 / sigma2);
        T exponant = -0.5 * (S1+S2);
        T scale = (A * 0.5 * M_1_PI) / (sigma1*sigma2);

        return scale * exp(exponant);
    }
}

template<class T>
void dose3d(std::vector<T> X, std::vector<T> Y, T* para, T* dose3d, int Nx, int Ny, int Nz)
{
#pragma omp parallel for firstprivate(X,Y,Nx,Ny,Nz)
    for (int nz = 0; nz < Nz; ++nz)
    {
        T A1 = para[nz * 4];
        T mux1 = para[nz * 4 + 1];
        T muy1 = para[nz * 4 + 2];
        T sigma1 = para[nz * 4 + 3];
        for (int ny = 0; ny < Ny; ++ny)
        {
            T y{ Y[ny] };
            for (int nx = 0; nx < Nx; ++nx)
            {
                T x{ X[nx] };
                int idx3d{ nx + ny * Nx + nz * (Nx * Ny) };
                dose3d[idx3d] = gauss2d(x, y, A1, mux1, muy1, sigma1);
            }
        }
    }
}

template<class T>
void dose3d_N_iso(std::vector<T> X, std::vector<T> Y, std::vector<T>  para, T* dose3d, int Nx, int Ny, int Nz, int N_gaussian)
{
#pragma omp parallel for firstprivate(X,Y,Nx,Ny,Nz,para,N_gaussian)

    for (int ny = 0; ny < Ny; ++ny)
    {
        T y{Y[ny]};
        for (int nx = 0; nx < Nx; ++nx)
        {
            T x{X[nx]};
            for (int nz = 0; nz < Nz; ++nz)
            {
                int idx3d{nx + ny * Nx + nz * (Nx * Ny)};
                T temp{};
                for (int ng = 0; ng < N_gaussian; ++ng)
                {
                    T A1 = para[nz * N_gaussian * 4 + 4 * ng];
                    T mux1 = para[nz * N_gaussian * 4 + 4 * ng + 1];
                    T muy1 = para[nz * N_gaussian * 4 + 4 * ng + 2];
                    T sigma1 = para[nz * N_gaussian * 4 + 4 * ng + 3];
                    T dist2d = (x-mux1)*(x-mux1) + (y-muy1)*(y-muy1);
                    if (dist2d < 16*sigma1*sigma1)
                    {
                        temp += gauss2d(x, y, A1, mux1, muy1, sigma1);
                    }
                }
                dose3d[idx3d] = temp;
            }
        }
    }
}
template<class T>
void dose3d_N_iso_Gradient(std::vector<T> X, std::vector<T> Y, std::vector<T>  para, T* output, int Nx, int Ny, int Nz, int N_gaussian)
{
    int64_t N_gauss_para = int64_t(Nz) * int64_t(N_gaussian) * int64_t(4) ;
    
#pragma omp parallel for firstprivate(X,Y,Nx,Ny,Nz,para,N_gaussian,N_gauss_para)

    for (int ny = 0; ny < Ny; ++ny)
    {
        T y{Y[ny]};
        for (int nx = 0; nx < Nx; ++nx)
        {
            T x{X[nx]};
            for (int nz = 0; nz < Nz; ++nz)
            {
                int64_t idx3d{nx + ny * Nx + nz * (Nx * Ny)};
                for (int ng = 0; ng < N_gaussian; ++ng)
                {
                    T A1 = para[nz * N_gaussian * 4 + 4 * ng];
                    T mux1 = para[nz * N_gaussian * 4 + 4 * ng + 1];
                    T muy1 = para[nz * N_gaussian * 4 + 4 * ng + 2];
                    T sigma1 = para[nz * N_gaussian * 4 + 4 * ng + 3];
                    T sigma_sqr = sigma1 * sigma1;
                    T x_mux = x - mux1;
                    T y_muy = y - muy1;
                    T dist2d = x_mux*x_mux + y_muy*y_muy;
                    if (dist2d < 16*sigma_sqr)
                    {
                        T GA = gauss2d(x, y, T(1.0), mux1, muy1, sigma1); // avoid G/A when A is small
                        T G = A1 * GA;
                        T One_sigma_sqr = 1.0 / (sigma_sqr);
                        T w = (x_mux * x_mux + y_muy * y_muy - 2.0 * sigma_sqr) / (sigma1 * sigma_sqr);
                        output[idx3d * N_gauss_para + nz * N_gaussian * 4 + 4 * ng] = GA;                      // dG/dA
                        output[idx3d * N_gauss_para + nz * N_gaussian * 4 + 4 * ng + 1] = x_mux * One_sigma_sqr * G; // dG/dmux
                        output[idx3d * N_gauss_para + nz * N_gaussian * 4 + 4 * ng + 2] = y_muy * One_sigma_sqr * G; // dG/dmuy
                        output[idx3d * N_gauss_para + nz * N_gaussian * 4 + 4 * ng + 3] = w * G;                     // dG/ds
                    }
                }
            }
        }
    }
}

template<class T>
void dose3d_N(std::vector<T> X, std::vector<T> Y, std::vector<T>  para, T* dose3d, int Nx, int Ny, int Nz, int N_gaussian)
{
#pragma omp parallel for firstprivate(X,Y,Nx,Ny,Nz,para,N_gaussian)

    for (int ny = 0; ny < Ny; ++ny)
    {
        T y{Y[ny]};
        for (int nx = 0; nx < Nx; ++nx)
        {
            T x{X[nx]};
            for (int nz = 0; nz < Nz; ++nz)
            {
                int idx3d{nx + ny * Nx + nz * (Nx * Ny)};
                T temp{};
                for (int ng = 0; ng < N_gaussian; ++ng)
                {
                    T A = para[nz * N_gaussian * 6 + 6 * ng];
                    T mux = para[nz * N_gaussian * 6 + 6 * ng + 1];
                    T muy = para[nz * N_gaussian * 6 + 6 * ng + 2];
                    T sigma1 = para[nz * N_gaussian * 6 + 6 * ng + 3];
                    T sigma2 = para[nz * N_gaussian * 6 + 6 * ng + 4];
                    T beta = para[nz * N_gaussian * 6 + 6 * ng + 5];
                    temp += mvn2d(x, y, A, mux, muy, sigma1, sigma2, beta);
                }
                dose3d[idx3d] = temp;
            }
        }
    }
}
template<class T>
void dose3d_N_Gradient(std::vector<T> X, std::vector<T> Y, std::vector<T>  para, T* output, int Nx, int Ny, int Nz, int N_gaussian)
{
    int64_t N_gauss_para = int64_t(Nz) * int64_t(N_gaussian) * int64_t(6) ;
#pragma omp parallel for firstprivate(X,Y,Nx,Ny,Nz,para,N_gaussian)
    for (int ny = 0; ny < Ny; ++ny)
    {
        T y{Y[ny]};
        for (int nx = 0; nx < Nx; ++nx)
        {
            T x{X[nx]};
            for (int nz = 0; nz < Nz; ++nz)
            {
                int idx3d{nx + ny * Nx + nz * (Nx * Ny)};
                for (int ng = 0; ng < N_gaussian; ++ng)
                {
                    T A = para[nz * N_gaussian * 6 + 6 * ng];
                    T mux = para[nz * N_gaussian * 6 + 6 * ng + 1];
                    T muy = para[nz * N_gaussian * 6 + 6 * ng + 2];
                    T sigma1 = para[nz * N_gaussian * 6 + 6 * ng + 3];
                    T sigma2 = para[nz * N_gaussian * 6 + 6 * ng + 4];
                    T beta = para[nz * N_gaussian * 6 + 6 * ng + 5];
                    T GA = mvn2d(x, y, T(1.0), mux, muy, sigma1, sigma2, beta);
                    T G = A * GA;
                    T sigma1_sqr = sigma1 * sigma1;
                    T sigma2_sqr = sigma2 * sigma2;
                    T One_sigma1_sqr = 1.0 / (sigma1_sqr);
                    T One_sigma2_sqr = 1.0 / (sigma2_sqr);
                    T x_mux = x - mux;
                    T y_muy = y - muy;
                    T cosb = cos(beta);
                    T sinb = sin(beta);
                    T Y1 = x_mux * cosb - y_muy * sinb;
                    T Y2 = x_mux * sinb + y_muy * cosb;
                    T S1 = (Y1 * Y1) / (sigma1_sqr);
                    T S2 = (Y2 * Y2) / (sigma2_sqr);
                    T Y1_sigma1_sqr = Y1 * One_sigma1_sqr;
                    T Y2_sigma2_sqr = Y2 * One_sigma2_sqr;
                    T w1 = (Y1_sigma1_sqr * cosb + Y2_sigma2_sqr * sinb);
                    T w2 = (-Y1_sigma1_sqr * sinb + Y2_sigma2_sqr * cosb);
                    T w3 = (S1 - 1.0) / (sigma1);
                    T w4 = (S2 - 1.0) / (sigma2);
                    T w5 = Y1 * Y2 * (One_sigma1_sqr - One_sigma2_sqr);
                    output[idx3d * N_gauss_para + nz * N_gaussian * 6 + 6 * ng] = GA;    // dG/dA
                    output[idx3d * N_gauss_para + nz * N_gaussian * 6 + 6 * ng + 1] = w1 * G; // dG/dmux
                    output[idx3d * N_gauss_para + nz * N_gaussian * 6 + 6 * ng + 2] = w2 * G; // dG/dmuy
                    output[idx3d * N_gauss_para + nz * N_gaussian * 6 + 6 * ng + 3] = w3 * G; // dG/ds1
                    output[idx3d * N_gauss_para + nz * N_gaussian * 6 + 6 * ng + 4] = w4 * G; // dG/ds2
                    output[idx3d * N_gauss_para + nz * N_gaussian * 6 + 6 * ng + 5] = w5 * G; // dG/db
                }
            }
        }
    }
}
template<class T>
void Gauss2d::interface(std::vector<T> X, std::vector<T> Y, std::vector<T> para, T* dose3d_ptr, int Nx, int Ny, int Nz, int N_para, int N_gaussian)
{
    if (N_para == Nz * N_gaussian * 6)
    {
        dose3d_N(X, Y, para, dose3d_ptr, Nx, Ny, Nz, N_gaussian);
    }
    else if (N_para == Nz * N_gaussian * 4)
    {
        dose3d_N_iso(X, Y, para, dose3d_ptr, Nx, Ny, Nz, N_gaussian);
    }
}

template<class T>
void Gauss2d::interface_gradient(std::vector<T> X, std::vector<T> Y, std::vector<T> para, T* grad_ptr, int Nx, int Ny, int Nz, int N_para, int N_gaussian)
{
    if (N_para == Nz * N_gaussian * 6)
    {
        dose3d_N_Gradient(X, Y, para, grad_ptr, Nx, Ny, Nz, N_gaussian);
    }
    else if (N_para == Nz * N_gaussian * 4)
    {
        dose3d_N_iso_Gradient(X, Y, para, grad_ptr, Nx, Ny, Nz, N_gaussian);
    }
}

template void Gauss2d::interface(std::vector<double> X, std::vector<double> Y, std::vector<double> para, double* dose3d_ptr, int Nx, int Ny, int Nz, int N_para, int N_gaussian);
template void Gauss2d::interface(std::vector<float> X, std::vector<float> Y, std::vector<float> para, float* dose3d_ptr, int Nx, int Ny, int Nz, int N_para, int N_gaussian);

template void Gauss2d::interface_gradient(std::vector<double> X, std::vector<double> Y, std::vector<double> para, double* dose3d_ptr, int Nx, int Ny, int Nz, int N_para, int N_gaussian);
template void Gauss2d::interface_gradient(std::vector<float> X, std::vector<float> Y, std::vector<float> para, float* dose3d_ptr, int Nx, int Ny, int Nz, int N_para, int N_gaussian);

