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
    return A * c * exp(-0.5 * pow((x - mu) / sigma, 2.0));
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
    if ((sigma < 1e-10))
    {
        return 0.0;
    }
    else
    {
        T half_1_sigma2 = 0.5 * pow(1.0 / sigma, 2.0);
        return A * M_1_PI * half_1_sigma2 * exp(-half_1_sigma2 * (pow((x - mux), 2.0) + pow((y - muy), 2.0)));
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
void dose3d_2(std::vector<T> X, std::vector<T> Y, T* para, T* dose3d, int Nx, int Ny, int Nz)
{
#pragma omp parallel for firstprivate(X,Y,Nx,Ny,Nz)
    for (int nz = 0; nz < Nz; ++nz)
    {
        T A1 = para[nz * 8];
        T mux1 = para[nz * 8 + 1];
        T muy1 = para[nz * 8 + 2];
        T sigma1 = para[nz * 8 + 3];
        T A2 = para[nz * 8 + 4];
        T mux2 = para[nz * 8 + 5];
        T muy2 = para[nz * 8 + 6];
        T sigma2 = para[nz * 8 + 7];
        for (int ny = 0; ny < Ny; ++ny)
        {
            T y{ Y[ny] };
            for (int nx = 0; nx < Nx; ++nx)
            {
                T x{ X[nx] };
                int idx3d{ nx + ny * Nx + nz * (Nx * Ny) };
                dose3d[idx3d] = gauss2d_2(x, y, A1, mux1, muy1, sigma1, A2, mux2, muy2, sigma2);
            }
        }
    }
}
template<class T>
void dose3d_N_iso(std::vector<T> X, std::vector<T> Y, std::vector<T>  para, T* dose3d, int Nx, int Ny, int Nz, int N_gaussian)
{
#pragma omp parallel for firstprivate(X,Y,Nx,Ny,Nz,para,N_gaussian)
    for (int nz = 0; nz < Nz; ++nz)
    {
        T A{};
        for (int ng = 0; ng < N_gaussian; ++ng)
        {
            A += para[nz * N_gaussian * 4 + 4 * ng];
        }
        if (A > 1e-10)
        {
            for (int ny = 0; ny < Ny; ++ny)
            {
                T y{ Y[ny] };
                for (int nx = 0; nx < Nx; ++nx)
                {
                    T x{ X[nx] };
                    int idx3d{ nx + ny * Nx + nz * (Nx * Ny) };
                    T temp{};
                    for (int ng = 0; ng < N_gaussian; ++ng)
                    {
                        T A1 = para[nz * N_gaussian * 4 + 4 * ng];
                        T mux1 = para[nz * N_gaussian * 4 + 4 * ng + 1];
                        T muy1 = para[nz * N_gaussian * 4 + 4 * ng + 2];
                        T sigma1 = para[nz * N_gaussian * 4 + 4 * ng + 3];
                        temp += gauss2d(x, y, A1, mux1, muy1, sigma1);
                    }
                    dose3d[idx3d] = temp;
                }
            }
        }
    }
}

template<class T>
void dose3d_N(std::vector<T> X, std::vector<T> Y, std::vector<T>  para, T* dose3d, int Nx, int Ny, int Nz, int N_gaussian)
{
#pragma omp parallel for firstprivate(X,Y,Nx,Ny,Nz,para,N_gaussian)
    for (int nz = 0; nz < Nz; ++nz)
    {
        T A{};
        for (int ng = 0; ng < N_gaussian; ++ng)
        {
            A += para[nz * N_gaussian * 4 + 4 * ng];
        }
        if (A > 1e-10)
        {
            for (int ny = 0; ny < Ny; ++ny)
            {
                T y{ Y[ny] };
                for (int nx = 0; nx < Nx; ++nx)
                {
                    T x{ X[nx] };
                    int idx3d{ nx + ny * Nx + nz * (Nx * Ny) };
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

template void Gauss2d::interface(std::vector<double> X, std::vector<double> Y, std::vector<double> para, double* dose3d_ptr, int Nx, int Ny, int Nz, int N_para, int N_gaussian);
template void Gauss2d::interface(std::vector<float> X, std::vector<float> Y, std::vector<float> para, float* dose3d_ptr, int Nx, int Ny, int Nz, int N_para, int N_gaussian);



