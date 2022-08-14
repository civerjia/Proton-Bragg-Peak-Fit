#define _USE_MATH_DEFINES
#include <algorithm>
#include <cmath>
#include <array>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <string>
#include <fstream>      // std::ofstream
#include <omp.h>

#include "bp.h"
#include "mex.h"

namespace Gauss2d
{
    template<class T>
    void interface(std::vector<T> X, std::vector<T> Y, std::vector<T> para, T* dose3d_ptr, int Nx, int Ny, int Nz, int N_para, int N_gaussian);

    template <class T>
    void linspace(T* L, T dmin, T dmax, int N)
    {
        for (int i = 0; i < N; ++i)
        {
            L[i] = dmin + T(i) * (dmax - dmin) / T(N - 1);
        }
    }
}

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

// template void Gauss2d::interface(std::vector<double> X, std::vector<double> Y, std::vector<double> para, double* dose3d_ptr, int Nx, int Ny, int Nz, int N_para, int N_gaussian);
// template void Gauss2d::interface(std::vector<float> X, std::vector<float> Y, std::vector<float> para, float* dose3d_ptr, int Nx, int Ny, int Nz, int N_para, int N_gaussian);

template<class T>
void updata_Ai(T* gauss_para, T* idd_o_ptr, int Nz, int N_gaussian, int Ng)
{
    // Ng : number of parameters in gaussian2d model , 4 or 6
    #pragma omp parallel for firstprivate(Nz,Ng,N_gaussian)
    for (int nz = 0; nz < Nz; nz = nz + 1)
	{
		T area = idd_o_ptr[nz];
		for (int i = 0; i < N_gaussian; ++i) 
		{
			// update Ai(z) = Ai(z) * IDD[z]
			gauss_para[nz * N_gaussian * Ng + i * Ng] = gauss_para[nz * N_gaussian * Ng + i * Ng] * area;
		}
	}
}

// disable matlab entry function
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    float *X;
    float *Y;
    float *Z;
    float *gauss_para;
    float *bf_para;
    X = (float*)mxGetPr(prhs[0]);
    Y = (float*)mxGetPr(prhs[1]);
    Z = (float*)mxGetPr(prhs[2]);

    gauss_para = (float*)mxGetPr(prhs[3]);
    bf_para = (float*)mxGetPr(prhs[4]);

    int N_gaussian = static_cast<int>(*mxGetPr(prhs[5]));
    int isGPU = static_cast<int>(*mxGetPr(prhs[6]));

    const mwSize *dim_X = mxGetDimensions(prhs[0]);
    const mwSize *dim_Y = mxGetDimensions(prhs[1]);
    const mwSize *dim_Z = mxGetDimensions(prhs[2]);

    const mwSize *dim_gauss_para = mxGetDimensions(prhs[3]);
    const mwSize *dim_bf_para = mxGetDimensions(prhs[4]);
    int Nx = static_cast<int>(dim_X[0]*dim_X[1]);
    int Ny = static_cast<int>(dim_Y[0]*dim_Y[1]);
    int Nz = static_cast<int>(dim_Z[0]*dim_Z[1]);
    int N_gauss_para = static_cast<int>(dim_gauss_para[0]*dim_gauss_para[1]);
    int N_bf_para = static_cast<int>(dim_bf_para[0]*dim_bf_para[1]);

    // gauss2d 
    std::vector<float> X_vec(Nx), Y_vec(Ny);
    std::copy(X, X + Nx, X_vec.begin());
    std::copy(Y, Y + Ny, Y_vec.begin());
    std::vector<float> gauss_para_vec(N_gauss_para);
    std::copy(gauss_para, gauss_para + N_gauss_para, gauss_para_vec.begin());
    // output 
    const mwSize size[3]{ mwSize(Nx), mwSize(Ny), mwSize(Nz) };
    plhs[0] = mxCreateNumericArray(3, size, mxSINGLE_CLASS, mxREAL);//mxDOUBLE_CLASS
    float* dose3d_ptr{};
    dose3d_ptr = (float*)mxGetPr(plhs[0]);
    
    // run calculation
    std::vector<float> idd_o(Nz);
    // step 1 get IDD
	BP::IDD_array_N( Z, idd_o.data(), bf_para, Nz, N_bf_para);
    if((N_gauss_para == Nz*N_gaussian*6) | (N_gauss_para == Nz*N_gaussian*4))
    {
        int Ng = N_gauss_para / (Nz * N_gaussian);
        // step 2 update weight Ai in Gauss2d model
	    updata_Ai(gauss_para_vec.data(), idd_o.data(), Nz, N_gaussian, Ng);
        // step 3 get 3D dose layer by layer
        Gauss2d::interface( X_vec, Y_vec, gauss_para_vec, dose3d_ptr, Nx, Ny, Nz, N_gauss_para, N_gaussian);// cpu version
    }
    else
    {
        mexPrintf("N_gaussian(%d),Nz(%d) not match para size(%d)!\n",N_gaussian,Nz,N_gauss_para);
    }
}