// //#define _USE_MATH_DEFINES
// #include "cuda_runtime.h"
// #include "device_launch_parameters.h"

// #include <thrust/host_vector.h>
// #include <thrust/device_vector.h>
// #include <thrust/generate.h>
// #include <thrust/sort.h>
// #include <thrust/copy.h>
// #include <algorithm>

// #include <cmath>
// #include <array>
// #include <vector>
// //#include <cmath>// cannot use in cu files
// #include <iostream>
// #include <stdio.h>
// // #include "mex.h"
// // #include <omp.h>
// #include <string>
// #include <fstream>      // std::ofstream
// #include <filesystem>
#include "gauss2d.h"

//https://stackoverflow.com/questions/37566987/cuda-atomicadd-for-doubles-definition-error
#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600
#else
__device__ double atomicAdd(double* a, double b) { return b; }
#endif

__constant__ const float const_1_PI_f = 0.318309886183791;
__constant__ const float const_1_2PI_f = 0.159154943091895;
__constant__ const double const_1_PI = 0.318309886183791;
__constant__ const double const_1_2PI = 0.159154943091895;


// 1/(2*pi), 2D gauss constant
__constant__ int constmem_Nsize[5];
// declare constant memory, Nx,Ny,Nz,N_gaussian,N_para
// max 16384 float numbers
// __constant__ float constmem_para_f[6*2*1024];// N_para = 6*N_gaussian*Nz < 6*2*512, N_gaussian = 2, Nz <=512


static void HandleError(cudaError_t err, const char *file, int line)
{
	if(err!= cudaSuccess)
	{
		printf("%s in %s at line %d\n",cudaGetErrorString(err),file,line);
		exit(EXIT_FAILURE);
	}
}
#define HANDLE_ERROR(err) (HandleError(err,__FILE__,__LINE__))


void set_constant_mem_Nsize(int Nx, int Ny, int Nz, int N_gaussian, int N_para)
{
    int cNsize[5] = { Nx,Ny,Nz,N_gaussian,N_para }; // copy host data to constant memory
    //cudaError_t mem_err;
    //mem_err = cudaMemcpyToSymbol(constmem_Nsize, &cNsize, sizeof(int) * 5);
    HANDLE_ERROR(cudaMemcpyToSymbol(constmem_Nsize, &cNsize, sizeof(int) * 5));
}
// void set_constant_mem_para(float * para, int N_para)
// {
//     HANDLE_ERROR(cudaMemcpyToSymbol(constmem_para_f, &para, sizeof(float) * N_para));
// }


__inline__ __device__ float gauss2d(float x, float y, float A, float mux, float muy, float sigma)
{// isotropic 2d gaussian function
    if ((sigma < 1e-7f))
    {
        return 0.0f;
    }
    else
    {
        float half_1_sigma2 = 1.0f / (2.0f * sigma * sigma);
        float c = const_1_PI_f * half_1_sigma2;
        float xnew = x - mux;
        float ynew = y - muy;
        return A * c * expf(-half_1_sigma2 * (xnew * xnew + ynew * ynew));
    }
}
__inline__ __device__ double gauss2d(double x, double y, double A, double mux, double muy, double sigma)
{// isotropic 2d gaussian function
    if ((sigma < 1e-7f))
    {
        return 0.0f;
    }
    else
    {
        float half_1_sigma2 = 1.0 / (2.0 * sigma * sigma);
        float c = const_1_PI * half_1_sigma2;
        float xnew = x - mux;
        float ynew = y - muy;
        return A * c * exp(-half_1_sigma2 * (xnew * xnew + ynew * ynew));
    }
}
__inline__ __device__ float mvn2d(float x, float y, float A, float mux, float muy, float sigma1, float sigma2, float beta)
{
    if ((sigma1 < 1e-7f) | (sigma2 < 1e-7f) | (A < 1e-7f))
    {
        // float number effective decimal digits is 7
        return 0.0f;
    }
    else
    {
        // Rotation matrix R, coordinate X = [x,y]-[mux,muy], Y = RX
        float sinb = sinf(beta);
        float cosb = cosf(beta);
        float Y1 = cosb * (x - mux) - sinb * (y - muy);// beta is angle in rad
        float Y2 = sinb * (x - mux) + cosb * (y - muy);
        // S = Y' * M^-1 * Y
        float S1 = (Y1 / sigma1) * (Y1 / sigma1);
        float S2 = (Y2 / sigma2) * (Y2 / sigma2);
        float exponant = -0.5f * (S1 + S2);
        float scale = (A * const_1_2PI_f) / (sigma1 * sigma2);

        return scale * expf(exponant);
    }
}
__inline__ __device__ double mvn2d(double x, double y, double A, double mux, double muy, double sigma1, double sigma2, double beta)
{
    if ((sigma1 < 1e-7) | (sigma2 < 1e-7) | (A < 1e-7))
    {
        // float number effective decimal digits is 7
        return 0.0;
    }
    else
    {
        // Rotation matrix R, coordinate X = [x,y]-[mux,muy], Y = RX
        double sinb = sin(beta);
        double cosb = cos(beta);
        double Y1 = cosb * (x - mux) - sinb * (y - muy);// beta is angle in rad
        double Y2 = sinb * (x - mux) + cosb * (y - muy);
        // S = Y' * M^-1 * Y
        double S1 = (Y1 / sigma1) * (Y1 / sigma1);
        double S2 = (Y2 / sigma2) * (Y2 / sigma2);
        double exponant = -0.5 * (S1 + S2);
        double scale = (A * const_1_2PI) / (sigma1 * sigma2);

        return scale * exp(exponant);
    }
}
template<class T>
__global__ void dose3d_N_iso(T* X, T* Y, T* para, T* dose3d)
{// isotropic 2d gaussian function
    int nxy = blockIdx.x * blockDim.x + threadIdx.x;
    int ng = blockIdx.y * blockDim.y + threadIdx.y;
    int nz = blockIdx.z * blockDim.z + threadIdx.z;
    int Nx = constmem_Nsize[0];
    int Ny = constmem_Nsize[1];
    int Nz = constmem_Nsize[2];
    int N_gaussian = constmem_Nsize[3];
    int nx = nxy / Ny;
    int ny = nxy % Ny;

    T x, y, A1, mux1, muy1, sigma1;
    int64_t idx3d;
    if (nxy < Nx * Ny & ng < N_gaussian & nz < Nz)
    {
        idx3d = nxy + nz * (Nx * Ny);
        x = X[nx];
        y = Y[ny];
        A1 = para[nz * N_gaussian * 4 + 4 * ng];
        mux1 = para[nz * N_gaussian * 4 + 4 * ng + 1];
        muy1 = para[nz * N_gaussian * 4 + 4 * ng + 2];
        sigma1 = para[nz * N_gaussian * 4 + 4 * ng + 3];

        T dist2d = (x - mux1) * (x - mux1) + (y - muy1) * (y - muy1);
        if (dist2d < 16 * sigma1 * sigma1)
        {
            atomicAdd(dose3d + idx3d, gauss2d(x, y, A1, mux1, muy1, sigma1));
        }
    }

}
// template<class T>
// __global__ void dose3d_N_iso(T* X, T* Y, T* para, T* dose3d)
// {// isotropic 2d gaussian function
//     int nx = blockIdx.x * blockDim.x + threadIdx.x;
//     int ny = blockIdx.y * blockDim.y + threadIdx.y;
//     int nz = blockIdx.z * blockDim.z + threadIdx.z;
//     int Nx = constmem_Nsize[0];
//     int Ny = constmem_Nsize[1];
//     int Nz = constmem_Nsize[2];
//     int N_gaussian = constmem_Nsize[3];
//     if (nz < Nz & ny < Ny & nx < Nx)
//     {
//         int64_t idx3d{ nx + ny * Nx + nz * (Nx * Ny) };
//         T x{ X[nx] };
//         T y{ Y[ny] };
//         T temp{};
//         for (int ng = 0; ng < N_gaussian; ++ng)
//         {
//             T A1 = para[nz * N_gaussian * 4 + 4 * ng];
//             T mux1 = para[nz * N_gaussian * 4 + 4 * ng + 1];
//             T muy1 = para[nz * N_gaussian * 4 + 4 * ng + 2];
//             T sigma1 = para[nz * N_gaussian * 4 + 4 * ng + 3];
//             T dist2d = (x-mux1)*(x-mux1) + (y-muy1)*(y-muy1);
//             if (dist2d < 16*sigma1*sigma1)
//             {
//                 temp += gauss2d(x, y, A1, mux1, muy1, sigma1);
//             }
//         }
//         dose3d[idx3d] = temp;
//     }
// }
template<class T>
__global__ void dose3d_N_iso_Gradient(T* X, T* Y, T* para, T* output)
{// isotropic 2d gaussian function
    int nx = blockIdx.x * blockDim.x + threadIdx.x;
    int ny = blockIdx.y * blockDim.y + threadIdx.y;
    int nz = blockIdx.z * blockDim.z + threadIdx.z;
    int Nx = constmem_Nsize[0];
    int Ny = constmem_Nsize[1];
    int Nz = constmem_Nsize[2];
    int N_gaussian = constmem_Nsize[3];
    if (nz < Nz & ny < Ny & nx < Nx)
    {
        int64_t N_gauss_para = int64_t(Nz) * int64_t(N_gaussian) * int64_t(4) ;
        int64_t idx3d{ nx + ny * Nx + nz * (Nx * Ny) };
        T x{ X[nx] };
        T y{ Y[ny] };
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
                T G = gauss2d(x, y, A1, mux1, muy1, sigma1); // avoid G/A when A is small
                
                T One_sigma_sqr = 1.0 / (sigma_sqr);
                T w = (x_mux * x_mux + y_muy * y_muy - 2.0 * sigma_sqr) / (sigma1 * sigma_sqr);
                output[idx3d * N_gauss_para + nz * N_gaussian * 4 + 4 * ng] = G / (A1);                      // dG/dA
                output[idx3d * N_gauss_para + nz * N_gaussian * 4 + 4 * ng + 1] = x_mux * One_sigma_sqr * G; // dG/dmux
                output[idx3d * N_gauss_para + nz * N_gaussian * 4 + 4 * ng + 2] = y_muy * One_sigma_sqr * G; // dG/dmuy
                output[idx3d * N_gauss_para + nz * N_gaussian * 4 + 4 * ng + 3] = w * G;                     // dG/ds
            }
        }
    }
}

template<class T>
__global__ void dose3d_N(T* X, T* Y, T* para, T* dose3d)
{// general 2d gaussian function
    int nx = blockIdx.x * blockDim.x + threadIdx.x;
    int ny = blockIdx.y * blockDim.y + threadIdx.y;
    int nz = blockIdx.z * blockDim.z + threadIdx.z;
    int Nx = constmem_Nsize[0];
    int Ny = constmem_Nsize[1];
    int Nz = constmem_Nsize[2];
    int N_gaussian = constmem_Nsize[3];
    bool inBoundary = (nz < Nz & ny < Ny & nx < Nx);
    if (inBoundary)
    {
        int64_t idx3d{ nx + ny * Nx + nz * (Nx * Ny) };
        T x{ X[nx] };
        T y{ Y[ny] };
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

__global__ void dose3d_N_Gradient(float* X, float* Y, float* para, float* output)
{// general 2d gaussian function
    int nx = blockIdx.x * blockDim.x + threadIdx.x;
    int ny = blockIdx.y * blockDim.y + threadIdx.y;
    int nz = blockIdx.z * blockDim.z + threadIdx.z;
    int Nx = constmem_Nsize[0];
    int Ny = constmem_Nsize[1];
    int Nz = constmem_Nsize[2];
    int N_gaussian = constmem_Nsize[3];
    bool inBoundary = (nz < Nz & ny < Ny & nx < Nx);
    if (inBoundary)
    {
        int64_t N_gauss_para = int64_t(Nz) * int64_t(N_gaussian) * int64_t(6) ;
        int64_t idx3d{ nx + ny * Nx + nz * (Nx * Ny) };
        float x{ X[nx] };
        float y{ Y[ny] };
        for (int ng = 0; ng < N_gaussian; ++ng)
        {
            float A = para[nz * N_gaussian * 6 + 6 * ng];
            float mux = para[nz * N_gaussian * 6 + 6 * ng + 1];
            float muy = para[nz * N_gaussian * 6 + 6 * ng + 2];
            float sigma1 = para[nz * N_gaussian * 6 + 6 * ng + 3];
            float sigma2 = para[nz * N_gaussian * 6 + 6 * ng + 4];
            float beta = para[nz * N_gaussian * 6 + 6 * ng + 5];
            float G = mvn2d(x, y, A, mux, muy, sigma1, sigma2, beta);
            float sigma1_sqr = sigma1 * sigma1;
            float sigma2_sqr = sigma2 * sigma2;
            float One_sigma1_sqr = 1.0f / (sigma1_sqr);
            float One_sigma2_sqr = 1.0f / (sigma2_sqr);
            float x_mux = x - mux;
            float y_muy = y - muy;
            float cosb = cosf(beta);
            float sinb = sinf(beta);
            float Y1 = x_mux * cosb - y_muy * sinb;
            float Y2 = x_mux * sinb + y_muy * cosb;
            float S1 = (Y1 * Y1) / (sigma1_sqr);
            float S2 = (Y2 * Y2) / (sigma2_sqr);
            float Y1_sigma1_sqr = Y1 * One_sigma1_sqr;
            float Y2_sigma2_sqr = Y2 * One_sigma2_sqr;
            float w1 = (Y1_sigma1_sqr * cosb + Y2_sigma2_sqr * sinb);
            float w2 = (-Y1_sigma1_sqr * sinb + Y2_sigma2_sqr * cosb);
            float w3 = (S1 - 1.0f) / (sigma1);
            float w4 = (S2 - 1.0f) / (sigma2);
            float w5 = Y1 * Y2 * (One_sigma1_sqr - One_sigma2_sqr);
            output[idx3d * N_gauss_para + nz * N_gaussian * 6 + 6 * ng] = G / (A);    // dG/dA
            output[idx3d * N_gauss_para + nz * N_gaussian * 6 + 6 * ng + 1] = w1 * G; // dG/dmux
            output[idx3d * N_gauss_para + nz * N_gaussian * 6 + 6 * ng + 2] = w2 * G; // dG/dmuy
            output[idx3d * N_gauss_para + nz * N_gaussian * 6 + 6 * ng + 3] = w3 * G; // dG/ds1
            output[idx3d * N_gauss_para + nz * N_gaussian * 6 + 6 * ng + 4] = w4 * G; // dG/ds2
            output[idx3d * N_gauss_para + nz * N_gaussian * 6 + 6 * ng + 5] = w5 * G; // dG/db
        }
    }

    
    

}
__global__ void dose3d_N_Gradient(double* X, double* Y, double* para, double* output)
{// general 2d gaussian function
    int nx = blockIdx.x * blockDim.x + threadIdx.x;
    int ny = blockIdx.y * blockDim.y + threadIdx.y;
    int nz = blockIdx.z * blockDim.z + threadIdx.z;
    int Nx = constmem_Nsize[0];
    int Ny = constmem_Nsize[1];
    int Nz = constmem_Nsize[2];
    int N_gaussian = constmem_Nsize[3];
    bool inBoundary = (nz < Nz & ny < Ny & nx < Nx);
    if (inBoundary)
    {
        int64_t N_gauss_para = int64_t(Nz) * int64_t(N_gaussian) * int64_t(6) ;
        int64_t idx3d{ nx + ny * Nx + nz * (Nx * Ny) };
        double x{ X[nx] };
        double y{ Y[ny] };
        for (int ng = 0; ng < N_gaussian; ++ng)
        {
            double A = para[nz * N_gaussian * 6 + 6 * ng];
            double mux = para[nz * N_gaussian * 6 + 6 * ng + 1];
            double muy = para[nz * N_gaussian * 6 + 6 * ng + 2];
            double sigma1 = para[nz * N_gaussian * 6 + 6 * ng + 3];
            double sigma2 = para[nz * N_gaussian * 6 + 6 * ng + 4];
            double beta = para[nz * N_gaussian * 6 + 6 * ng + 5];
            double G = mvn2d(x, y, A, mux, muy, sigma1, sigma2, beta);
            double sigma1_sqr = sigma1 * sigma1;
            double sigma2_sqr = sigma2 * sigma2;
            double One_sigma1_sqr = 1.0 / (sigma1_sqr);
            double One_sigma2_sqr = 1.0 / (sigma2_sqr);
            double x_mux = x - mux;
            double y_muy = y - muy;
            double cosb = cos(beta);
            double sinb = sin(beta);
            double Y1 = x_mux * cosb - y_muy * sinb;
            double Y2 = x_mux * sinb + y_muy * cosb;
            double S1 = (Y1 * Y1) / (sigma1_sqr);
            double S2 = (Y2 * Y2) / (sigma2_sqr);
            double Y1_sigma1_sqr = Y1 * One_sigma1_sqr;
            double Y2_sigma2_sqr = Y2 * One_sigma2_sqr;
            double w1 = (Y1_sigma1_sqr * cosb + Y2_sigma2_sqr * sinb);
            double w2 = (-Y1_sigma1_sqr * sinb + Y2_sigma2_sqr * cosb);
            double w3 = (S1 - 1.0) / (sigma1);
            double w4 = (S2 - 1.0) / (sigma2);
            double w5 = Y1 * Y2 * (One_sigma1_sqr - One_sigma2_sqr);
            output[idx3d * N_gauss_para + nz * N_gaussian * 6 + 6 * ng] = G / (A);    // dG/dA
            output[idx3d * N_gauss_para + nz * N_gaussian * 6 + 6 * ng + 1] = w1 * G; // dG/dmux
            output[idx3d * N_gauss_para + nz * N_gaussian * 6 + 6 * ng + 2] = w2 * G; // dG/dmuy
            output[idx3d * N_gauss_para + nz * N_gaussian * 6 + 6 * ng + 3] = w3 * G; // dG/ds1
            output[idx3d * N_gauss_para + nz * N_gaussian * 6 + 6 * ng + 4] = w4 * G; // dG/ds2
            output[idx3d * N_gauss_para + nz * N_gaussian * 6 + 6 * ng + 5] = w5 * G; // dG/db
        }
    }
}
template<class T>
void Gauss2d::cuda_interface(std::vector<T> X, std::vector<T> Y, std::vector<T> para, T* dose3D, int Nx, int Ny, int Nz, int N_para, int N_gaussian)
{
    int64_t dose_size = int64_t(Nx) * int64_t(Ny) * int64_t(Nz);
    set_constant_mem_Nsize(Nx, Ny, Nz, N_gaussian, N_para);
    // copy data to device
    thrust::device_vector<T> X_dev = X;
    thrust::device_vector<T> Y_dev = Y;
    thrust::device_vector<T> para_dev = para;
	thrust::device_vector<T> dose3D_dev(dose_size);
    // cast to raw pointer
    T* X_dev_ptr = thrust::raw_pointer_cast(X_dev.data());
    T* Y_dev_ptr = thrust::raw_pointer_cast(Y_dev.data());
    T* para_dev_ptr = thrust::raw_pointer_cast(para_dev.data());
    T* dose3D_dev_ptr = thrust::raw_pointer_cast(dose3D_dev.data());

    
    if (Nz * N_gaussian * 6 == N_para)
    {
        dim3 threadsPerBlock(8, 8, 8);
        dim3 numBlocks((Nx - 1 + threadsPerBlock.x) / threadsPerBlock.x, (Ny - 1 + threadsPerBlock.y) / threadsPerBlock.y,
            (Nz - 1 + threadsPerBlock.z) / threadsPerBlock.z);
        dose3d_N << <numBlocks, threadsPerBlock >> > (X_dev_ptr, Y_dev_ptr, para_dev_ptr, dose3D_dev_ptr);
    }
    else if (Nz * N_gaussian * 4 == N_para)
    {
        dim3 threadsPerBlock(4, 64, 2);
        dim3 numBlocks((Nx*Ny - 1 + threadsPerBlock.x) / threadsPerBlock.x, (N_gaussian - 1 + threadsPerBlock.y) / threadsPerBlock.y,
            (Nz - 1 + threadsPerBlock.z) / threadsPerBlock.z);
        dose3d_N_iso << <numBlocks, threadsPerBlock >> > (X_dev_ptr, Y_dev_ptr, para_dev_ptr, dose3D_dev_ptr);
    }
    // copy data back to host
    //thrust::copy(dose3D_dev.begin(), dose3D_dev.end(), dose3D.begin());// dose3D is std::vector or host vector
    // direct copy data to host pointer
    HANDLE_ERROR(cudaMemcpy(dose3D, dose3D_dev_ptr, sizeof(T)* dose_size, cudaMemcpyDeviceToHost));
}

void Gauss2d::cuda_interface_gradient(std::vector<float> X, std::vector<float> Y, std::vector<float> para, float* grad, int Nx, int Ny, int Nz, int N_para, int N_gaussian)
{
    int64_t grad_size = int64_t(N_para) * int64_t(Nx) * int64_t(Ny) * int64_t(Nz);
    set_constant_mem_Nsize(Nx, Ny, Nz, N_gaussian, N_para);
    // copy data to device
    thrust::device_vector<float> X_dev = X;
    thrust::device_vector<float> Y_dev = Y;
    thrust::device_vector<float> para_dev = para;
	thrust::device_vector<float> grad_dev(grad_size);
    // cast to raw pointer
    float* X_dev_ptr = thrust::raw_pointer_cast(X_dev.data());
    float* Y_dev_ptr = thrust::raw_pointer_cast(Y_dev.data());
    float* para_dev_ptr = thrust::raw_pointer_cast(para_dev.data());
    float* grad_dev_ptr = thrust::raw_pointer_cast(grad_dev.data());

    dim3 threadsPerBlock(8, 8, 8);
    dim3 numBlocks((Nx - 1 + threadsPerBlock.x) / threadsPerBlock.x, (Ny - 1 + threadsPerBlock.y) / threadsPerBlock.y,
        (Nz - 1 + threadsPerBlock.z) / threadsPerBlock.z);
    if (Nz * N_gaussian * 6 == N_para)
    {
        dose3d_N_Gradient << <numBlocks, threadsPerBlock >> > (X_dev_ptr, Y_dev_ptr, para_dev_ptr, grad_dev_ptr);
    }
    else if (Nz * N_gaussian * 4 == N_para)
    {
        dose3d_N_iso_Gradient << <numBlocks, threadsPerBlock >> > (X_dev_ptr, Y_dev_ptr, para_dev_ptr, grad_dev_ptr);
    }
    // copy data back to host
    HANDLE_ERROR(cudaMemcpy(grad, grad_dev_ptr, sizeof(float)* grad_size, cudaMemcpyDeviceToHost));
}
void Gauss2d::cuda_interface_gradient(std::vector<double> X, std::vector<double> Y, std::vector<double> para, double* grad, int Nx, int Ny, int Nz, int N_para, int N_gaussian)
{
    int64_t grad_size = int64_t(N_para) * int64_t(Nx) * int64_t(Ny) * int64_t(Nz);
    set_constant_mem_Nsize(Nx, Ny, Nz, N_gaussian, N_para);
    // copy data to device
    thrust::device_vector<double> X_dev = X;
    thrust::device_vector<double> Y_dev = Y;
    thrust::device_vector<double> para_dev = para;
	thrust::device_vector<double> grad_dev(grad_size);
    // cast to raw pointer
    double* X_dev_ptr = thrust::raw_pointer_cast(X_dev.data());
    double* Y_dev_ptr = thrust::raw_pointer_cast(Y_dev.data());
    double* para_dev_ptr = thrust::raw_pointer_cast(para_dev.data());
    double* grad_dev_ptr = thrust::raw_pointer_cast(grad_dev.data());

    dim3 threadsPerBlock(8, 8, 8);
    dim3 numBlocks((Nx - 1 + threadsPerBlock.x) / threadsPerBlock.x, (Ny - 1 + threadsPerBlock.y) / threadsPerBlock.y,
        (Nz - 1 + threadsPerBlock.z) / threadsPerBlock.z);
    if (Nz * N_gaussian * 6 == N_para)
    {
        dose3d_N_Gradient << <numBlocks, threadsPerBlock >> > (X_dev_ptr, Y_dev_ptr, para_dev_ptr, grad_dev_ptr);
    }
    else if (Nz * N_gaussian * 4 == N_para)
    {
        dose3d_N_iso_Gradient << <numBlocks, threadsPerBlock >> > (X_dev_ptr, Y_dev_ptr, para_dev_ptr, grad_dev_ptr);
    }
    // copy data back to host
    HANDLE_ERROR(cudaMemcpy(grad, grad_dev_ptr, sizeof(double)* grad_size, cudaMemcpyDeviceToHost));
}

template void Gauss2d::cuda_interface(std::vector<double> X, std::vector<double> Y, std::vector<double> para, double* dose3D, int Nx, int Ny, int Nz, int N_para, int N_gaussian);
template void Gauss2d::cuda_interface(std::vector<float> X, std::vector<float> Y, std::vector<float> para, float* dose3D, int Nx, int Ny, int Nz, int N_para, int N_gaussian);
// template void Gauss2d::cuda_interface_gradient(std::vector<double> X, std::vector<double> Y, std::vector<double> para, double* grad, int Nx, int Ny, int Nz, int N_para, int N_gaussian);
// template void Gauss2d::cuda_interface_gradient(std::vector<float> X, std::vector<float> Y, std::vector<float> para, float* grad, int Nx, int Ny, int Nz, int N_para, int N_gaussian);


// disable matlab entry function
// void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
//     float *X;
//     float *Y;
//     float *para;
//     X = (float*)mxGetPr(prhs[0]);
//     Y = (float*)mxGetPr(prhs[1]);
//     para = (float*)mxGetPr(prhs[2]);

//     const mwSize *dim_X = mxGetDimensions(prhs[0]);
//     const mwSize *dim_Y = mxGetDimensions(prhs[1]);
//     const mwSize *dim_para = mxGetDimensions(prhs[2]);
//     int Nx = static_cast<int>(dim_X[0]*dim_X[1]);
//     int Ny = static_cast<int>(dim_Y[0]*dim_Y[1]);
//     int N_para = static_cast<int>(dim_para[0]*dim_para[1]);

//     int Nz = static_cast<int>(*mxGetPr(prhs[3]));
//     int N_gaussian = static_cast<int>(*mxGetPr(prhs[4]));

//     set_constant_mem_Nsize(Nx, Ny, Nz, N_gaussian, N_para);

//     thrust::host_vector<float> X_vec(Nx), Y_vec(Ny);
//     thrust::copy(X, X + Nx, X_vec.begin());
//     thrust::copy(Y, Y + Ny, Y_vec.begin());

//     thrust::host_vector<float> para_vec(N_para);
//     thrust::copy(para, para + N_para, para_vec.begin());

//     const mwSize size[3]{ mwSize(Nx), mwSize(Ny), mwSize(Nz) };
//     int64_t size3d = int64_t(Nx)* int64_t(Ny)* int64_t(Nz);
//     plhs[0] = mxCreateNumericArray(3, size, mxSINGLE_CLASS, mxREAL);
//     float* dose3d_ptr{};
//     dose3d_ptr = (float*)mxGetPr(plhs[0]);
//     thrust::host_vector<float> dose3d(size3d);


//     cpu_interface(X_vec, Y_vec, para_vec, dose3d, Nx, Ny, Nz, N_para, N_gaussian);
//     thrust::copy(dose3d.begin(), dose3d.end(), dose3d_ptr);

// //     std::string filename = "dose3D";
// //     save_dat(filename,dose3d);
// }


//int main()
//{
//    int Nx = 128;
//    int Ny = 128;
//    int Nz = 1;
//    int N_gaussian = 1;
//    int N_para = 6;
//    set_constant_mem_Nsize(Nx, Ny, Nz, N_gaussian, N_para);
//    thrust::host_vector<float> X(Nx);
//    thrust::host_vector<float> Y(Ny);
//    thrust::host_vector<float> dose3D(Nx * Ny);
//    thrust::host_vector<float> para(6);// A, mux, muy, sigma1, simga2, beta(in rad)
//    para[0] = 1.0f;
//    para[1] = 0.0f;
//    para[2] = 0.0f;
//    para[3] = 1.0f;
//    para[4] = 2.0f;
//    para[5] = 0.0f;
//
//    linspace(X.data(), -1.0f, 1.0f, Nx);
//    linspace(Y.data(), -1.0f, 1.0f, Nx);
//    cpu_interface(X, Y, para, dose3D, Nx, Ny, Nz, N_para, N_gaussian);
//    for (int ix = 0; ix < Nx; ++ix)
//    {
//        for (int iy = 0; iy < Ny; ++iy)
//        {
//            std::cout << dose3D[iy + ix * Ny] << ' ';
//        }
//        std::cout << '\n';
//    }
//    
//    // cudaDeviceReset must be called before exiting in order for profiling and
//    // tracing tools such as Nsight and Visual Profiler to show complete traces.
//    cudaError_t cudaStatus = cudaDeviceReset();
//    if (cudaStatus != cudaSuccess) {
//        fprintf(stderr, "cudaDeviceReset failed!");
//        return 1;
//    }
//    return 0;
//}

