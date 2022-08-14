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

using T = float;
using host_vec = thrust::host_vector<float>;
using device_vec = thrust::device_vector<float>;

__constant__ const T const_1_PI = 0.318309886183791;
__constant__ const T const_1_SQRT_2PI = 0.398942280401433;
// 1/sqrt(2*pi), 1D gauss constant
__constant__ const T const_1_2PI = 0.159154943091895;
// 1/(2*pi), 2D gauss constant
__constant__ int constmem_Nsize[5];
// declare constant memory, Nx,Ny,Nz,N_gaussian,N_para

void set_constant_mem(int Nx, int Ny, int Nz, int N_gaussian, int N_para)
{
    int cNsize[5] = { Nx,Ny,Nz,N_gaussian,N_para }; // copy host data to constant memory
    //cudaError_t mem_err;
    //mem_err = cudaMemcpyToSymbol(constmem_Nsize, &cNsize, sizeof(int) * 5);
    cudaMemcpyToSymbol(constmem_Nsize, &cNsize, sizeof(int) * 5);
}

__inline__ __device__ T gauss1d(T x, T A, T mu, T sigma)
{
    T c{ const_1_SQRT_2PI / sigma };
    T xnew = (x - mu) / sigma;
    return A * c * expf(-0.5f * xnew * xnew);
}

__inline__ __device__ T gauss2d(T x, T y, T A, T mux, T muy, T sigma)
{// isotropic 2d gaussian function
    if ((sigma < 1e-7f))
    {
        return 0.0f;
    }
    else
    {
        T half_1_sigma2 = 1.0f / (2.0f * sigma * sigma);
        T c = const_1_PI * half_1_sigma2;
        T xnew = x - mux;
        T ynew = y - muy;
        return A * c * expf(-half_1_sigma2 * (xnew * xnew + ynew * ynew));
    }
}
__inline__ __device__ T mvn2d(T x, T y, T A, T mux, T muy, T sigma1, T sigma2, T beta)
{
    if ((sigma1 < 1e-7f) | (sigma2 < 1e-7f) | (A < 1e-7f))
    {
        // float number effective decimal digits is 7
        return 0.0f;
    }
    else
    {
        // Rotation matrix R, coordinate X = [x,y]-[mux,muy], Y = RX
        T sinb = sinf(beta);
        T cosb = cosf(beta);
        T Y1 = cosb * (x - mux) - sinb * (y - muy);// beta is angle in rad
        T Y2 = sinb * (x - mux) + cosb * (y - muy);
        // S = Y' * M^-1 * Y
        T S1 = (Y1 / sigma1) * (Y1 / sigma1);
        T S2 = (Y2 / sigma2) * (Y2 / sigma2);
        T exponant = -0.5f * (S1 + S2);
        T scale = (A * const_1_2PI) / (sigma1 * sigma2);

        return scale * expf(exponant);
    }
}
__global__ void dose3d_N_iso(T* X, T* Y, T* para, T* dose3d)
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
        int idx3d{ nx + ny * Nx + nz * (Nx * Ny) };
        T x{ X[nx] };
        T y{ Y[ny] };
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


__global__ void dose3d_N(T* X, T* Y, T* para, T* dose3d)
{// general 2d gaussian function
    int nx = blockIdx.x * blockDim.x + threadIdx.x;
    int ny = blockIdx.y * blockDim.y + threadIdx.y;
    int nz = blockIdx.z * blockDim.z + threadIdx.z;
    int Nx = constmem_Nsize[0];
    int Ny = constmem_Nsize[1];
    int Nz = constmem_Nsize[2];
    int N_gaussian = constmem_Nsize[3];
    // __shared__ X_shared[512];
    // __shared__ Y_shared[512];
    if (nx >= Nx || ny >= Ny || nz >= Nz) return;

    // X_shared[nx] = X[nx];
    // Y_shared[ny] = Y[ny];
    // __syncthreads();
    // 
    int idx3d{ nx + ny * Nx + nz * (Nx * Ny) };
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

void Gauss2d::cuda_interface(std::vector<float> X, std::vector<float> Y, std::vector<float> para, float* dose3D, int Nx, int Ny, int Nz, int N_para, int N_gaussian)
{
    int64_t dose_size = int64_t(Nx) * int64_t(Ny) * int64_t(Nz);
    set_constant_mem(Nx, Ny, Nz, N_gaussian, N_para);
    // copy data to device
    device_vec X_dev = X;
    device_vec Y_dev = Y;
    device_vec para_dev = para;
	device_vec dose3D_dev(dose_size);
    // cast to raw pointer
    T* X_dev_ptr = thrust::raw_pointer_cast(X_dev.data());
    T* Y_dev_ptr = thrust::raw_pointer_cast(Y_dev.data());
    T* para_dev_ptr = thrust::raw_pointer_cast(para_dev.data());
    T* dose3D_dev_ptr = thrust::raw_pointer_cast(dose3D_dev.data());

    dim3 threadsPerBlock(8, 8, 8);
    dim3 numBlocks((Nx - 1 + threadsPerBlock.x) / threadsPerBlock.x, (Ny - 1 + threadsPerBlock.y) / threadsPerBlock.y,
        (Nz - 1 + threadsPerBlock.z) / threadsPerBlock.z);
    if (Nz * N_gaussian * 6 == N_para)
    {
        dose3d_N << <numBlocks, threadsPerBlock >> > (X_dev_ptr, Y_dev_ptr, para_dev_ptr, dose3D_dev_ptr);
    }
    else if (Nz * N_gaussian * 4 == N_para)
    {
        dose3d_N_iso << <numBlocks, threadsPerBlock >> > (X_dev_ptr, Y_dev_ptr, para_dev_ptr, dose3D_dev_ptr);
    }
    // copy data back to host
    //thrust::copy(dose3D_dev.begin(), dose3D_dev.end(), dose3D.begin());// dose3D is std::vector or host vector
    // direct copy data to host pointer
    cudaMemcpy(dose3D, dose3D_dev_ptr, sizeof(float)* dose_size, cudaMemcpyDeviceToHost);
}

// disable matlab entry function
// void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
//     T *X;
//     T *Y;
//     T *para;
//     X = (T*)mxGetPr(prhs[0]);
//     Y = (T*)mxGetPr(prhs[1]);
//     para = (T*)mxGetPr(prhs[2]);

//     const mwSize *dim_X = mxGetDimensions(prhs[0]);
//     const mwSize *dim_Y = mxGetDimensions(prhs[1]);
//     const mwSize *dim_para = mxGetDimensions(prhs[2]);
//     int Nx = static_cast<int>(dim_X[0]*dim_X[1]);
//     int Ny = static_cast<int>(dim_Y[0]*dim_Y[1]);
//     int N_para = static_cast<int>(dim_para[0]*dim_para[1]);

//     int Nz = static_cast<int>(*mxGetPr(prhs[3]));
//     int N_gaussian = static_cast<int>(*mxGetPr(prhs[4]));

//     set_constant_mem(Nx, Ny, Nz, N_gaussian, N_para);

//     host_vec X_vec(Nx), Y_vec(Ny);
//     thrust::copy(X, X + Nx, X_vec.begin());
//     thrust::copy(Y, Y + Ny, Y_vec.begin());

//     host_vec para_vec(N_para);
//     thrust::copy(para, para + N_para, para_vec.begin());

//     const mwSize size[3]{ mwSize(Nx), mwSize(Ny), mwSize(Nz) };
//     int64_t size3d = int64_t(Nx)* int64_t(Ny)* int64_t(Nz);
//     plhs[0] = mxCreateNumericArray(3, size, mxSINGLE_CLASS, mxREAL);
//     T* dose3d_ptr{};
//     dose3d_ptr = (T*)mxGetPr(plhs[0]);
//     host_vec dose3d(size3d);


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
//    set_constant_mem(Nx, Ny, Nz, N_gaussian, N_para);
//    host_vec X(Nx);
//    host_vec Y(Ny);
//    host_vec dose3D(Nx * Ny);
//    host_vec para(6);// A, mux, muy, sigma1, simga2, beta(in rad)
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

