#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/sort.h>
#include <thrust/copy.h>

#define _USE_MATH_DEFINES
#include <algorithm>
#include <cmath>
#include <array>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <string>
#include <fstream>      // std::ofstream
// #include <filesystem> // C++17, save_dat function

#include <omp.h>

#include "bp.h"
namespace Gauss2d
{
    
    template<class T>
    void cuda_interface(std::vector<T> X, std::vector<T> Y, std::vector<T> para, T* dose3d_ptr, int Nx, int Ny, int Nz, int N_para, int N_gaussian);
    // template<class T>
    // void cuda_interface_gradient(std::vector<T> X, std::vector<T> Y, std::vector<T> para, T* grad, int Nx, int Ny, int Nz, int N_para, int N_gaussian);
    void cuda_interface_gradient(std::vector<double> X, std::vector<double> Y, std::vector<double> para, double* grad, int Nx, int Ny, int Nz, int N_para, int N_gaussian);
    void cuda_interface_gradient(std::vector<float> X, std::vector<float> Y, std::vector<float> para, float* grad, int Nx, int Ny, int Nz, int N_para, int N_gaussian);
    
    template<class T>
    void interface(std::vector<T> X, std::vector<T> Y, std::vector<T> para, T* dose3d_ptr, int Nx, int Ny, int Nz, int N_para, int N_gaussian);

    template<class T>
    void interface_gradient(std::vector<T> X, std::vector<T> Y, std::vector<T> para, T* grad_ptr, int Nx, int Ny, int Nz, int N_para, int N_gaussian);

    template <class T>
    void linspace(T* L, T dmin, T dmax, int N)
    {
        for (int i = 0; i < N; ++i)
        {
            L[i] = dmin + T(i) * (dmax - dmin) / T(N - 1);
        }
    }
    // save data to binary file
    // template <class T>
    // void save_dat(std::string filename, thrust::host_vector<T> vec)
    // {
    //     std::ofstream outdata; // outdata is like cin
    //     std::string fmt{ ".dat" };
    //     outdata.open(filename + fmt, std::ofstream::binary | std::ofstream::out); // opens the file
    //     outdata.write((char*)vec.data(), vec.size() * sizeof(T));
    //     outdata.close();
    // }
    // template <class T>
    // void save_dat(std::string filename, std::vector<T> vec)
    // {
    //     std::ofstream outdata; // outdata is like cin
    //     std::string fmt{ ".dat" };
    //     outdata.open(filename + fmt, std::ofstream::binary | std::ofstream::out); // opens the file
    //     outdata.write((char*)vec.data(), vec.size() * sizeof(T));
    //     outdata.close();
    // }
}

