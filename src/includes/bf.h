#pragma once
#define _USE_MATH_DEFINES
#include <omp.h>
#include <cmath>
#include <vector>

namespace BP_fast
{
	template <class T>
	void IDD_array_N(T* depth, T* dose_o, T* para_i, int size, int para_size);
	template <class T>
	void get_jacobian(T* depth, T* grad_o, T* para_i, int size, int para_size);
}