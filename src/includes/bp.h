#pragma once
#define _USE_MATH_DEFINES
#include "parabolic_cylinder_function.h"
#include <omp.h>
#include <cmath>
#include <vector>
namespace BP
{
	template <class T>
	T IDD(T depth, T R0, T sigma, T epsilon, T Phi);
	template <class T>
	T IDDV2(T depth, T R0, T sigma, T epsilon, T Phi, T a, T b, T c, T d);
	template <class T>
	T Bragg_hat(T z, T r, T s, T e, T p);
	template <class T>
	T Bragg(T z, T r, T s, T e, T p);
	template <class T>
	void IDD_array(T* depth, T* dose_o, int size, T R0, T sigma, T epsilon, T Phi);
	template <class T>
	void IDD_array_v2(T* depth, T* dose_o, int size, T R0, T sigma, T epsilon, T Phi, T a, T b, T c, T d);
	template <class T>
	void grad(T z, T* grad4_o, T r, T s, T e, T p);
	template <class T>
	void grad_array(T* z, T* grad_o, int size, T* para);
	template <class T>
	void fitBragg(T* z, T* dose, int size, T* para, T* lambda, int niter, T* idd_o, T* loss_o, T* para_o, T& Lmin);
	template <class T>
	void IDD_array_N(T* depth, T* dose_o, T* para_i, int size, int para_size);
	template <class T>
	void IDD_array_N_v2(T* depth, T* dose_o, T* para_i, int size, int para_size);
	template <class T>
	void get_mean_grad(T* depth, T* grad_o, T* para_i, int size, int para_size);
	template <class T>
	void get_jacobian(T* depth, T* grad_o, T* para_i, int size, int para_size);
}

