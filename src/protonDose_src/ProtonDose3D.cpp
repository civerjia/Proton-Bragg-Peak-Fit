#include <omp.h>
#include "gauss2d.h"
#include "matrix.h"
#include "mex.h"

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
    mxClassID id0 = mxGetClassID(prhs[0]);
    mxClassID id1 = mxGetClassID(prhs[1]);
    mxClassID id2 = mxGetClassID(prhs[2]);
    mxClassID id3 = mxGetClassID(prhs[2]);
    mxClassID id4 = mxGetClassID(prhs[2]);
    // common parts
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
    if(!((N_gauss_para == Nz*N_gaussian*6) | (N_gauss_para == Nz*N_gaussian*4)))
    {
        mexPrintf("N_gaussian(%d),Nz(%d) not match N_gauss_para(%d)!\n",N_gaussian,Nz,N_gauss_para);
        mexErrMsgTxt("Size not match!\n");
    }
    mxClassID mex_type_class = mxSINGLE_CLASS;
    if(id0 + id1 + id2 + id3 + id4 == 5*mxDOUBLE_CLASS)
    {
        using mex_type1 = double;
        mex_type_class = mxDOUBLE_CLASS; 
        mex_type1 *X;
        mex_type1 *Y;
        mex_type1 *Z;
        mex_type1 *gauss_para;
        mex_type1 *bf_para;
        X = (mex_type1*)mxGetPr(prhs[0]);
        Y = (mex_type1*)mxGetPr(prhs[1]);
        Z = (mex_type1*)mxGetPr(prhs[2]);
        gauss_para = (mex_type1*)mxGetPr(prhs[3]);
        bf_para = (mex_type1*)mxGetPr(prhs[4]);

        // gauss2d 
        std::vector<mex_type1> X_vec(Nx), Y_vec(Ny);
        std::vector<mex_type1> gauss_para_vec(N_gauss_para);
        std::copy(X, X + Nx, X_vec.begin());
        std::copy(Y, Y + Ny, Y_vec.begin());
        std::copy(gauss_para, gauss_para + N_gauss_para, gauss_para_vec.begin());
        // output 
        const mwSize size[3]{ mwSize(Nx), mwSize(Ny), mwSize(Nz) };
        plhs[0] = mxCreateNumericArray(3, size, mex_type_class, mxREAL);
        mex_type1* dose3d_ptr{};
        dose3d_ptr = (mex_type1*)mxGetPr(plhs[0]);
        // run calculation
        std::vector<mex_type1> idd_o(Nz);
        // step 1 get IDD
	    BP::IDD_array_N( Z, idd_o.data(), bf_para, Nz, N_bf_para);
        int Ng = N_gauss_para / (Nz * N_gaussian);
        // step 2 update weight Ai in Gauss2d model
	    updata_Ai(gauss_para_vec.data(), idd_o.data(), Nz, N_gaussian, Ng);
        // step 3 get 3D dose layer by layer
        if(isGPU == 0)
        {
            Gauss2d::interface( X_vec, Y_vec, gauss_para_vec, dose3d_ptr, Nx, Ny, Nz, N_gauss_para, N_gaussian);// cpu version
        }
	    else
        {
            Gauss2d::cuda_interface( X_vec,  Y_vec,  gauss_para_vec, dose3d_ptr, Nx, Ny, Nz, N_gauss_para, N_gaussian);// gpu version
        }
    }
    else if(id0 + id1 + id2 + id3 + id4 == 5*mxSINGLE_CLASS)
    {
        using mex_type2 = float;
        mex_type_class = mxSINGLE_CLASS; 
        mex_type2 *X;
        mex_type2 *Y;
        mex_type2 *Z;
        mex_type2 *gauss_para;
        mex_type2 *bf_para;
        X = (mex_type2*)mxGetPr(prhs[0]);
        Y = (mex_type2*)mxGetPr(prhs[1]);
        Z = (mex_type2*)mxGetPr(prhs[2]);
        gauss_para = (mex_type2*)mxGetPr(prhs[3]);
        bf_para = (mex_type2*)mxGetPr(prhs[4]);

        // gauss2d 
        std::vector<mex_type2> X_vec(Nx), Y_vec(Ny);
        std::vector<mex_type2> gauss_para_vec(N_gauss_para);
        std::copy(X, X + Nx, X_vec.begin());
        std::copy(Y, Y + Ny, Y_vec.begin());
        std::copy(gauss_para, gauss_para + N_gauss_para, gauss_para_vec.begin());
        // output 
        const mwSize size[3]{ mwSize(Nx), mwSize(Ny), mwSize(Nz) };
        plhs[0] = mxCreateNumericArray(3, size, mex_type_class, mxREAL);
        mex_type2* dose3d_ptr{};
        dose3d_ptr = (mex_type2*)mxGetPr(plhs[0]);
        // run calculation
        std::vector<mex_type2> idd_o(Nz);
        // step 1 get IDD
	    BP::IDD_array_N( Z, idd_o.data(), bf_para, Nz, N_bf_para);
        int Ng = N_gauss_para / (Nz * N_gaussian);
        // step 2 update weight Ai in Gauss2d model
	    updata_Ai(gauss_para_vec.data(), idd_o.data(), Nz, N_gaussian, Ng);
        // step 3 get 3D dose layer by layer
        if(isGPU == 0)
        {
            Gauss2d::interface( X_vec, Y_vec, gauss_para_vec, dose3d_ptr, Nx, Ny, Nz, N_gauss_para, N_gaussian);// cpu version
        }
	    else
        {
            Gauss2d::cuda_interface( X_vec,  Y_vec,  gauss_para_vec, dose3d_ptr, Nx, Ny, Nz, N_gauss_para, N_gaussian);// gpu version
        }
    }
    else{
        mexErrMsgTxt("Inconsistency float point input data type!, Should be all double or single\n");
    }
}