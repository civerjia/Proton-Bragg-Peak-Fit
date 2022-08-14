#include "gauss2d.h"
#include "mex.h"

// disable matlab entry function
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    float *X;
    float *Y;
    float *para;
    X = (float*)mxGetPr(prhs[0]);
    Y = (float*)mxGetPr(prhs[1]);
    para = (float*)mxGetPr(prhs[2]);

    const mwSize *dim_X = mxGetDimensions(prhs[0]);
    const mwSize *dim_Y = mxGetDimensions(prhs[1]);
    const mwSize *dim_para = mxGetDimensions(prhs[2]);
    int Nx = static_cast<int>(dim_X[0]*dim_X[1]);
    int Ny = static_cast<int>(dim_Y[0]*dim_Y[1]);
    int N_para = static_cast<int>(dim_para[0]*dim_para[1]);

    int Nz = static_cast<int>(*mxGetPr(prhs[3]));
    int N_gaussian = static_cast<int>(*mxGetPr(prhs[4]));
    int isGPU = static_cast<int>(*mxGetPr(prhs[5]));

    std::vector<float> X_vec(Nx), Y_vec(Ny);
    std::copy(X, X + Nx, X_vec.begin());
    std::copy(Y, Y + Ny, Y_vec.begin());

    const mwSize size[3]{ mwSize(Nx), mwSize(Ny), mwSize(Nz) };
    plhs[0] = mxCreateNumericArray(3, size, mxSINGLE_CLASS, mxREAL);//mxDOUBLE_CLASS
    float* dose3d_ptr{};
    dose3d_ptr = (float*)mxGetPr(plhs[0]);

    std::vector<float> para_vec(N_para);
    std::copy(para, para + N_para, para_vec.begin());

    
    if((N_para == Nz*N_gaussian*6) | (N_para == Nz*N_gaussian*4))
    {
        //get 3D dose layer by layer
        if(isGPU == 0)
        {
            Gauss2d::interface( X_vec, Y_vec, para_vec, dose3d_ptr, Nx, Ny, Nz, N_para, N_gaussian);// cpu version
        }
	    else
        {
            Gauss2d::cuda_interface( X_vec, Y_vec, para_vec, dose3d_ptr, Nx, Ny, Nz, N_para, N_gaussian);// gpu version
        }
    }
    else
    {
        mexPrintf("N_gaussian(%d),Nz(%d) not match para size(%d)!\n",N_gaussian,Nz,N_para);
    }
}