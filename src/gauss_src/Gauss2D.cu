#include "gauss2d.h"
#include "mex.h"

// disable matlab entry function
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    mxClassID id0 = mxGetClassID(prhs[0]);
    mxClassID id1 = mxGetClassID(prhs[1]);
    mxClassID id2 = mxGetClassID(prhs[2]);
    // common parts
    const mwSize *dim_X = mxGetDimensions(prhs[0]);
    const mwSize *dim_Y = mxGetDimensions(prhs[1]);
    const mwSize *dim_para = mxGetDimensions(prhs[2]);
    int Nx = static_cast<int>(dim_X[0]*dim_X[1]);
    int Ny = static_cast<int>(dim_Y[0]*dim_Y[1]);
    int N_para = static_cast<int>(dim_para[0]*dim_para[1]);
    int Nz = static_cast<int>(*mxGetPr(prhs[3]));
    int N_gaussian = static_cast<int>(*mxGetPr(prhs[4]));
    int isGPU = static_cast<int>(*mxGetPr(prhs[5]));// keep api the same
    int index = static_cast<int>(*mxGetPr(prhs[6]));// indicator

    if(!((N_para == Nz*N_gaussian*6) | (N_para == Nz*N_gaussian*4)))
    {
        mexPrintf("N_gaussian(%d),Nz(%d) not match para size(%d)!\n",N_gaussian,Nz,N_para);
        mexErrMsgTxt("Size not match!\n");
    }
    mxClassID mex_type_class = mxSINGLE_CLASS;
    if(id0 + id1 + id2 == 3*mxDOUBLE_CLASS)
    {
        //mexPrintf("Double input!\n");
        using mex_type1 = double;
        mex_type_class = mxDOUBLE_CLASS; // mxDOUBLE_CLASS=6 mxSINGLE_CLASS

        mex_type1 *X;
        mex_type1 *Y;
        mex_type1 *para;
        X = (mex_type1*)mxGetPr(prhs[0]);
        Y = (mex_type1*)mxGetPr(prhs[1]);
        para = (mex_type1*)mxGetPr(prhs[2]);
        std::vector<mex_type1> X_vec(Nx), Y_vec(Ny);
        std::vector<mex_type1> para_vec(N_para);
        std::copy(X, X + Nx, X_vec.begin());
        std::copy(Y, Y + Ny, Y_vec.begin());
        std::copy(para, para + N_para, para_vec.begin());
        if(index == 0)
        {
            const mwSize size[3]{ mwSize(Nx), mwSize(Ny), mwSize(Nz) };
            plhs[0] = mxCreateNumericArray(3, size, mex_type_class, mxREAL);//mxDOUBLE_CLASS mxSINGLE_CLASS
            mex_type1* dose3d_ptr{};
            dose3d_ptr = (mex_type1*)mxGetPr(plhs[0]);
            //get 3D dose layer by layer
            Gauss2d::interface( X_vec, Y_vec, para_vec, dose3d_ptr, Nx, Ny, Nz, N_para, N_gaussian);// cpu version
        }
        else
        {
            const mwSize size[2]{ mwSize(N_para) ,mwSize(Nx)*mwSize(Ny)*mwSize(Nz) };
            plhs[0] = mxCreateNumericArray(2, size, mex_type_class, mxREAL);//mxDOUBLE_CLASS mxSINGLE_CLASS
            mex_type1* grad_ptr{};
            grad_ptr = (mex_type1*)mxGetPr(plhs[0]);
            //get gradient
            Gauss2d::interface_gradient( X_vec, Y_vec, para_vec, grad_ptr, Nx, Ny, Nz, N_para, N_gaussian);// cpu version
        }
    }
    else if(id0 + id1 + id2 == 3*mxSINGLE_CLASS)
    {
        //mexPrintf("Single input!\n");
        using mex_type2 = float;
        mex_type_class = mxSINGLE_CLASS;

        mex_type2 *X;
        mex_type2 *Y;
        mex_type2 *para;
        X = (mex_type2*)mxGetPr(prhs[0]);
        Y = (mex_type2*)mxGetPr(prhs[1]);
        para = (mex_type2*)mxGetPr(prhs[2]);
        std::vector<mex_type2> X_vec(Nx), Y_vec(Ny);
        std::vector<mex_type2> para_vec(N_para);
        std::copy(X, X + Nx, X_vec.begin());
        std::copy(Y, Y + Ny, Y_vec.begin());
        std::copy(para, para + N_para, para_vec.begin());
        if(index == 0)
        {
            const mwSize size[3]{ mwSize(Nx), mwSize(Ny), mwSize(Nz) };
            plhs[0] = mxCreateNumericArray(3, size, mex_type_class, mxREAL);//mxDOUBLE_CLASS mxSINGLE_CLASS
            mex_type2* dose3d_ptr{};
            dose3d_ptr = (mex_type2*)mxGetPr(plhs[0]);
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
            const mwSize size[2]{ mwSize(N_para) ,mwSize(Nx)*mwSize(Ny)*mwSize(Nz) };
            plhs[0] = mxCreateNumericArray(2, size, mex_type_class, mxREAL);//mxDOUBLE_CLASS mxSINGLE_CLASS
            mex_type2* grad_ptr{};
            grad_ptr = (mex_type2*)mxGetPr(plhs[0]);
            //get gradient
            if(isGPU == 0)
            {
                Gauss2d::interface_gradient( X_vec, Y_vec, para_vec, grad_ptr, Nx, Ny, Nz, N_para, N_gaussian);// cpu version
            }
            else
            {
                Gauss2d::cuda_interface_gradient( X_vec, Y_vec, para_vec, grad_ptr, Nx, Ny, Nz, N_para, N_gaussian);// gpu version
            }
        }
    }
    else{
        mexErrMsgTxt("Inconsistency float point input data type!, Should be all double or single\n");
    }

}