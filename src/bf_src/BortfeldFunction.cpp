#include <cmath>
#include "mex.h"
#include "matrix.h"
#include "parabolic_cylinder_function.h"
#include "bp.h"
#include "bf.h"
#include <vector>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mxClassID id0 = mxGetClassID(prhs[0]);
    mxClassID id1 = mxGetClassID(prhs[1]);
    // common parts
    const mwSize *dim_Z = mxGetDimensions(prhs[0]);
    const mwSize *dim_para = mxGetDimensions(prhs[1]);
    int Nz = static_cast<int>(dim_Z[0] * dim_Z[1]);
    int N_para = static_cast<int>(dim_para[0] * dim_para[1]);
    if (nrhs == 3)
    {
        if (N_para % 4 != 0)
        {
            mexPrintf("mod(N_para(%d),4) != 0!\n", N_para);
            mexErrMsgTxt("Size not match!\n");
        }

        mxClassID mex_type_class = mxSINGLE_CLASS;
        if (id0 + id1 == 2 * mxDOUBLE_CLASS)
        {
            // double precision use accurate version
            // mexPrintf("Double input!\n");
            using mex_type1 = double;
            mex_type_class = mxDOUBLE_CLASS; // mxDOUBLE_CLASS=6 mxSINGLE_CLASS=7

            mex_type1 *Z{}, *para{};
            Z = (mex_type1 *)mxGetPr(prhs[0]);
            para = (mex_type1 *)mxGetPr(prhs[1]);
            int idx = int(*mxGetPr(prhs[2]));
            if (idx == 0)
            {
                // get dose
                const mwSize size[2]{dim_Z[0], dim_Z[1]};
                plhs[0] = mxCreateNumericArray(2, size, mex_type_class, mxREAL);
                mex_type1 *idd_o_ptr;
                idd_o_ptr = (mex_type1 *)mxGetPr(plhs[0]);
                BP::IDD_array_N(Z, idd_o_ptr, para, Nz, N_para);
            }
            else if (idx == 1)
            {
                // get Jacobian
                const mwSize size[2]{mwSize(Nz), mwSize(N_para)};
                plhs[0] = mxCreateNumericArray(2, size, mex_type_class, mxREAL);
                mex_type1 *grad_o_ptr;
                grad_o_ptr = (mex_type1 *)mxGetPr(plhs[0]);
                BP::get_jacobian(Z, grad_o_ptr, para, Nz, N_para);
            }
        }
        else if (id0 + id1 == 2 * mxSINGLE_CLASS)
        {
            // mexPrintf("Single input!\n");
            using mex_type2 = float;
            mex_type_class = mxSINGLE_CLASS; // mxDOUBLE_CLASS=6 mxSINGLE_CLASS=7

            mex_type2 *Z{}, *para{};
            Z = (mex_type2 *)mxGetPr(prhs[0]);
            para = (mex_type2 *)mxGetPr(prhs[1]);
            int idx = int(*mxGetPr(prhs[2]));
            if (idx == 0)
            {
                // get dose
                const mwSize size[2]{dim_Z[0], dim_Z[1]};
                plhs[0] = mxCreateNumericArray(2, size, mex_type_class, mxREAL);
                mex_type2 *idd_o_ptr;
                idd_o_ptr = (mex_type2 *)mxGetPr(plhs[0]);
                BP::IDD_array_N(Z, idd_o_ptr, para, Nz, N_para);
            }
            else if (idx == 1)
            {
                // get Jacobian
                const mwSize size[2]{mwSize(Nz), mwSize(N_para)};
                plhs[0] = mxCreateNumericArray(2, size, mex_type_class, mxREAL);
                mex_type2 *grad_o_ptr;
                grad_o_ptr = (mex_type2 *)mxGetPr(plhs[0]);
                BP::get_jacobian(Z, grad_o_ptr, para, Nz, N_para);
            }
        }
        else
        {
            mexErrMsgTxt("Inconsistency float point input data type!, Should be all double or single\n");
        }
    }
    else if(nrhs == 4){
        if (N_para % 4 != 0)
        {
            mexPrintf("mod(N_para(%d),4) != 0!\n", N_para);
            mexErrMsgTxt("Size not match!\n");
        }

        mxClassID mex_type_class = mxSINGLE_CLASS;
        if (id0 + id1 == 2 * mxDOUBLE_CLASS)
        {
            // mexPrintf("Double input!\n");
            using mex_type1 = double;
            mex_type_class = mxDOUBLE_CLASS; // mxDOUBLE_CLASS=6 mxSINGLE_CLASS=7

            mex_type1 *Z{}, *para{};
            Z = (mex_type1 *)mxGetPr(prhs[0]);
            para = (mex_type1 *)mxGetPr(prhs[1]);
            int idx = int(*mxGetPr(prhs[2]));
            if (idx == 0)
            {
                // get dose
                const mwSize size[2]{dim_Z[0], dim_Z[1]};
                plhs[0] = mxCreateNumericArray(2, size, mex_type_class, mxREAL);
                mex_type1 *idd_o_ptr;
                idd_o_ptr = (mex_type1 *)mxGetPr(plhs[0]);
                BP_fast::IDD_array_N(Z, idd_o_ptr, para, Nz, N_para);
            }
            else if (idx == 1)
            {
                // get Jacobian
                const mwSize size[2]{mwSize(Nz), mwSize(N_para)};
                plhs[0] = mxCreateNumericArray(2, size, mex_type_class, mxREAL);
                mex_type1 *grad_o_ptr;
                grad_o_ptr = (mex_type1 *)mxGetPr(plhs[0]);
                BP_fast::get_jacobian(Z, grad_o_ptr, para, Nz, N_para);
            }
        }
        else if (id0 + id1 == 2 * mxSINGLE_CLASS)
        {
            // mexPrintf("Single input!\n");
            using mex_type2 = float;
            mex_type_class = mxSINGLE_CLASS; // mxDOUBLE_CLASS=6 mxSINGLE_CLASS=7

            mex_type2 *Z{}, *para{};
            Z = (mex_type2 *)mxGetPr(prhs[0]);
            para = (mex_type2 *)mxGetPr(prhs[1]);
            int idx = int(*mxGetPr(prhs[2]));
            if (idx == 0)
            {
                // get dose
                const mwSize size[2]{dim_Z[0], dim_Z[1]};
                plhs[0] = mxCreateNumericArray(2, size, mex_type_class, mxREAL);
                mex_type2 *idd_o_ptr;
                idd_o_ptr = (mex_type2 *)mxGetPr(plhs[0]);
                BP_fast::IDD_array_N(Z, idd_o_ptr, para, Nz, N_para);
            }
            else if (idx == 1)
            {
                // get Jacobian
                const mwSize size[2]{mwSize(Nz), mwSize(N_para)};
                plhs[0] = mxCreateNumericArray(2, size, mex_type_class, mxREAL);
                mex_type2 *grad_o_ptr;
                grad_o_ptr = (mex_type2 *)mxGetPr(plhs[0]);
                BP_fast::get_jacobian(Z, grad_o_ptr, para, Nz, N_para);
            }
        }
        else
        {
            mexErrMsgTxt("Inconsistency float point input data type!, Should be all double or single\n");
        }
    }
    else{
        mexErrMsgTxt("Wrong number of input arguments!\n");
    }
}

// void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
//     mexPrintf("hello, world\n");
// }
