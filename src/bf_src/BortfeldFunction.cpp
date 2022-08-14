#include <cmath>
#include "mex.h"
#include "parabolic_cylinder_function.h"
#include "bp.h"
#include <vector>


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    float * Z{}, * para{};
    int idx{ 1 };
    Z = (float*)mxGetPr(prhs[0]);
    para = (float*)mxGetPr(prhs[1]);
    idx = int(*mxGetPr(prhs[2]));

    const mwSize *dim_Z = mxGetDimensions(prhs[0]);
    const mwSize *dim_para = mxGetDimensions(prhs[1]);
    int Nz = static_cast<int>(dim_Z[0]*dim_Z[1]);
    int N_para = static_cast<int>(dim_para[0]*dim_para[1]);

    if (N_para % 4 == 0)
    {
        // check Nz 
        if (idx == 0)
        {
            // get dose
            const mwSize size[2]{ dim_Z[0], dim_Z[1] };
            plhs[0] = mxCreateNumericArray(2, size, mxSINGLE_CLASS, mxREAL);
            float* idd_o_ptr;
            idd_o_ptr = (float*)mxGetPr(plhs[0]);
            BP::IDD_array_N(Z, idd_o_ptr, para, Nz, N_para);
        }
        else if(idx == 1)
        {
            // get Jacobian
            const mwSize size[2]{ mwSize(Nz), mwSize(N_para) };
            plhs[0] = mxCreateNumericArray(2, size, mxSINGLE_CLASS, mxREAL);
            float * grad_o_ptr;
            grad_o_ptr = (float*)mxGetPr(plhs[0]);
            BP::get_jacobian(Z, grad_o_ptr, para, Nz, N_para);
        }
    }

}

//void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
//    mexPrintf("hello, world\n");
//}
