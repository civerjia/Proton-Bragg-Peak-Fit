#include <cmath>
#include "mex.h"
#include "parabolic_cylinder_function.h"
#include "bp.h"
#include <vector>


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    double * z{}, * para{}, * idd_o{}, * grad_o{};
    int size{ 4 }, idx{ 1 }, para_size{4};
    z = mxGetPr(prhs[0]);
    para = mxGetPr(prhs[1]);
    idx = int(*mxGetPr(prhs[2]));

    if (int(mxGetN(prhs[0])) <= int(mxGetM(prhs[0])))
    {
        size = int(mxGetM(prhs[0]));
    }
    else
    {
        size = int(mxGetN(prhs[0]));
    }
    if (int(mxGetN(prhs[1])) <= int(mxGetM(prhs[1])))
    {
        para_size = int(mxGetM(prhs[1]));
    }
    else
    {
        para_size = int(mxGetN(prhs[1]));
    }

    int num_bps = para_size / 4; // number of bragg peaks
    if (para_size % 4 == 0)
    {
        // check size 
        if (idx == 0)
        {
            // get dose
            plhs[0] = mxCreateDoubleMatrix(size, 1, mxREAL);
            idd_o = mxGetPr(plhs[0]);
            BP::IDD_array_new(z, idd_o, para, size, para_size);
        }
        else if(idx == 1)
        {
            // get Jacobian
            plhs[0] = mxCreateDoubleMatrix(size, para_size, mxREAL);
            grad_o = mxGetPr(plhs[0]);
            BP::get_jacobian(z, grad_o, para, size, para_size);
        }
        else
        {
            // get mean gradient
            plhs[0] = mxCreateDoubleMatrix(para_size, 1, mxREAL);
            grad_o = mxGetPr(plhs[0]);
            BP::get_mean_grad(z, grad_o, para, size, para_size);
        }
    }

}

//void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
//    mexPrintf("hello, world\n");
//}
