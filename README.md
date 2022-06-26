# Proton Bragg Peak Fit

- Provide 3 different IDD measurement data, data acquired from Zebra and other MLIC which can be used in Proton Radiographics
- Bortfeld function implemented in C++, provide IDD, mean gradient and jacobian outputs.
- Compiled with Visual Studio + Intel OneAPI, faster than mex in MATLAB

How to use:  
- `output = bf_mex((1:64)*0.3,[15,0.3,1e-3,0.4, 12,0.4,1e-3,0.4],'idd')`
- `[x,idd_o] = precise_fit(z,idd_i,num_bp,strict);`
- `x = fast_fit(z,idd_i,num_bp);`

Details can be found in demo.m, column 1D array is prefered such as, zeros(n,1)

Zebra data:

<img src="./Zebra_fit.png" width="500" height="500">

Multi-Bragg Peaks Fit:

<img src="./IDD_fit.png" width="500" height="500">

It take 160s to fit 11057 IDD curves with 2 bragg peak model @i9-9900k

It take 120s to fit 10498 IDD curves with 2 bragg peak model @i9-9900k

Reference :
- An analytical approximation of depth-dose distributions for therapeutic proton beams
- E. Cojocaru. Parabolic Cylinder Functions (https://www.mathworks.com/matlabcentral/fileexchange/22620-parabolic-cylinder-functions)

