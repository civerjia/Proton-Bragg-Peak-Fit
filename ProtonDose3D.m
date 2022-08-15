% help function of ProtonDose3D
% author : shuang zhou civerjia@gmail.com
% How to use: output = ProtonDose3D(x,y,z,gauss_para,bf_para,N_gaussian,isGPU);
% x         : 1D double/single array, size (1,n) or (n,1). unit is cm;
% y         : 1D double/single array, size (1,n) or (n,1). unit is cm;
% z         : 1D double/single array, size (1,n) or (n,1). unit is cm;
% gauss_para: 1D double/single array, parameter of gauss2d function 
%         size (4*Nz*N_gaussian,1) or (1,4*Nz*N_gaussian) or
%         (6*Nz*N_gaussian,1) or (1,6*Nz*N_gaussian)
%         4:(A,mux,muy,sigma), 6:(A,mux,muy,sigmax,sigmay,beta)
% bf_para   : 1D double/single array, parameter of BortfeldFunction 
% N_gaussian: int scalar, number of gaussian functions per layer
% isGPU     : int scalar, {0, 1} use GPU
% output: 2D double array or 3D array
%       : output is 3D Dose, size = (Nx,Ny,Nz)