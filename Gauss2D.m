% help function of Gauss2D
% author : shuang zhou civerjia@gmail.com
% How to use: output = Gauss2D(x,y,gauss_para,Nz,N_gaussian,isGPU,isGrad);
% x         : 1D double/single array, size (1,n) or (n,1). unit is cm;
% y         : 1D double/single array, size (1,n) or (n,1). unit is cm;
% gauss_para: 1D double/single array, parameter of gauss2d function 
%         size (4*Nz*N_gaussian,1) or (1,4*Nz*N_gaussian) or
%         (6*Nz*N_gaussian,1) or (1,6*Nz*N_gaussian)
%         4:(A,mux,muy,sigma), 6:(A,mux,muy,sigmax,sigmay,beta)
% Nz        : int scalar, number of Z layer
% N_gaussian: int scalar, number of gaussian functions per layer
% isGPU     : int scalar, {0, 1} use GPU
% isGrad    : int scalar, {0, 1} return Jacobian
% output: 2D double array or 3D array
%       : if isGrad = 0, output is Dose, size = (Nx,Ny,Nz)
%       : if isGrad = 1, output is jacobian, size = (N_gauss_para,Nx*Ny*Nz)