addpath('../data/');
addpath('..');
%% speed test BortfeldFunction 
% BortfeldFunction time: single(0.249554s), double(0.248533s)
% BortfeldFunction Grad time: single(2.758517s), double(2.772066s)
t = zeros(4,1);
N = 1e4;
isSingle = 1;
t(1) = test_bf(N,isSingle,0);
isSingle = 0;
t(2) = test_bf(N,isSingle,0);
fprintf("BortfeldFunction time: single(%fs), double(%fs)\n",t(1),t(2));
isSingle = 1;
t(3) = test_bf(N,isSingle,1);
isSingle = 0;
t(4) = test_bf(N,isSingle,1);
fprintf("BortfeldFunction Grad time: single(%fs), double(%fs)\n",t(3),t(4));
disp("No significant difference, because they use the same PCF part(double), PCF cannot work with single.")
%% speed test Gauss2D 
% Gauss2D Grad GPU time: single(0.942952s), double(1.575574s)
% Gauss2D Grad CPU time: single(0.131198s), double(0.201214s)
% Gauss2D GPU time: single(0.196284s), double(0.450232s)
% Gauss2D CPU time: single(0.842755s), double(0.665228s)
N = 1;
t = zeros(2,1);
config = [1,1,1];% isSingle,isGPU,isGrad
t(1) = test_gauss2d(N,config(1),config(2),config(3));
config = [0,1,1];% isSingle,isGPU,isGrad
t(2) = test_gauss2d(N,config(1),config(2),config(3));
fprintf("Gauss2D Grad GPU time: single(%fs), double(%fs)\n",t(1),t(2));

config = [1,0,1];% isSingle,isGPU,isGrad
t(1) = test_gauss2d(N,config(1),config(2),config(3));
config = [0,0,1];% isSingle,isGPU,isGrad
t(2) = test_gauss2d(N,config(1),config(2),config(3));
fprintf("Gauss2D Grad CPU time: single(%fs), double(%fs)\n",t(1),t(2));

N = 1e2;
config = [1,1,0];% isSingle,isGPU,isGrad
t(1) = test_gauss2d(N,config(1),config(2),config(3));
config = [0,1,0];% isSingle,isGPU,isGrad
t(2) = test_gauss2d(N,config(1),config(2),config(3));
fprintf("Gauss2D GPU time: single(%fs), double(%fs)\n",t(1),t(2));

config = [1,0,0];% isSingle,isGPU,isGrad
t(1) = test_gauss2d(N,config(1),config(2),config(3));
config = [0,0,0];% isSingle,isGPU,isGrad
t(2) = test_gauss2d(N,config(1),config(2),config(3));
fprintf("Gauss2D CPU time: single(%fs), double(%fs)\n",t(1),t(2));
%% speed test ProtonDose3D
% ProtonDose3D GPU time: single(0.493263s), double(0.537123s)
% ProtonDose3D CPU time: single(0.751850s), double(0.496489s)
N = 1e2;
t = zeros(2,1);
config = [1,1];% isSingle,isGPU
t(1) = test_protondose3d(N,config(1),config(2));
config = [0,1];% isSingle,isGPU
t(2) = test_protondose3d(N,config(1),config(2));
fprintf("ProtonDose3D GPU time: single(%fs), double(%fs)\n",t(1),t(2));

config = [1,0];% isSingle,isGPU
t(1) = test_protondose3d(N,config(1),config(2));
config = [0,0];% isSingle,isGPU
t(2) = test_protondose3d(N,config(1),config(2));
fprintf("ProtonDose3D CPU time: single(%fs), double(%fs)\n",t(1),t(2));
%% test functions
function [t,idd_o] = test_bf(N,isSingle,isGrad)
    if isSingle
        z = single(linspace(0,19,64));
        bf_para = single([15,0.3,1e-3,0.4, 12,0.4,1e-3,0.4]);
    else
        z = (linspace(0,19,64));
        bf_para = ([15,0.3,1e-3,0.4, 12,0.4,1e-3,0.4]);
    end
    tic;
    for i = 1:N
        idd_o = BortfeldFunction(z,bf_para,isGrad);
    end

    t = toc;
end
function [t,dose] = test_gauss2d(N,isSingle,isGPU,isGrad)
    Nz = 64;
    if isSingle
        x = single(((1:128)-64.5)*0.2);
        y = x;
        gauss_para = single(repmat([0.5,-2,-3,1,3,45*pi/180, 0.5,2,3,3,1,15*pi/280],1,Nz));
    else
        x = ((1:128)-64.5)*0.2;
        y = x;
        gauss_para = repmat([0.5,-2,-3,1,3,45*pi/180, 0.5,2,3,3,1,15*pi/280],1,Nz);
    end
    N_gaussian = 2;
    tic;
    for i = 1:N
        dose = Gauss2D(x,y,gauss_para,Nz,N_gaussian,isGPU,isGrad);
    end
    t =  toc;
end
function [t,dose3d] = test_protondose3d(N,isSingle,isGPU)
    Nz = 64;
    if isSingle
        x = single(((1:128)-64.5)*0.2);
        y = x;
        gauss_para = single(repmat([0.5,-2,-3,1,3,45*pi/180, 0.5,2,3,3,1,15*pi/180],1,Nz));
        z = single(linspace(0,19,Nz));
        bf_para = single([15,0.3,1e-3,0.4, 12,0.4,1e-3,0.4]);
    else
        x = ((1:128)-64.5)*0.2;
        y = x;
        gauss_para = repmat([0.5,-2,-3,1,3,45*pi/180, 0.5,2,3,3,1,15*pi/180],1,Nz);
        z = (linspace(0,19,Nz));
        bf_para = ([15,0.3,1e-3,0.4, 12,0.4,1e-3,0.4]);
    end
    N_gaussian = 2;
    tic;
    for i = 1:N
        dose3d = ProtonDose3D(x,y,z,gauss_para,bf_para,N_gaussian,isGPU);
    end
    t = toc;
end