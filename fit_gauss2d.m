x = (((1:128)-64.5)*0.2)';
y = x;
xy = [x,y];
% gauss_para = ([0.5,-2,-3,1,3,45*pi/180, 0.5,2,3,3,1,15*pi/280]);
gauss_para = ([0.5,-2,-3,1, 0.5,2,3,1]);
Nz = 1;
N_gaussian = 2;
isGPU = 1;
isGrad = 0;
dose = gauss2d_optim(gauss_para,xy);

x0 = gauss_para + 0.5*rand(1);
% lb = (repmat([0,-20,-20,1e-7,1e-7,0],1,N_gaussian));
% ub = (repmat([1, 20, 20,  10,  10,pi]',1,N_gaussian));
lb = (repmat([0,-20,-20,1e-7],1,N_gaussian));
ub = (repmat([1, 20, 20,  10]',1,N_gaussian));

options = optimoptions('lsqcurvefit','display','none','SpecifyObjectiveGradient',true);
[x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@gauss2d_optim,x0,xy,dose(:),lb,ub,options);
[dose_x,J] = gauss2d_optim(x,xy);
resnorm
e = J - full(jacobian);
imagesc(reshape(dose,128,128) - reshape(dose_x,128,128));
%%
% gradient check
eps = 1e-6;
for i = 1:4
    delta = zeros(1,8);
    delta(i) = eps;
    [dose_x2,J] = gauss2d_optim(x+delta,xy);
    [dose_x1,J] = gauss2d_optim(x,xy);
    err = J(:,i) - (dose_x2 - dose_x1)/eps;
    max_err = max(abs(err));
    assert(max_err < 1e-7,'Max grad err = %f, %dth parameter failed',max_err,i);
end
function [F,J] = gauss2d_optim(gauss_para,xy)%[F,J]
    % depth : 1D double array, size (1,n) or (n,1). unit is cm;
    % para  : 1D double array, parameter of bortfeld function
    %         size must be (4*m,1) or (1,4*m)  [range, sigma, slope, Phi]
    %         m = {1,2,3...}, m denotes the number of bragg peaks
    isGPU = 1;
    Nz = evalin('base', 'Nz');
    N_gaussian = evalin('base', 'N_gaussian');
    % single(xy(:,1)),single(xy(:,2)),single(gauss_para)
    F = double(Gauss2D(single(xy(:,1)),single(xy(:,2)),single(gauss_para),Nz,N_gaussian,isGPU,0));% objective function values at x
    F = F(:);
    if nargout > 1   % two output arguments
        % return Jacobian size (Nx*Ny*Nz,N_gauss_para) 
        J = double(Gauss2D(single(xy(:,1)),single(xy(:,2)),single(gauss_para),Nz,N_gaussian,isGPU,1))';
    end
end