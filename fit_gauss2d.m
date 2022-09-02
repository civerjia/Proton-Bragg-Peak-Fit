x = (((1:128)-64.5)*0.2)';
y = x;
Nz = 1;
N_gaussian = 2;
% calculate dose CPU version
doseFun = @(gauss_para,xy) Gauss2D(xy(:,1),xy(:,2),gauss_para,Nz,N_gaussian,0,0);
% gauss_para = ([0.5,-2,-3,1,3,45*pi/180, 0.5,2,3,3,1,15*pi/280]);
gauss_para = ([0.5,0,0,1, 0.5,0,0,3]);
dose = doseFun(gauss_para,xy);

[para,loss]=fit_gauss2d(x,y,N_gaussian,dose);
function [para,loss]=fit_gauss2d(x,y,N_gaussian,dose)
Nz = 1;
isGPU = int32(gpuDeviceCount>0);% if GPU valid, use GPU

dx = abs(x(2)-x(1));
dy = abs(y(2)-y(1));
xy = [x,y];

xy_meas = F(dose);
x_measure = xy_meas(1:128);
y_measure = xy_meas(129:256);
[~,arg_x] = max(x_measure);
[~,arg_y] = max(y_measure);
mux = x(arg_x);
muy = y(arg_y);
S = sum(dose,"all")*dx*dy;
sigma = max(dose,[],"all")/2/S;
gauss_para = repmat([S/N_gaussian,mux,muy,sigma],1,N_gaussian);

% lb = (repmat([0,-20,-20,1e-7,1e-7,0],1,N_gaussian));
% ub = (repmat([1, 20, 20,  10,  10,pi]',1,N_gaussian));
lb = (repmat([0,-64,-64,1e-6],1,N_gaussian));
ub = (repmat([1e8, 64, 64,  30]',1,N_gaussian));

options = optimoptions('lsqcurvefit','display','none','SpecifyObjectiveGradient',true);
[para,loss] = lsqcurvefit(@fun,gauss_para,xy,dose,lb,ub,options);

function [f,j] = fun(gauss_para,xy)
    f = Gauss2D(xy(:,1),xy(:,2),gauss_para,Nz,N_gaussian,isGPU,0);% objective function values at x
    if nargout > 1
        j = Gauss2D(xy(:,1),xy(:,2),gauss_para,Nz,N_gaussian,isGPU,1)';
    end
end
end
%% forward and backward of measure function
function xy_meas = F(in)
% in : 2D image, (128,128)
% xy_meas : (256,1)
x_measure = sum(in,1)';
y_measure = sum(in,2);
xy_meas = [x_measure;y_measure];
end