function [x_best,idd_o] = fit_adam(depth,idd_i,num_bp)
% idd_i = z2d(4981,:)'; % 214 454 1587 1724 4981 6342
% depth = z';
% num_bp = 3;

z_max = max(depth);
[v,i] = maxk(abs(diff(medfilt1(idd_i,3))),num_bp);
%para0 = [zr,0.07*zr,1e-3,zv*zr*1e-2];% prior
x0 = zeros(4*num_bp,1);
lb = repmat([0,0,-10,0],1,num_bp)';
ub = repmat([1.2*z_max 10, 10, 10],1,num_bp)';
zv = v;
zr = depth(i);
x0(1:4:end) = zr;
x0(2:4:end) = 0.07*zr;
x0(3:4:end) = 1e-3;
x0(4:4:end) = zv.*zr.*1e-2;

optimOptions.Niter = 2000;
optimOptions.alpha = 1e-2;% learning rate
optimOptions.beta1 = 0.9;%0.9
optimOptions.beta2 = 0.999;%0.999
optimOptions.epsilon = 1e-8;% avoid 1/0
optimOptions.tol = 1e-6;% stopping criterion
optimOptions.lb = lb;% lower bound
optimOptions.ub = ub;% upper bound

f = @(x) bf_mex(depth,x,'idd');
g = @(x) bf_mex(depth,x,'jacobian');
[x_best,loss] = AdamOptim(x0,idd_i,optimOptions,f,g);
idd_o = f(x_best);
% figure;
% plot(f(x_best))
% hold on
% plot(idd_i)
end