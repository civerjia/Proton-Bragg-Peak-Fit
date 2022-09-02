x = (((1:128)-64.5)*0.2)';
y = x;
Nz = 1;
N_gaussian = 2;
% calculate dose CPU version
dose2d = @(gauss_para,xy) Gauss2D(xy(:,1),xy(:,2),gauss_para,Nz,N_gaussian,0,0);
% gauss_para = ([0.5,-2,-3,1,3,45*pi/180, 0.5,2,3,3,1,15*pi/280]);
gauss_para = ([0.5,0,0,1, 0.5,0,0,3]);
dose = dose2d(gauss_para,xy);
dose_n = dose + 1e-1*max(dose,[],'all')*(rand(128)-0.5);
xy_meas = F(dose_n);
img = Ft(xy_meas);
[gauss_para_best,loss] = fit2Ddose2(x,y,N_gaussian,xy_meas);
dose_recon = dose2d(gauss_para_best,xy);
% mse(dose_recon,dose2d(gauss_para,xy));
figure
plot(dose_n(64,:))
hold on
plot(dose(64,:))
plot(dose_recon(64,:))
figure
plot(xy_meas)
hold on
plot(F(dose_recon))
err = abs(dose2d(gauss_para,xy)- dose2d(gauss_para_best,xy));
max(err,[],'all')./max(dose,[],'all')
function [gauss_para_best,loss] = fit2Ddose2(x,y,N_gaussian,xy_meas)
Nz = 1;
isGPU = int32(gpuDeviceCount>0);% if GPU valid, use GPU

dx = abs(x(2)-x(1));
dy = abs(y(2)-y(1));
xy = [x,y];

x_measure = xy_meas(1:128);
y_measure = xy_meas(129:256);
[xmax,arg_x] = max(x_measure);
[ymax,arg_y] = max(y_measure);
mux = x(arg_x);
muy = y(arg_y);
S1 = sum(x_measure)*dy*dx;
S2 = sum(y_measure)*dx*dy;
S = (S1+S2)/2;
sigma = (xmax + ymax)/2/S;
gauss_para = repmat([S/N_gaussian,mux,muy,sigma],1,N_gaussian)';
gauss_para = gauss_para + 0.1*rand(size(gauss_para));
optimOptions.Niter = 1000;
optimOptions.alpha = 1e-2;
optimOptions.beta1 = 0.9;%0.9
optimOptions.beta2 = 0.999;%0.999
optimOptions.epsilon = 1e-8;%1e-8
optimOptions.tol = 1e-6;%1e-6;
optimOptions.lb = (repmat([0,-0,-0,1e-6],1,N_gaussian))';
optimOptions.ub = (repmat([1e8, 0, 0, 30],1,N_gaussian))';

f = @(gauss_para) Gauss2D(xy(:,1),xy(:,2),gauss_para,Nz,N_gaussian,isGPU,0);% objective function values at x
g = @(gauss_para) Gauss2D(xy(:,1),xy(:,2),gauss_para,Nz,N_gaussian,isGPU,1);

[gauss_para_best,loss] = Adam(gauss_para,xy_meas,optimOptions,f,g);


end
%% forward and backward of measure function
function xy_meas = F(in)
% in : 2D image, (128,128)
% xy_meas : (256,1)
x_measure = sum(in,1)';
y_measure = sum(in,2);
xy_meas = [x_measure;y_measure];
end
function img = Ft(xy_meas)
% xy_meas : (256,1)
% img : 2D image, (128,128)
x_measure = xy_meas(1:128)';
y_measure = xy_meas(129:256);
img = repmat(x_measure,128,1) + repmat(y_measure,1,128);
img = img(:);
end
%% optimizer
function [theta_best,loss] = Adam(para,measure,optimOptions,f,g)

    T = optimOptions.Niter;
    alpha = optimOptions.alpha;
    beta1 = optimOptions.beta1;%0.9
    beta2 = optimOptions.beta2;%0.999
    epsilon = optimOptions.epsilon;%1e-8
    tol = optimOptions.tol;%1e-6;
    lb = optimOptions.lb;% lower bound
    ub = optimOptions.ub;% upper bound
    loss = zeros(T,1);
    m_tm1 = 0;
    v_tm1 = 0;
    theta_tm1 = para;
    
    theta_best = para;
    loss_best = 1e9;
    loss(1) = norm((F(f(theta_tm1)) - measure),'fro');
    for t = 2:T
        % get gradient = jacobian*error
        g_t = 2*g(theta_tm1)*Ft(F(f(theta_tm1)) - measure);
        % Update biased first moment estimate
        m_t = beta1*m_tm1 + (1-beta1)*g_t;
        % Update biased second raw moment estimate
        v_t = beta2*v_tm1 + (1-beta2)*g_t.^2;
        % Compute bias-corrected first moment estimate
        m_t_hat = m_t / (1-beta1^(t-1));
        % Compute bias-corrected second raw moment estimate
        v_t_hat = v_t / (1-beta2^(t-1));
        % Update parameters
        theta_t = theta_tm1 - alpha*m_t_hat./(sqrt(v_t_hat)+epsilon);
        
        % constrain
        theta_t(theta_t < lb) = lb(theta_t < lb);
        theta_t(theta_t > ub) = ub(theta_t > ub);
        
        theta_tm1 = theta_t;
        m_tm1 = m_t;
        v_tm1 = v_t;
            
        pred = F(f(theta_t));
        loss(t) = norm((pred - measure),'fro');

        if mod(t,100) == 0
            alpha = alpha*0.9;
        end
        
        if loss(t) < loss_best
           loss_best = loss(t);
           theta_best = theta_t;
        end
        if (abs(loss(t) - loss(t-1)) < tol)
            loss = loss(1:t);
            break;
        end
    end
    
end