dcm_name = "data\WaterDose_E110.24_gridMatchCTgrid_N1000000.dcm";
dose = flip(squeeze(dicomread(dcm_name)),3);
info = dicominfo(dcm_name);
dose = double(dose)*info.DoseGridScaling*1e4;%*info.DoseGridScaling
% load('data\pRG2.mat');

idd = squeeze(sum(dose,[1,2]));
[idd_max,argmax] = max(idd);
Nz = round(argmax/2)*2+20;
in = zeros(128,128,Nz);
in(31:91,31:91,:) = dose(:,:,1:Nz);
[xy_meas,zx_meas,zy_meas] = F(in);
%% show measurement
f = figure;
subplot(3,1,1)
plot(xy_meas);xlim([1,256])
title('XY measurement')
subplot(3,1,2)
imagesc(zx_meas);colorbar;colormap('gray');
title('ZX measurement')
subplot(3,1,3)
imagesc(zx_meas);colorbar;colormap('gray');
title('ZY measurement')
exportgraphics(f,'./images/MLSIC_measurement.png','Resolution',600)
%%
x = (((1:128)-64.5)*0.2)';
y = x;
xyz_meas = F2(in);
% Ft2(F2(in));
%%
N_gaussian = 3;
[gauss_para_best,loss] = fit2LayerXY(x,y,N_gaussian,xy_meas);
%% show results for xy
Nz = 2;
% calculate dose CPU version
dose2d = @(gauss_para,x,y) Gauss2D(x,y,gauss_para,Nz,N_gaussian,0,0);
dose_recon = dose2d(gauss_para_best,x,y);
figure
subplot(2,3,1)
imagesc(dose_recon(:,:,1));colorbar;colormap('gray');
subplot(2,3,2)
imagesc(in(:,:,1));colorbar;colormap('gray');
subplot(2,3,3)
imagesc(in(:,:,1) - dose_recon(:,:,1));colorbar;colormap('gray');
subplot(2,3,4)
imagesc(dose_recon(:,:,2));colorbar;colormap('gray');
subplot(2,3,5)
imagesc(in(:,:,2));colorbar;colormap('gray');
subplot(2,3,6)
imagesc(in(:,:,2) - dose_recon(:,:,2));colorbar;colormap('gray');
figure
plot(Fxy(dose_recon));
hold on
plot(xy_meas)
disp(squeeze(sum(dose_recon,[1,2]))-idd(1:2));
%% fit z 2 layer
prior_para = gauss_para_best;
[gauss_para_best_z,loss_z1] = fit2LayerZ(x,y,N_gaussian,zx_meas(:,1),zy_meas(:,1),prior_para);

%% show results for z
Nz = 2;
% calculate dose CPU version
dose2d = @(gauss_para,x,y) Gauss2D(x,y,gauss_para,Nz,N_gaussian,0,0);
dose_recon = dose2d(gauss_para_best_z,x,y);
figure
subplot(2,3,1)
imagesc(dose_recon(:,:,1));colorbar;colormap('gray');
subplot(2,3,2)
imagesc(in(:,:,3));colorbar;colormap('gray');
subplot(2,3,3)
imagesc(in(:,:,3) - dose_recon(:,:,1));colorbar;colormap('gray');
subplot(2,3,4)
imagesc(dose_recon(:,:,2));colorbar;colormap('gray');
subplot(2,3,5)
imagesc(in(:,:,4));colorbar;colormap('gray');
subplot(2,3,6)
imagesc(in(:,:,4) - dose_recon(:,:,2));colorbar;colormap('gray');
figure
plot(Fz2(dose_recon));
hold on
plot([zx_meas(:,1);zy_meas(:,1)])

%% fit z layer by layer
loss_z = zeros(1000,40);
loss_z(1:length(loss),1) = loss;
dose3d_para = zeros(2*N_gaussian*4,40);
dose3d_para(:,1) = gauss_para_best;
prior_para = gauss_para_best;
for i = 1:size(zx_meas,2)
    [gauss_para_best_i,loss_i] = fit2LayerZ(x,y,N_gaussian,zx_meas(:,i),zy_meas(:,i),prior_para);
    prior_para = gauss_para_best_i;
    dose3d_para(:,i+1) = gauss_para_best_i;
    loss_z(1:length(loss_i),i+1) = loss_i;
end
%%
dose2d = @(gauss_para,x,y) Gauss2D(x,y,gauss_para,80,N_gaussian,0,0);
dose_recon3d = dose2d(dose3d_para(:),x,y);
f = figure;
plot(squeeze(sum(dose_recon3d,[1,2])),'-o')
hold on
plot(idd)
xlim([1,80])
legend('Reconstructed IDD','TOPAS Simulation','Location','northwest')
exportgraphics(f,'./images/Recon3DIDD.png','Resolution',600);
f = figure;
subplot(3,1,1)
imagesc(squeeze(dose_recon3d(64,:,:)));colorbar;colormap('gray');
title('Reconstructed')
subplot(3,1,2)
imagesc(squeeze(in(64,:,:)));colorbar;colormap('gray');
title('TOPAS Simulation')
subplot(3,1,3)
imagesc(squeeze(in(64,:,:) - dose_recon3d(64,:,:)));colorbar;colormap('gray');
title('Difference')
exportgraphics(f,'./images/Recon3DSlice.png','Resolution',600);
%% fit 2 layer dose
function [gauss_para_best,loss] = fit2LayerXY(x,y,N_gaussian,xyz_meas)
Nz = 2;
% isGPU = int32(gpuDeviceCount>0);% if GPU valid, use GPU
isGPU = 0;

dx = abs(x(2)-x(1));
dy = abs(y(2)-y(1));
xy = [x,y];

x_measure = xyz_meas(1:128);
y_measure = xyz_meas(129:256);
[xmax,arg_x] = max(x_measure);
[ymax,arg_y] = max(y_measure);
mux = x(arg_x);
muy = y(arg_y);
S1 = sum(x_measure)*dy*dx;
S2 = sum(y_measure)*dx*dy;
sigma1 = xmax/S1;
sigma2 = ymax/S2;
% layer1
gauss_para1 = rand(4*N_gaussian,1);
gauss_para1(1) = S1;
gauss_para1(4) = sigma1;
gauss_para1(2:4:end) = mux;
gauss_para1(3:4:end) = muy;
gauss_para1(5:4:end) = 0.1*S1+0.1*rand(N_gaussian-1,1);
% layer2
gauss_para2 = rand(4*N_gaussian,1);
gauss_para2(1) = S2;
gauss_para2(4) = sigma2;
gauss_para2(2:4:end) = mux;
gauss_para2(3:4:end) = muy;
gauss_para2(5:4:end) = 0.1*S2+0.1*rand(N_gaussian-1,1);

gauss_para = [gauss_para1;gauss_para2];
gauss_para = gauss_para + 0.1*rand(size(gauss_para));
% optimization setting
optimOptions.N_gaussian = N_gaussian;
optimOptions.Niter = 1000;% max iteration
optimOptions.alpha = 1e-1;% step length
optimOptions.beta1 = 0.9;%0.9
optimOptions.beta2 = 0.999;%0.999
optimOptions.epsilon = 1e-8;%1e-8
optimOptions.tol = 1e-6;%1e-6; stopping criterion
optimOptions.lb = reshape(repmat([1e-7,-(12.7),-(12.7),1e-6],Nz,N_gaussian)',[],1);
optimOptions.ub = reshape(repmat([1e8, (12.7), (12.7), 30],Nz,N_gaussian)',[],1);

f = @(gauss_para) Gauss2D(xy(:,1),xy(:,2),gauss_para,Nz,N_gaussian,isGPU,0);% objective function values at x
g = @(gauss_para) Gauss2D(xy(:,1),xy(:,2),gauss_para,Nz,N_gaussian,isGPU,1);

[gauss_para_best,loss] = Adamxy(gauss_para,xyz_meas,optimOptions,f,g,S1+0.005,S2+0.005);

end
%% fit 3D dose
function [gauss_para_best,loss] = fit2LayerZ(x,y,N_gaussian,zx_meas,zy_meas,prior_para)
Nz = 2;
% isGPU = int32(gpuDeviceCount>0);% if GPU valid, use GPU
isGPU = 0;

dx = abs(x(2)-x(1));
dy = abs(y(2)-y(1));
xy = [x,y];

S1 = sum(zx_meas)*dy*dx;
S2 = sum(zy_meas)*dx*dy;
gauss_para = prior_para;
gauss_para(1:4:4*N_gaussian) = gauss_para(1:4:4*N_gaussian).*S1/sum(gauss_para(1:4:4*N_gaussian));
gauss_para(4*N_gaussian+1:4:end) = gauss_para(4*N_gaussian+1:4:end).*S2/sum(gauss_para(4*N_gaussian+1:4:end));
% optimization setting
optimOptions.N_gaussian = N_gaussian;
optimOptions.Niter = 1000;% max iteration
optimOptions.alpha = 1e-4;% step length
optimOptions.beta1 = 0.9;%0.9
optimOptions.beta2 = 0.999;%0.999
optimOptions.epsilon = 1e-8;%1e-8
optimOptions.tol = 1e-6;%1e-6; stopping criterion
optimOptions.lb = reshape(repmat([1e-7,-(12.7),-(12.7),1e-6],Nz,N_gaussian)',[],1);
optimOptions.ub = reshape(repmat([1e8, (12.7), (12.7), 30],Nz,N_gaussian)',[],1);

f = @(gauss_para) Gauss2D(xy(:,1),xy(:,2),gauss_para,Nz,N_gaussian,isGPU,0);% objective function values at x
g = @(gauss_para) Gauss2D(xy(:,1),xy(:,2),gauss_para,Nz,N_gaussian,isGPU,1);
z_meas = [zx_meas;zy_meas];
[gauss_para_best,loss] = Adamz(gauss_para,z_meas,optimOptions,f,g,S1+5e-3,S2+5e-3);

end
%% measure function
function out = Fxy(in)
x_meas = sum(in(:,:,1),1);
y_meas = sum(in(:,:,2),2);
out = [x_meas';y_meas];%(256,1)
end
function out = Fxy_trans(xy_meas)
% xy_meas : (256,1)
x_meas = xy_meas(1:128);
y_meas = xy_meas(129:end);
out = cat(3,repmat(x_meas,1,128)',repmat(y_meas,1,128));
out = out(:);
end

function [zx_meas,zy_meas] = Fz(in)
% in : (128,128,Nz)
space = 4;
zx_meas = zeros(8,size(in,3)/2-1);
zy_meas = zx_meas;
temp = reshape(1:128,space,[]);
temp1 = temp';
zxidx = reshape(temp1,8,[]);
temp = reshape([127:128,1:126],space,[]);
temp1 = temp';
zyidx = reshape(temp1,8,[]);
for i = 3:2:size(in,3)
    zx = sum(in(:,:,i),1);
    zy = sum(in(:,:,i+1),2);
    zx_meas(:,(i+1)/2-1) = sum(zx(zxidx),2);
    zy_meas(:,(i+1)/2-1) = sum(zy(zyidx),2);
end
end
function z = Fz_trans(zx_meas,zy_meas)
% z1_meas : (8,Nz/2-1)
% z2_meas : (8,Nz/2-1)
z = zeros(128,128,size(z1_meas,2)*2+2);
temp = reshape(1:128,space,[]);
temp1 = temp';
zxidx = reshape(temp1,8,[]);
temp = reshape([3:128,1:2],space,[]);
temp1 = temp';
zyidx = reshape(temp1,8,[]);
zx_temp = zeros(128,size(zx_meas,2));
zy_temp = zx_temp;
for i = 1:8
    zx_temp(zxidx(i,:),:) = repmat(zx_meas(i,:),16,1);
    zy_temp(zyidx(i,:),:) = repmat(zy_meas(i,:),16,1);
end
zx = repmat(reshape(zx_temp,1,128,[]),128,1,1);
zy = repmat(reshape(zy_temp,128,1,[]),1,128,1);
z(:,:,3:2:end) = zx;
z(:,:,4:2:end) = zy;
end
function z_meas = Fz2(in)
% in : (128,128,Nz)
space = 4;
zx_meas = zeros(8,size(in,3)/2);
zy_meas = zx_meas;
temp = reshape(1:128,space,[]);
temp1 = temp';
zxidx = reshape(temp1,8,[]);
temp = reshape([127:128,1:126],space,[]);
temp1 = temp';
zyidx = reshape(temp1,8,[]);
for i = 1
    zx = sum(in(:,:,i),1);
    zy = sum(in(:,:,i+1),2);
    zx_meas(:,(i+1)/2) = sum(zx(zxidx),2);
    zy_meas(:,(i+1)/2) = sum(zy(zyidx),2);
end
z_meas = [zx_meas(:);zy_meas(:)];
end
function z = Fz_trans2(z_meas)
% z_meas : (16,1)
zx_meas = z_meas(1:8);
zy_meas = z_meas(9:end);
z = zeros(128,128,2);
space = 4;
temp = reshape(1:128,space,[]);
temp1 = temp';
zxidx = reshape(temp1,8,[]);
temp = reshape([3:128,1:2],space,[]);
temp1 = temp';
zyidx = reshape(temp1,8,[]);
zx_temp = zeros(128,size(zx_meas,2));
zy_temp = zx_temp;
for i = 1:8
    zx_temp(zxidx(i,:),:) = repmat(zx_meas(i,:),16,1);
    zy_temp(zyidx(i,:),:) = repmat(zy_meas(i,:),16,1);
end
zx = repmat(reshape(zx_temp,1,128,[]),128,1,1);
zy = repmat(reshape(zy_temp,128,1,[]),1,128,1);
z(:,:,1) = zx;
z(:,:,2) = zy;
z = z(:);
end
function [xy_meas,z1_meas,z2_meas] = F(in)
    x_meas = sum(in(:,:,1),1);
    y_meas = sum(in(:,:,2),2);
    xy_meas = [x_meas';y_meas];%(256,1)
    [z1_meas,z2_meas] = Fz(in);
end
function out = Ft(xy_meas,z1_meas,z2_meas)
    Nz = size(z1_meas,2)*2+2;
    out = zeros(128,128,Nz);
    x_meas = xy_meas(1:128);
    y_meas = xy_meas(129:end);
    out(:,:,1) = repmat(x_meas,1,128)';
    out(:,:,2) = repmat(y_meas,1,128);

    temp = reshape(1:128,4,[]);
    temp1 = temp';
    idx = reshape(temp1,8,[]);
    z1_temp = zeros(128,size(z1_meas,2));
    z2_temp = z1_temp;
    for i = 1:8
        z1_temp(idx(i,:),:) = repmat(z1_meas(i,:),16,1);
        z2_temp(idx(i,:),:) = repmat(z2_meas(i,:),16,1);
    end
    z1 = repmat(reshape(z1_temp,1,128,[]),128,1,1);
    z2 = repmat(reshape(z2_temp,128,1,[]),1,128,1);
    out(:,:,3:2:end) = z1;
    out(:,:,4:2:end) = z2;
end
function [xyz_meas] = F2(in)
    x_meas = sum(in(:,:,1),1);
    y_meas = sum(in(:,:,2),2);
    xy_meas = [x_meas';y_meas];%(256,1)
    [zx_meas,zy_meas] = Fz(in);
    xyz_meas = [xy_meas;zx_meas(:);zy_meas(:)];
end
function out = Ft2(xyz_meas)
    N = length(xyz_meas);
    xy_meas = xyz_meas(1:256);
    zx_meas = reshape(xyz_meas(257:256+(N-256)/2),8,[]);
    zy_meas = reshape(xyz_meas(257+(N-256)/2:end),8,[]);

    Nz = size(zx_meas,2)*2+2;
    out = zeros(128,128,Nz);
    x_meas = xy_meas(1:128);
    y_meas = xy_meas(129:end);
    out(:,:,1) = repmat(x_meas,1,128)';
    out(:,:,2) = repmat(y_meas,1,128);

    temp = reshape(1:128,4,[]);
    temp1 = temp';
    idx = reshape(temp1,8,[]);
    z1_temp = zeros(128,size(zx_meas,2));
    z2_temp = z1_temp;
    for i = 1:8
        z1_temp(idx(i,:),:) = repmat(zx_meas(i,:),16,1);
        z2_temp(idx(i,:),:) = repmat(zy_meas(i,:),16,1);
    end
    z1 = repmat(reshape(z1_temp,1,128,[]),128,1,1);
    z2 = repmat(reshape(z2_temp,128,1,[]),1,128,1);
    out(:,:,3:2:end) = z1;
    out(:,:,4:2:end) = z2;
    out = out(:);
end
%% optimizer
function [theta_best,loss] = Adamxy(para,xyz_meas,optimOptions,f,g,IDD1,IDD2)
    N_gaussian = optimOptions.N_gaussian;
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
    loss(1) = norm(Fxy(f(theta_tm1)) - xyz_meas,'fro');
    for t = 2:T
        % get gradient = jacobian*error
        g_t = 2*g(theta_tm1)*Fxy_trans(Fxy(f(theta_tm1)) - xyz_meas);
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
        % concentricity 
        theta_t(2:4:end) = mean(theta_t(2:4:end));
        theta_t(3:4:end) = mean(theta_t(3:4:end));
        % force the total dose equal IDD at current depth
        theta_t(1:4:4*N_gaussian) = theta_t(1:4:4*N_gaussian).*IDD1/sum(theta_t(1:4:4*N_gaussian));
        theta_t(4*N_gaussian+1:4:end) = theta_t(4*N_gaussian+1:4:end).*IDD2/sum(theta_t(4*N_gaussian+1:4:end));

        

        
        theta_tm1 = theta_t;
        m_tm1 = m_t;
        v_tm1 = v_t;
            
        pred = Fxy(f(theta_t));
        loss(t) = norm((pred - xyz_meas),'fro');

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
function [theta_best,loss] = Adamz(para,xyz_meas,optimOptions,f,g,IDD1,IDD2)
    N_gaussian = optimOptions.N_gaussian;
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
    loss(1) = norm(Fz2(f(theta_tm1)) - xyz_meas,'fro');
    for t = 2:T
        % get gradient = jacobian*error
        g_t = 2*g(theta_tm1)*Fz_trans2(Fz2(f(theta_tm1)) - xyz_meas);
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
        % concentricity 
        theta_t(2:4:end) = mean(theta_t(2:4:end));
        theta_t(3:4:end) = mean(theta_t(3:4:end));
        % force the total dose equal IDD at current depth
        theta_t(1:4:4*N_gaussian) = theta_t(1:4:4*N_gaussian).*IDD1/sum(theta_t(1:4:4*N_gaussian));
        theta_t(4*N_gaussian+1:4:end) = theta_t(4*N_gaussian+1:4:end).*IDD2/sum(theta_t(4*N_gaussian+1:4:end));

        
        theta_tm1 = theta_t;
        m_tm1 = m_t;
        v_tm1 = v_t;
            
        pred = Fz2(f(theta_t));
        loss(t) = norm((pred - xyz_meas),'fro');

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
