addpath('../data/');
addpath('..');
%% fit and plot IBA Zebra data
load('../data/zebra_idd.mat');
f = figure;
c = lines();
for i = 1:10
    idd_i = gains(:,i);
    idd_i = idd_i/max(idd_i);
    z = depth*0.1;%unit is cm
    strict = 0;% strict = 1, max number of bp is fixed to num_bp
    num_bp = 2;
    [x,idd_o] = precise_fit(z,idd_i,num_bp,strict);
    plot(z,idd_i,'Color',c(i,:),'LineWidth',1);hold on
    plot(z,idd_o,'-.','Color',c(i,:),'LineWidth',1)
end
ylabel('Dose(a.u.)')
xlabel('Depth(cm)')
grid on;
grid minor;
exportgraphics(f,'../images/Zebra_fit.png','Resolution',600)
%% fit and plot a few samples of head phantom data measured by MLSIC
load('../data/head_idd.mat');
idx = [214 454 1587 1724 4981];
f = figure;
c = lines();
cnt = 1;
for i = idx
    idd_i = z2d(i,:)'; %214 454 1587 1724 4981
    num_bp = 3;
    x = fast_fit(z',idd_i,num_bp);
    idd_o = bf_mex(z',x,'idd');
    plot(z,idd_i,'Color',c(cnt,:),'LineWidth',1);hold on
    plot(z,idd_o,'-.','Color',c(cnt,:),'LineWidth',1)
    cnt = cnt + 1;
end
ylabel('Dose(a.u.)')
xlabel('Depth(cm)')
grid on;
grid minor;
exportgraphics(f,'../images/IDD_fit.png','Resolution',600)
%% time benchmark 1 about 160s @i9-9900k
load('./data/head_idd.mat');
num_bp = 2;
tic;
x_out = zeros(num_bp*4,size(z2d,1));
parfor i = 1:size(z2d,1)
    idd_i = z2d(i,:)';
    x = fast_fit(z',idd_i,num_bp);
    x_out(:,i) = x;
end
toc;
%% time benchmark 2 about 123.8s @i9-9900k
load('../data/electron_idd.mat');
num_bp = 2;
tic;
x_out = zeros(num_bp*4,size(z2d,1));
parfor i = 1:size(z2d,1)
    idd_i = z2d(i,:)';
    x = fast_fit(z',idd_i,num_bp);
    x_out(:,i) = x;
end
toc;
%% show results
load('../data/electron_idd.mat');
num_bp = 2;
c = lines();
for i = 1674%1:size(z2d,1)
    idd_i = z2d(i,:)';
    x = fast_fit(z',idd_i,num_bp);
    idd_o = bf_mex(z',x,'idd');
    plot(z,idd_i,'Color',c(mod(i,256),:),'LineWidth',1);hold on
    plot(z,idd_o,'-.','Color',c(mod(i,256),:),'LineWidth',1)
end
%% convert results to image
zrange = medfilt1(x_out(1,:),3)';
[vo1] = protonspot2image(x_para,y_para,zrange);
zrange = medfilt1(x_out(5,:),3)';
[vo2] = protonspot2image(x_para,y_para,zrange);
imagesc([vo1,vo2])
