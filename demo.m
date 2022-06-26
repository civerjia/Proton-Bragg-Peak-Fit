
load('zebra_idd.mat');
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
exportgraphics(f,'Zebra_fit.png','Resolution',600)
%%
load('head_idd.mat');
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
exportgraphics(f,'IDD_fit.png','Resolution',600)
%% time benchmark
load('head_idd.mat');
num_bp = 2;
tic;
x_out = zeros(num_bp*4,size(z2d,1));
parfor i = 1:size(z2d,1)
    idd_i = z2d(i,:)';
    x = fast_fit(z',idd_i,num_bp);
    x_out(:,i) = x;
end
toc;