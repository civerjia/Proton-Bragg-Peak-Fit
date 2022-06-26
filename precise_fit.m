function [x,idd_o,resnorm] = precise_fit(depth,idd_i,num_bp,strict)
    fun = @(para,z) bf_mex(z,para,'idd');
    z_max = max(depth);
    [v,i] = findpeaks(medfilt1(idd_i,3),'MinPeakProminence',0.005,...
        'MinPeakDistance',2,'MinPeakWidth',1,'Annotate','extents');
    [v,idx] = sort(v,'descend');
    i = i(idx);
    if length(v) < num_bp
        [v,i] = maxk(abs(diff(medfilt1(idd_i,3))),num_bp);
    elseif strict~=1
        num_bp = length(v);
    end
    %para0 = [zr,0.07*zr,1e-3,zv*zr*1e-2];% prior
    x0 = zeros(4*num_bp,1);
    lb = repmat([0,0,-10,0]',1,num_bp);
    ub = repmat([1.2*z_max 10, 10, 10]',1,num_bp);
    zv = v(1:num_bp);
    zr = depth(i(1:num_bp));
    x0(1:4:end) = zr;
    x0(2:4:end) = 0.07*zr;
    x0(3:4:end) = 1e-3;
    x0(4:4:end) = zv.*zr.*1e-2;
    options = optimoptions('lsqcurvefit','display','none');
    [x,resnorm] = lsqcurvefit(fun,x0,depth,idd_i,lb,ub,options);
    idd_o = fun(x,depth);
end

