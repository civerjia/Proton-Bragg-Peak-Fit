function output = bf_mex(depth,para,func_name)
    % depth : 1D double array, size (1,n) or (n,1). unit is cm;
    % para  : 1D double array, parameter of bortfeld function
    %         size must be (4*m,1) or (1,4*m)  [range, sigma, slope, Phi]
    %         m = {1,2,3...}, m denotes the number of bragg peaks
    if strcmp(func_name,'idd')
        % return IDD
        output = double(BortfeldFunction(single(depth),single(para),0));
    elseif strcmp(func_name,'jacobian')
        % return mean gradient
        output = double(BortfeldFunction(single(depth),single(para),1));
    else
        error('Undefined function name. Plz use idd or jacobian');
    end
end

