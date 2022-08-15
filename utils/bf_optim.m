function [F,J] = bf_optim(para,depth)
    % depth : 1D double array, size (1,n) or (n,1). unit is cm;
    % para  : 1D double array, parameter of bortfeld function
    %         size must be (4*m,1) or (1,4*m)  [range, sigma, slope, Phi]
    %         m = {1,2,3...}, m denotes the number of bragg peaks
    
    F = double(BortfeldFunction(single(depth),single(para),0));% objective function values at x
    if nargout > 1   % two output arguments
        % return gradient
        J = double(BortfeldFunction(single(depth),single(para),1));
    end
end

