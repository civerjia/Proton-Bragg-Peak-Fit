% help function of BortfeldFunction.mexw64
% author : shuang zhou civerjia@gmail.com
% How to use: BortfeldFunction(depth,para,idx);
% depth : 1D double array, size (1,n) or (n,1). unit is cm;
% para  : 1D double array, parameter of bortfeld function 
%         size must be (4*m,1) or (1,4*m)
%         m = {1,2,3...}, m denotes the number of bragg peaks
% idx   : int scalar, idx = 0, 1, or other number
% output: 1D double array or 2D array
%       : if idx = 0, output is Integrated Depth Dose(IDD), size(n,1)
%       : if idx = other, output is mean gradient of para, size(4*m,1)
%       : if idx = 1, output is jacobian, size(m,4)