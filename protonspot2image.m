function [vo] = protonspot2image(x_para,y_para,zrange)
    xloc = linspace(-100,100,100)-1;
    yloc = linspace(-100,100,100)-2;
    [Xloc,Yloc] = meshgrid(xloc,yloc);
    xi = 2*x_para(1:end-1,2)-129;
    yi = 2*y_para(1:end-1,2)-129;
    vi = double(zrange(1:end-1));
    F = scatteredInterpolant(xi,yi,vi);
    vo = F(Xloc,Yloc);
end

