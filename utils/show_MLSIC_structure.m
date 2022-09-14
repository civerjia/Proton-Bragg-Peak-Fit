dz = 1.6;% PCB thickness mm
gap = 2.91;% WET thickness
objCenter = [0;0;-dz/2];
objDim3 = [304.8;304.8;dz];% mm
color = 'w';
facealpha = 1;
fontsize = 10;
f = figure;
scatter3(0,0,0,'+');
% HV board
drawCube(objDim3,objCenter,fontsize,facealpha,color,0);
% XY board
drawCube(objDim3,objCenter-[0;0;gap],fontsize,facealpha,color,0);
drawXYstrip(-gap,dz,1,'w');
Nz = 4;
for i = 2:2:Nz
    % HV board
    drawCube(objDim3,objCenter-[0;0;i*gap],fontsize,facealpha,color,0);
    % Z board
    drawCube(objDim3,objCenter-[0;0;(i+1)*gap],fontsize,facealpha,color,0);
    drawZstrip(-(i+1)*gap,dz,1,'w');
end
% HV board
drawCube(objDim3,objCenter-[0;0;(Nz+2)*gap],fontsize,facealpha,color,0);
axis off;
function drawXYstrip(z,dz,facealpha,color)
    X = 256*([0 1 1 0]-0.5);
    Y = 1.874*([0 0 1 1]-0.5);
    Z = z*ones(1,4);
    hold on;
    for i = 1:128
        fill3(X,Y+(i-64)*2-1,Z,color,'FaceAlpha',facealpha);hold on;%X strip
        fill3(Y+(i-64)*2-1,X,Z-dz,color,'FaceAlpha',facealpha);hold on;%Y strip
    end
end
function drawZstrip(z,dz,facealpha,color)
    X = 256*([0 1 1 0]-0.5);
    Y = 7.8*([0 0 1 1]-0.5);
    Ys = 3.8*([0 0 1 1]-0.5);
    Z = z*ones(1,4);
    hold on;
    for i = 1:32
        fill3(X,Y+(i-16)*2*4-4,Z,color,'FaceAlpha',facealpha);hold on;%ZX strip
        if i == 1
            fill3(Ys+15.5*8+2,X,Z-dz,color,'FaceAlpha',facealpha);hold on;%ZY strip
            fill3(Ys-(15.5*8+2),X,Z-dz,color,'FaceAlpha',facealpha);hold on;%ZY strip
        else
            fill3(Y+(i-17)*2*4,X,Z-dz,color,'FaceAlpha',facealpha);hold on;%ZY strip
        end
    end
end
function drawCube(objDim3,objCenter,fontsize,facealpha,color,showLabel)
    
    % face: '$Down$','$Top$','$Rear$','$Front$','$Left$','$Right$'
    X = objDim3(1)*([0 1 1 0; 0 1 1 0; 0 0 0 0; 1 1 1 1; 0 1 1 0; 0 1 1 0]-0.5) + objCenter(1);
    Y = objDim3(2)*([0 0 1 1; 0 0 1 1; 0 1 1 0; 0 1 1 0; 0 0 0 0; 1 1 1 1]-0.5) + objCenter(2);
    Z = objDim3(3)*([0 0 0 0; 1 1 1 1; 0 0 1 1; 0 0 1 1; 0 0 1 1; 0 0 1 1]-0.5) + objCenter(3);
%     color = spring(6);
    hold on;
    for i = 1:6
        fill3(X(i,:),Y(i,:),Z(i,:),color,'FaceAlpha',facealpha);hold on;%color(i,:)
    end
    if showLabel
        % label 6 faces
        str = {'$Down$','$Top$','$Rear$','$Front$','$Left$','$Right$'};
        text(mean(X,2),mean(Y,2),mean(Z,2),str,'Color','black','FontSize',fontsize,'Interpreter','latex');
        % label 8 corners
        corner = 0.5*[-1  1 -1  1 -1  1 -1  1;
            1  1 -1 -1  1  1 -1 -1;
            1  1  1  1 -1 -1 -1 -1];
        vertex = objDim3.*corner + objCenter;
        str = {'$c_1$','$c_2$','$c_3$','$c_4$','$c_5$','$c_6$','$c_7$','$c_8$'};
        text(vertex(1,:),vertex(2,:),vertex(3,:),str,'Color','black','FontSize',fontsize,'Interpreter','latex');
    end
end

