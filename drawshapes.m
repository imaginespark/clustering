function drawshapes(points,class,ND)
% Aimed at ploting the points
% -------------------------------------------------------------------------
% Input:
% points - the data
% class - colored point mark
% ND - the number of points

%颜色变换步长
colorindent = 100/7;
%形状
shapes='*o^+sx<d.^ph>dv';
hold on;
for i=1:ND
    if (class(i)>0)
        %colr,colg,colb是画图的RGB颜色值
        color=colorindent*(class(i)-1);
        color_red= mod(color,10);
        color=fix(color/10);
        color_green=mod(color,10);
        color=fix(color/10);
        color_blue=mod(color,10);
        v=0.1*[color_red,color_green,color_blue];
        %选画图形状
        shapeindex= mod(class(i),14)+1;
        
        plot(points(i,1),points(i,2),shapes(shapeindex),'MarkerSize',5,'MarkerFaceColor',v,'MarkerEdgeColor',v);
    else
        %离群点用红色o画出
        plot(points(i,1),points(i,2),'o','MarkerSize',5,'MarkerEdgeColor','r');
    end
    set(gca,'box','on')
end
end