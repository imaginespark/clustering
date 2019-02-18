function [class, type] = dbscan(dataPath, k, Eps)
% Aimed at clustering the data with Density-Based Scan Algorithm with Noise (DBSCAN)
% -------------------------------------------------------------------------
% Input:
% dataPath - the file path of data
% k - number of objects in a neighborhood of an object
% (minimal number of objects considered as a cluster)
% Eps - neighborhood radius, if not known avoid this parameter or put []
% -------------------------------------------------------------------------
% Output:
% class - vector specifying assignment of the i-th object to certain
% cluster (m,1)
% type - vector specifying type of the i-th object
% (core: 1, border: 0, outlier: -1)

datafile=[dataPath,'/points.mat'];
load(datafile);
x = points;
[m,n]=size(x);

tic;
if nargin<3 | isempty(Eps)
    [Eps] = epsilon(x, k);
end


x=[[1:m]' x];
[m, n] = size(x);
type = zeros(1,m);
no = 1;
touched = zeros(m,1);


for i = 1 : m
    if touched(i) == 0
        ob = x(i,:);
        D = dist(ob(2:n), x(:,2:n));
        index = find(D <= Eps);
        
        if length(index) > 1 & length(index) < k+1
            type(i) = 0;
            class(i) = 0;
        end
        
        if length(index) == 1
            type(i) = -1;             %为什么是-1呢
            class(i) = -1;             %为什么是-1呢
            touched(i) = 1;
        end
        
        if length(index) >= k+1
            type(i) = 1;
            class(index) = ones(length(index),1)*max(no);
            
            while ~isempty(index)
                ob = x(index(1),:);
                touched(index(1)) = 1;
                index(1) = [];
                D = dist(ob(2:n), x(:,2:n));
                i1 = find(D <= Eps);
                
                if length(i1)>1
                    class(i1) = no;
                    if length(i1) >= k+1;
                        type(ob(1)) = 1;
                    else
                        type(ob(1)) = 0;
                    end
                    
                    for i = 1 : length(i1)
                        if touched(i1(i)) == 0
                            touched(i1(i)) = 1;
                            index = [index i1(i)];
                            class(i1(i)) = no;
                        end
                    end
                end
            end
            no = no + 1;
        end
    end
end

i1 = find(class == 0);
class(i1) = -1;
type(i1) = -1;
title('DBSCAN')
drawshapes(points, class, m);
toc;
end



