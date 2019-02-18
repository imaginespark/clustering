function dpeak_auto(dataPath,percent,kernel,k)
% Aimed at clustering the data with Density Peak Algorithm (DPeak) automatically
% -------------------------------------------------------------------------
% Input:
% dataPath - the file path of data
% percent - average percentage of neighbours
% kernel - gaussian or cutoff kernel
% k - number of clusters

tic;

disp('Description of distance.mat: [i, j, dist(i,j)]')

%加载distances矩阵
datafile = [dataPath,'/distances.mat'];
load(datafile);

xx = distances;

%ND或NL即为xx所表达的样本总数
ND = max(xx(:,2));
NL = max(xx(:,1));
if (NL > ND)
    ND = NL;
end

%计算xx的总行数，即为n*(n-1)/2
N = size(xx,1);

%完成距离矩阵初始化
if (exist([dataPath,'distMat.mat'],'file'))
    load([dataPath,'distMat.mat']);
    dist = distMat;
else
    for i = 1 : ND              %初始化dist 矩阵，所有元素全为0
        for j = 1 : ND
            dist(i,j) = 0;
        end
    end
    for i = 1 : N                   %根据xx矩阵，依次回填dist矩阵
        ii = xx(i,1);
        jj = xx(i,2);
        dist(ii,jj) = xx(i,3);
        dist(jj,ii) = xx(i,3);
    end
    save distance_matrix dist;
end

%计算截距dc
fprintf('average percentage of neighbours (hard coded): %5.6f\n', percent);

%计算前N*percent/100的个数
position = round(N*percent/100);

%将xx第3列，即按距离进行从小到大排序
sda = sort(xx(:,3));%get distance of all points

%处在position位置的元素数值，即是cutoff distance
dc = sda(position); %dc is a  distance from the  paires of some point which index in dataset is equal to value of variable position

fprintf('Computing Rho with gaussian kernel of radius: %5.6f\n', dc);

%计算局部密度
for i = 1 : ND
    rho(i) = 0.;
end

switch(kernel)
    case 'gaussian'
        for i = 1 : ND-1
            for j = i+1 : ND
                % i与j距离与dc相差越大,exp(-(dist(i,j)/dc)*(dist(i,j)/dc)值越小，
                % 最终累积得到rho(i)越小表示与样本点i距离较近的点少，否则表示样本点i周边有很多距离较近的点
                rho(i) = rho(i) + exp(-(dist(i,j)/dc) * (dist(i,j)/dc));
                rho(j) = rho(j) + exp(-(dist(i,j)/dc) * (dist(i,j)/dc));
            end
        end
    case 'cutoff'   % 论文里公式1的卡方函数
        neibors = zeros(ND,1);
        for i = 1 : ND-1
            count = 0;
            for j = 1 : ND
                if (i ~= j)
                    if (dist(i,j)<dc)
                        count = count + 1;
                        rho(i) = rho(i) + 1.;
                        neibors(i,count) = j;
                    end
                end
            end
        end
    otherwise
        disp('please input the correct kernel')
end

%取出dist矩阵中值最大元素
maxd = max(max(dist));

%rho_sorted是倒序排序后的向量，ordrho是rho_sorted各元素在原向量rho中的位置,即下标向量,
%即ordrho(1)是密度最大的样本点的位置，rho_sorted(1)是密度值，ordrho(i)是密度值排第i位的样本点
[rho_sorted, ordrho]=sort(rho, 'descend');

%计算距离
for ii = 2 : ND
    delta(ordrho(ii)) = maxd;
    for jj = 1 : ii-1
        if(dist(ordrho(ii),ordrho(jj)) < delta(ordrho(ii)))
            delta(ordrho(ii)) = dist(ordrho(ii), ordrho(jj));     %取距离的最小值
            nneigh(ordrho(ii)) = ordrho(jj);
        end
    end
end
delta(ordrho(1)) = max(delta(:));  %密度值最大的点对应的距离值


%方案一：绘制决策图（ρ-δ）
disp('Description of decision_graph: [density, delta]')

fid = fopen('decision_graph', 'w');
for i = 1 : ND
    fprintf(fid, '%6.2f %6.2f\n', rho(i), delta(i));
end

%方案二：绘制决策图（n-γ）
for i=1:ND
  ind(i)=i;
  gamma(i)=rho(i)*delta(i);
end

figure(1)
plot(rho(:),delta(:),'o','MarkerSize',3,'MarkerFaceColor','k','MarkerEdgeColor','k');
title('Decision Graph','FontSize',15.0)
xlabel('\rho')
ylabel('\delta')

%统计聚类中心个数
NCLUST = k;
for i = 1:ND
  cl(i) = -1;   %cl为分组标志数组，cl(i)=j表示第i个数据点分组到第j个cluster
end

[B, Index] = sort(gamma, 'descend');
disp(Index)

figure(2)
plot(ind(:),B(:),'o','MarkerSize',3,'MarkerFaceColor','k','MarkerEdgeColor','k');
title('Decision Graph','FontSize',15.0)
xlabel('n')
ylabel('\gamma')

% cl是每个数据点的所属类别号
% icl是所有聚类中心的序号
icl = Index(1:k);
cl(Index(1:k)) = 1:k;

%分组
disp('Performing assignation')
for i = 1 : ND
    if (cl(ordrho(i)) == -1)
        cl(ordrho(i)) = cl(nneigh(ordrho(i)));   %密度值不小于i点且距离i点最近
    end
end

%halo
for i = 1 : ND
    halo(i) = cl(i);
end

if (NCLUST > 1)
    for i = 1 : NCLUST
        bord_rho(i) = 0.;     %边界密度阈值
    end
    for i = 1 : ND-1
        for j = i+1 : ND
            if ((cl(i)~=cl(j)) && (dist(i,j)<=dc))  %距离足够小但不属于同一个cluster的i和j
                rho_aver = (rho(i)+rho(j))/2.;
                if (rho_aver > bord_rho(cl(i)))
                    bord_rho(cl(i)) = rho_aver;
                end
                if (rho_aver > bord_rho(cl(j)))
                    bord_rho(cl(j)) = rho_aver;
                end
            end
        end
    end
    for i = 1 : ND
        if (rho(i) < bord_rho(cl(i)))       %halo定义
            halo(i) = 0;
        end
    end
end

%统计核心点和光晕点个数
for i = 1 : NCLUST
    nc = 0;
    nh = 0;
    for j = 1 : ND
        if (cl(j) == i)
            nc = nc + 1;
        end
        if (halo(j) == i)
            nh = nh + 1;
        end
    end
    fprintf('CLUSTER: %i CENTER: %i ELEMENTS: %i CORE: %i HALO: %i \n', i, icl(i), nc, nh, nc-nh);
end

%可视化
cmap = colormap;
for i = 1 : NCLUST
    ic = int8((i*64.) / (NCLUST*1.));
    figure(3)
    hold on
    plot(rho(icl(i)),delta(icl(i)),'o','MarkerSize',10,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
    xlabel ('\rho');
    ylabel ('\delta');
end

if (exist([dataPath,'/points.mat'],'file'))
    load([dataPath,'/points.mat']);
else
    points = mdscale(dist, 2, 'criterion','metricstress');
    save 'points.mat' points;
end

for i = 1 : ND
    A(i,1) = 0.;
    A(i,2) = 0.;
end

for i = 1 : NCLUST
    nn = 0;
    ic = int8((i*64.)/(NCLUST*1.));
    for j = 1 : ND
        if (cl(j) == i)
            nn = nn + 1;
            A(nn,1) = points(j,1);
            A(nn,2) = points(j,2);
        end
    end
    hold on
    figure(4)
    title ('DPEAK')
    plot(A(1:nn,1),A(1:nn,2),'o','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
end


for i = 1 : ND
    if (halo(i) > 0)
        ic = int8((halo(i)*64.)/(NCLUST*1.));
        hold on
        plot(points(i,1),points(i,2),'o','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
    end
end

fr = fopen('cluster_assignation', 'w');
disp('Description of cluster_assignation: [id, cluster assignation without halo control, cluster assignation with halo control]');
for i = 1 : ND
    fprintf(fr, '%i %i %i\n',i,cl(i),halo(i));
end

result = cl';
save label result

%把peak画出来
for i = 1 : NCLUST
    ic = int8((i*64.)/(NCLUST*1.));
    figure(4)
    plot(points(icl(i),1),points(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','k');
    hold on
end
toc;