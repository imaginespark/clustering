function [Eps] = epsilon(x, k)
% Analytical way of estimating neighborhood radius for DBSCAN
% -------------------------------------------------------------------------
% Input:
% x - data matrix (m,n); m-objects, n-variables
% k - number of objects in a neighborhood of an object
% (minimal number of objects considered as a cluster)
% -------------------------------------------------------------------------
% Output:
% Eps - neighborhood radius, if not known avoid this parameter or put []

[m,n] = size(x);
Eps = ((prod(max(x)-min(x))*k*gamma(.5*n+1))/(m*sqrt(pi.^n))).^(1/n);

end