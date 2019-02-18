function [D] = dist(i, x)
% Aimed at Calculating the Euclidean distances between the i-th object and all objects in x
% -------------------------------------------------------------------------
% Input:
% i - an object (1,n)
% x - data matrix (m,n); m-objects, n-variables
% -------------------------------------------------------------------------
% Output:
% D - Euclidean distance (m,1)

[m,n]=size(x);
D = sqrt(sum((((ones(m,1)*i)-x).^2)'));

if n == 1
    D = abs((ones(m,1) * i - x))';
end

end