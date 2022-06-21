function [x1] = randomFieldSim(sz)
%Hard coded simulation of a randomfield of size sz with matern covariance
%and mean 1. Field is shifted so all values are positive and min(x1) = 0
%   size - size of field [x,y]

%size of image
m = sz(1);
n = sz(2);

%create coordinates
x = linspace(0,1,n);
y = linspace(0,1,m);
[X,Y] = meshgrid(x,y);


loc = [X(:) Y(:)]; %save coordinates
D = squareform(pdist(loc)); %distance matrix
r = matern_covariance(D,1,1,1); %covariance matrix

R = chol(r); %cholezky factorization
mu = ones(m*n,1); %mean = 1
z = randn(m*n,1); %normal displacement, sigma = 1

x1 = mu + R'*z; %random field
x1 = reshape(x1,[m n]); %reshape to mxn
x1 = x1+abs(min(min(x1))); %rescale to positive values
end

