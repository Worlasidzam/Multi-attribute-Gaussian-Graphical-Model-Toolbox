function [data] = dataGen(K,N,p,m)
%
% Objective: generate (mp)*N data sequence satisfying precsion matrix K
% 
% Needs functions: none
%
% INPUT
% p= number of nodes
% m = number of random variables per node
% K = pa*pa (pa=mp) precision matrix (inverse covariance); must be invertible
% N = number of data samples
%
% OUTPUT
% data has data samples -- rows are m variables per node with p nodes, 
%   columns are time samples. First m rows of every column are associated
%   with the first node, and so on ...
%
pa = m*p;
phi = sqrtm(K); % sqrtm is a function that computes matrix square root
% i.e.,  K = phi*phi where phi = phi' for K=K'
Htr = inv(phi'); % Will results in cov=Htr*Htr' such that inverse
% covariance cov^-1 = (Htr*Htr')^-1 = (Htr')^-1*Htr^-1 
%    = (inv(phi))^-1*(inv(phi'))^-1 = phi*phi' = K = precision matrix
%
w = randn(pa,N);
data = Htr*w; %data matrix: pa times N
end
