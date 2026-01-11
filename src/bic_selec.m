function [Theta, Z, edgest, lambda] = bic_selec(X,N,lambdav,pa,p, ...
    maxiter,rho,alpha,tol,varrho,niteradap,adapflag,printbic)
%
% nsamp = number of samples for stability selection
%
m= pa/p;
matr = m;
lambdanum = length(lambdav);
bic = zeros(lambdanum,1);
for nlambda = 1:1:lambdanum
    %
    lambda = lambdav(nlambda);
    lamsc = matr;
    [Theta, Z, edgest, iter] = optimize_admm_ma(X,pa,p, ...
        maxiter,rho,lambda,alpha,tol,varrho,niteradap,adapflag);
    sum(sum(sum(abs(Z) > 0)));
    eps = 10^(-10);
    % calculate BIC
    costo = -log(abs(det(Z))+eps) + (trace(X*Z));
    if (costo < -10^6)
        costo = 10^6;
    end
    addit2 = (log(N)/N)*(sum(sum(abs(Z) > 0))/2);
    costo = costo + addit2;
    %
    %
    bic(nlambda) = costo;
end
if (printbic == 1)
    bic'
end
index1bic = min(find(bic(:) == min(bic(:))));
lambda =lambdav(index1bic);
%
% run with selected lambda
%
rho = 2; % ADMM parameter
maxiter = 200; % maximum number of iterations allowed for optimization
tol = 0.0001;
varrho = 1;
[Theta, Z, edgest, iter] = optimize_admm_ma(X,pa,p, ...
    maxiter,rho,lambda,alpha,tol,varrho,niteradap,adapflag);
%
end