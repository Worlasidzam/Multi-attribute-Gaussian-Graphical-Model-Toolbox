%
% function opt_admm_ma.m
%
% Objective: Perform complete opimization of mult-attribute
% sparse-group lasso problem using ADMM algorithm
%
% Initiated: Dec 17, 2024  ---
% modifications of
% opt_admm_varlamA/jmlr_simul/nips2020/ssp2000/transferDrive
% Incorporate group lasso weights: mxlambda2 instead of lambda2
%
% INPUT:
% lambda1 = alpha*lambda 
% lambda2 = (1-alpha)*m*lambda where m=pa/p
% X = data covariance matrix
% p = number of nodes
% pa = dimension of (square) X
% maxiter = maximum number of iterartions of ADMM
% rho = penalty parameter of ADMM
%
% OUTPUT:
% Theta = estimate of inverse covariance
% Z = variable splitting of Theta
% iter = number of iterations to convergence
% edgest = adjacency matrix (edge-set): 
% 
%
function [Theta, Z, edgest, iter] = ...
    opt_admm_ma(X,pa,p,maxiter,rho,lambda1,lambda2,tol,varrho)
%
%
%  initialize concentration matrix Theta
%
Theta = zeros(pa,pa);
U = Theta; % dual variables in ADMM
Z = Theta; % variable splitting
Theta = diag(1./diag(X)); % diagonal matrix with inverse of diagonal of X
mn = pa/p; % number of attributes
ZOld = Z;
%
crit = 100; % some arbitrary high number
iter = 0; % initial iteration count
%
while (crit > tol)   % e.g. (crit > 0.0001)
    iter = iter + 1;
    %
    % Step 1 of ADMM
    %
    Ax = X - rho*(Z-U);
    [V, D] = eig(Ax); %eigen-decomposition
    for i=1:1:pa
        b = D(i,i);
        D(i,i) = -b + sqrt(b*b+4*rho);
        D(i,i) = D(i,i)/(2*rho);
    end
    Theta = V*D*V';
    %
    % Step 2 of ADMM
    % 
    % 
    Y = Theta + U;
    Zt = Z;
    %
    % element-wise thresholding : lasso
    %
    % kappa = lambda1/rho;
    for i=1:1:pa
        for j=1:1:pa
            if (i == j) % no penalty
                Zt(i,j) = Y(i,j);
            else % apply elementwise regularization
                atn = abs(Y(i,j));
                if (atn > 0)
                    kappa = lambda1(i,j)/rho;
                    ab = max(0,1-(kappa/atn));
                else
                    ab = 0;
                end
                Zt(i,j) = ab*Y(i,j);
            end
        end
    end
    %
    % group thresholding : group-lasso
    %
    %  kappa = lambda2/rho;
    if (pa > p) % more than one attribute
        for i=1:1:p
            for j=1:1:p
                ind1 = [];
                ind2 = [];
                for i1=1:1:mn
                    ind1 = [ind1 (i-1)*mn+i1];
                    ind2 = [ind2 (j-1)*mn+i1];
                end
                atemp = Zt(ind1, ind2);
                if (i == j) % no penalty
                    Z(ind1,ind2) = Zt(ind1,ind2);
                else % apply groupwise regularization
                    atn = norm(atemp,'fro');
                    if (atn > 0)
                        kappa = lambda2(i,j)/rho;
                        ab = max(0,1-(kappa/atn));
                    else
                        ab = 0;
                    end
                    Z(ind1,ind2) = ab*Zt(ind1,ind2);
                end
            end
        end
    else % no group-lasso if number of attributes equal 1
        Z = Zt;
    end
    %
    % Step 3 of ADMM: update dual variables
    %
    U = U + Theta - Z;
    %
    %
    % calculate residuals
    % primal residual
    primres = norm(Z - Theta,'fro');
    avec = [norm(Theta,'fro'), norm(Z,'fro')];
    threshr = pa*tol + tol*max(avec); 
    % dual residual
    dualres = rho*norm(Z - ZOld,'fro');
    threshd = pa*tol + tol*(1/rho)*norm(U,'fro'); 
    if ((primres < threshr) & (dualres < threshd))
        break;
    else
        crit = 100;
    end
    %
    if (varrho == 1)
        %
        % variable penalty parameter rho
        %
        mu = 10;
        if (primres > mu*dualres)
            rho = rho*2;
            U = U/2;
        elseif (dualres > mu*primres)
            rho = rho/2;
            U = 2*U;
        else
            rho = rho;
        end
    end
    %
    %
    %
    if (iter < 5) % do not quit too early; run at least 5 iterations
        crit = 100;
    end
    %
    if (iter > maxiter) % quit if iterated too long
        break;
    end
    % 
    ZOld = Z;
end
%
% DONE WITH ADMM Optimization
%
Thetastat = zeros(p,p); % edge detection
for i=1:1:p
    for j=i+1:1:p
        ind1 = [];
        ind2 = [];
        for i1=1:1:mn
            ind1 = [ind1 (i-1)*mn+i1];
            ind2 = [ind2 (j-1)*mn+i1];
        end
        xx = Z(ind1,ind2);
        Thetastat(i,j) = norm(xx,'fro');
    end
end
thresh = 0;
edgest = (Thetastat > thresh); % group norm > 0, there is an edge in 
% multi-attribute graph
%
end
