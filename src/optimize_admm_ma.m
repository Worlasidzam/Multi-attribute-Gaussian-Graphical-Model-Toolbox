%
% function optimize_admm_ma.m
%
function [Theta, Z, edgest, iter] = optimize_admm_ma(X,pa,p, ...
    maxiter,rho,lambda,alpha,tol,varrho,niteradap,adapflag)
%
m = pa/p;
matr = m;
lamsc = m;
if (m > 1)
    lambda1 = alpha*lambda*(ones(pa,pa)-eye(pa));
    lambda2 = (1-alpha)*lamsc*lambda*(ones(p,p)-eye(p));
    lambda1s = lambda*(ones(pa,pa)-eye(pa));
    lambda2s = lambda*(ones(p,p)-eye(p));
else
    lambda1 = lambda*(ones(p,p)-eye(p));
    lambda2 = 0;
end
%
% optimize using ADMM
%
rho = 2; % ADMM parameter
maxiter = 200; % maximum number of iterations allowed for optimization
tol = 0.0001;
varrho = 1;
[Theta, Z, edgest, iter] = ...
    opt_admm_ma(X,pa,p,maxiter,rho,lambda1,lambda2,tol,varrho);
if ( adapflag ~= 0 )
    count = 0;
    while (count < niteradap)
        count = count +1;
        Thetainit = Theta;
        inflag = 1; %1; initialize with lasso results
        if (count > 1)
            if (alpha > 0)
                lambda1s = lambda1a/alpha;
            end
            if ((1-alpha) > 0)
                lambda2s = lambda2a/((1-alpha)*lamsc);
            end
        end
        if (adapflag == 1) % logsum
            deleps = 0.0001;
            lambda1a = lambda1./(abs(Z)+deleps);
            divz = zeros(p,p);
            for i=1:1:p
                for j=1:1:p
                    if (i ~= j)
                        ind1 = [];
                        ind2 = [];
                        for i1=1:1:matr
                            ind1 = [ind1 (i-1)*matr+i1];
                            ind2 = [ind2 (j-1)*matr+i1];
                        end
                        atemp = Z(ind1, ind2);
                        divz(i,j) = norm(atemp,'fro');
                    end
                end
            end
            lambda2a = lambda2./(divz+deleps);
        elseif (adapflag == 2)
            ap = 3.7;
            % ap = 5;
            part1 = (abs(Z) <= lambda1s).*lambda1s;
            part2 = ((abs(Z) > lambda1s) & (abs(Z) < (ap*lambda1s)) ).* ...
                ((ap*lambda1s - abs(Z))/(ap-1));
            part3 = (abs(Z) >= (ap*lambda1s)).*zeros(pa,pa);
            lambda1a = alpha*(part1 + part2 + part3);
            divz = zeros(p,p);
            for i=1:1:p
                for j=1:1:p
                    if (i ~= j)
                        ind1 = [];
                        ind2 = [];
                        for i1=1:1:matr
                            ind1 = [ind1 (i-1)*matr+i1];
                            ind2 = [ind2 (j-1)*matr+i1];
                        end
                        atemp = Z(ind1, ind2);
                        divz(i,j) = norm(atemp,'fro');
                    end
                end
            end
            part1 = (divz <= lambda2s).*lambda2s;
            part2 = ( (divz > lambda2s) & (divz < (ap*lambda2s)) ).* ...
                ((ap*lambda2s - divz)/(ap-1));
            part3 = (divz >= (ap*lambda2s)).*zeros(p,p);
            lambda2a = (1-alpha)*lamsc*(part1 + part2 + part3);
        end
        %
        rho =2;
        varrho = 1;
        [Theta, Z, edgest, iter] = ...
            opt_admm_ma_adap(X,Thetainit,inflag,pa,p,maxiter,rho,lambda1a,lambda2a,tol,varrho);
    end
end
end