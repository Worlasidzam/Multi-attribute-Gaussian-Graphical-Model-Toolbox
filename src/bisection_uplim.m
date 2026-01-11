function [uplim, count, edmid] = bisection_uplim(X,N,pa,p, ...
    maxiter,rho,alpha,tol,varrho,niteradap,adapflag)
%
if (N <= 400)
    edthresh = round(0.002*p*(p-1)/2); % for N <= 400
else
    edthresh = round(0.005*p*(p-1)/2);  % for N =800
end
bisecflag = 1;
hiflag = 1;
loflag = 1;
bihi = 2;
bilo = 0.1;
%
lambda = bihi;
while (hiflag == 1)
    [Theta, Z, edgest, iter] = optimize_admm_ma(X,pa,p, ...
        maxiter,rho,lambda,alpha,tol,varrho,niteradap,adapflag);
    edhi = sum(sum(edgest));
    if (edhi <= edthresh)
        hiflag = 0;
        bihi = lambda;
        break;
    else
        lambda = lambda*2;
    end
end
%
lambda = bilo;
while (loflag == 1)
    [Theta, Z, edgest, iter] = optimize_admm_ma(X,pa,p, ...
        maxiter,rho,lambda,alpha,tol,varrho,niteradap,adapflag);
    edlo = sum(sum(edgest));
    if (edlo >= (10*edthresh))
        loflag = 0;
        bilo = lambda;
        break;
    else
        lambda = lambda/2;
    end
end
% start bisection mathod
count = 0;
while (bisecflag == 1)
    bimid = (bilo+bihi)/2;
    lambda = bimid;
    count = count + 1;
    [Theta, Z, edgest, iter] = optimize_admm_ma(X,pa,p, ...
        maxiter,rho,lambda,alpha,tol,varrho,niteradap,adapflag);
    edmid = sum(sum(edgest));
    % diff = abs(edmid - edthresh);
    diff = edmid - edthresh;
    if (( (diff <= (edthresh/2)) & (diff > 0) ) | (count > 20))
        bisecflag = 0;
        uplim = lambda;
        break;
    else
        if (edmid > edthresh)
            bilo = bimid;
        else
            bihi = bimid;
        end
    end
end
end