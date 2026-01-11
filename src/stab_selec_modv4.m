function [edgeopt, seledge0, lambda0, sl2] = stab_selec_modv4(data,N,lambdav,nsamp,pa,p, ...
    maxiter,rho,alpha,tol,varrho,niteradap,adapflag,pthresh, si2)
%
% nsamp = number of samples for stability selection
%
beta = 0.1;
b = floor(10*sqrt(N));
m= pa/p;
matr = m;
lambdanum = length(lambdav);
edgsel = zeros(p,p,nsamp,lambdanum);
seledge = zeros(p,p,lambdanum);
devlam = zeros(1,lambdanum);
devmod = devlam;
rsindexst = zeros(b,nsamp); % stores indices for all random samples
%
% collect all random sample indices first
%
sl2 = si2;
for isamp = 1:1:nsamp
    si2 = sl2;
    rng(si2);
    rsindex = randsample(N,b);
    rsindexst(:,isamp) = rsindex;
    % disp([' first number= ',' ',num2str(rsindex(1))]);
    %
    sl2 = rng;
end
%
% loop over lambda
%
for nlambda = lambdanum:-1:1
    lambda =lambdav(nlambda);
    for isamp = 1:1:nsamp
        rsindex = rsindexst(:,isamp);
        rsindex = sort(rsindex');
        datasamp = data(:,rsindex);
        covmat = (datasamp*datasamp')/b;
        X = covmat;
        [Theta, Z, edgest, iter] = optimize_admm_ma(X,pa,p, ...
            maxiter,rho,lambda,alpha,tol,varrho,niteradap,adapflag);
        edgsel(:,:,isamp,nlambda) = edgest;
    end
    seledge = squeeze(mean(edgsel,3)); % average over samples
    bernind = 2.*seledge.*(1-seledge);
    sum = 0;
    for i=1:1:p
        for j=i+1:1:p
            sum = sum + bernind(i,j,nlambda);
        end
    end
    devlam(nlambda) = 4*sum/(p*(p-1));
    devmod(nlambda) = max(devlam(nlambda:end));
    if (devmod(nlambda) > beta)
        lambda0 = lambdav(nlambda+1);
        break;
    end
end
devmod
%
% bisection search to refine lambda selection
%
bihi = lambda0;
bilo = lambdav(nlambda);
devhi = devmod(nlambda+1);
devlo = devmod(nlambda);
edgselbi = zeros(p,p,nsamp);
seledgebi = zeros(p,p);
bisecflag = 1;
% start bisection mathod
count = 0;
while (bisecflag == 1)
    bimid = (bilo+bihi)/2;
    lambda = bimid;
    count = count + 1;
    for isamp = 1:1:nsamp
        rsindex = rsindexst(:,isamp);
        rsindex = sort(rsindex');
        datasamp = data(:,rsindex);
        covmat = (datasamp*datasamp')/b;
        X = covmat;
        [Theta, Z, edgest, iter] = optimize_admm_ma(X,pa,p, ...
            maxiter,rho,lambda,alpha,tol,varrho,niteradap,adapflag);
        edgselbi(:,:,isamp) = edgest;
    end
    seledgebi = squeeze(mean(edgselbi,3)); % average over samples
    bernindbi = 2.*seledgebi.*(1-seledgebi);
    sum = 0;
    for i=1:1:p
        for j=i+1:1:p
            sum = sum + bernindbi(i,j);
        end
    end
    devmid = 4*sum/(p*(p-1));
    % diff = abs(edmid - edthresh);
    diff = devmid - beta;
    if (abs(diff) <= 0.01 | count > 6)
        bisecflag = 0;
        lambda0 = lambda;
        break;
    else
        if (devmid > beta)
            bilo = bimid;
        else
            bihi = bimid;
        end
    end
end
disp([' selected lambda= ',' ',num2str(lambda0),' ',' devmid = ',' ', ...
            num2str(devmid),' ',' count (no of iter) = ',num2str(count)]);
%
%
%indx0 = min(find(lambdav(:) == lambda0));
seledge0 = seledgebi; % squeeze(seledge(:,:,indx0));
edgeopt = (seledgebi >= pthresh); % (seledge(:,:,indx0) >= pthresh);
%
end