% run_stability_bic.m
%
% Initiated May 17, 2025, modified Jun 22 to include BIC based approach
% modified May 14, 2025  -- choose uplim via bisection
%
% main program to use the ADMM algorithm with BIC for penalty parameter
% selection
%
% Functions called: GenGraphPrec.m, dataGen.m, opt_admm_ma.m,
%                   performance.m, opt_admm_ma_adap.m

% Add project paths so MATLAB can find all functions
thisFile = mfilename('fullpath');
thisDir  = fileparts(thisFile);
rootDir  = fileparts(thisDir);

addpath(genpath(fullfile(rootDir, 'src')));
%
%
%
% Example settings :
% User set parameters for optimization
%
mcrun = 2; % number of Monte Carlo runs
nsamp = 20; % number of subsamples
pthresh = 0.6; % pruning threshold
Nv = [200 400 800];
p = 100;
m = 4; % number of attributes per node
matr = m;
%
bicflag = 1; % 1 if want to run BIC method too, else 0
adapflag = 0; % 1 if you want to run LSP also, 2 if SCAD , 0 for lasso
niteradap = 1; %2;
model = 2;  % 4 is BA graph, 3 is Erdos-Renyi graph with connection probability fracp
fracp = 0.01; % needed for ER graph
print = 1;
printbic =0;
%
%
rho = 2; % ADMM parameter initial setting
maxiter = 200; % maximum number of iterations allowed for optimization
tol = 0.0001;
varrho = 1; % varibale ADMM parameter rho
%
%
%
alpha = 0.05;
alpha0 = alpha;
%
%
% processing starts here
%
%
datnum = length(Nv);
% the following variables will store performance meaures: initialization
hamming = zeros(mcrun,datnum); % Hamming distance
f1 = zeros(mcrun,datnum); % F1 score
hammingp = hamming;
f1p = f1;
relfrob = zeros(mcrun,datnum); %the normalized Frobenius norm of
%inverse covariance estimation error
timing = zeros(mcrun,datnum);
iterations = zeros(mcrun,datnum);
if (bicflag == 1)
    hammingbic = hamming;
    f1bic = f1;
    relfrobbic = f1;
end
%
%
for ndat = 1:1:datnum
    N = Nv(ndat);
    N
    %
    rng('default'); % rng(0,"twister")
    si = rng;
    sl = si;
    %
    rng(1); % for stability selection
    si2 = rng;
    sl2 = si2;
    %
    tStart = tic;
    Tim = zeros(1,mcrun);
    %
    for run=1:1:mcrun
        if (rem(run,20) == 0)
            disp([' run= ',' ',num2str(run),' ',' N= ',' ',num2str(N)]);
            % run
        end
        si = sl;
        rng(si);
        [Kprec, edgindex, spfac] = ...
            GenGraphPrec(p,m,fracp,model,print);
        K = Kprec; % chosen precision matrix
        [data] = dataGen(Kprec,N,p,m);
        %
        % end of data generation
        %
        sl = rng;
        %
        % sample covariance matrix
        %
        covmat = (data*data')/N;
        pa = m*p;
        %
        %
        lamsc = m;
        X = covmat;
        %
        %
        %  try figuring out the upper limit on lambda
        %
        [uplim, count, edmid] = bisection_uplim(X,N,pa,p, ...
            maxiter,rho,alpha,tol,varrho,niteradap,adapflag);
        if (print == 1)
            disp([' uplim= ',' ',num2str(uplim),' ',' count= ',' ',num2str(count), ...
                ' ',' edmid= ',' ',num2str(edmid),' ',' run= ',' ',num2str(run), ...
                ' ',' sample size = ',' ',num2str(N)]);
        end
        %
        % uplim = uplim/2;
        uplim = uplim;
        % lolim = uplim/10;
        lolim = uplim/5;
        % lambdav = logspace(log10(lolim),log10(uplim),8);
        lambdav = logspace(log10(lolim),log10(uplim),8);
        if (print == 1)
            lambdav
        end
        %
        si2 = sl2;
        [edgeopt, seledge0, lambda0, sl2] = stab_selec_modv4(data,N,lambdav,nsamp,pa,p, ...
            maxiter,rho,alpha,tol,varrho,niteradap,adapflag,pthresh,si2);
        %
        %
        X = covmat;
        %
        % optimize using ADMM
        rho = 2; % ADMM parameter
        maxiter = 200; % maximum number of iterations allowed for optimization
        tol = 0.0001;
        varrho = 1;
        lambda = lambda0;
        tic;
        [Theta, Z, edgest, iter] = optimize_admm_ma(X,pa,p, ...
            maxiter,rho,lambda,alpha,tol,varrho,niteradap,adapflag);
        sum(sum(sum(abs(Z) > 0)));
        % calculate performance measures
        Tim(run) = toc;
        timing(run,ndat) = toc;
        iterations(run,ndat) = iter;
        [f1score, hammingdist, errfrob] = performance(edgest,edgindex,p,K,Theta);
        [f1scorep, hammingdistp, errfrob] = performance(edgeopt,edgindex,p,K,Theta);
        disp([' f1score= ',' ',num2str(f1score),' ',' (pruned) f1scorep= ',' ', ...
            num2str(f1scorep),' ',' pruning threshold = ',num2str(pthresh)]);
        hamming(run,ndat) = hammingdist;
        f1(run,ndat) = f1score;
        hammingp(run,ndat) = hammingdistp;
        f1p(run,ndat) = f1scorep;
        relfrob(run,ndat) = errfrob;
        %
        if (bicflag == 1) % implement BIC method
            uplim = uplim/2; % /2;
            lolim = uplim/10;
            lambdav = logspace(log10(lolim),log10(uplim),8);
            if (print == 1)
                lambdav
            end
            [Theta, Z, edgest, lambda] = bic_selec(X,N,lambdav,pa,p, ...
                maxiter,rho,alpha,tol,varrho,niteradap,adapflag,printbic);
            [f1score, hammingdist, errfrob] = performance(edgest,edgindex,p,K,Theta);
            disp([' BIC f1score= ',' ',num2str(f1score),' ', ...
                'selected lambda= ',' ',num2str(lambda)]);
            hammingbic(run,ndat) = hammingdist;
            f1bic(run,ndat) = f1score;
            relfrobbic(run,ndat) = errfrob;
        end
    end
    tcum = sum(Tim);
    disp([' time per run ',' ',num2str(tcum/mcrun)])
end
%
tEnd = toc(tStart)
%disp([' time per run ',' ',num2str(tcum/mcrun)])
%
display(' results for STABILITY without pruning')
f1scorestat = [mean(f1,1) ; std(f1,1)] % calculate mean and standard deviation
hammingdiststat = [mean(hamming,1) ; std(hamming,1)]
display(' results for STABILITY with pruning')
f1scorestatp = [mean(f1p,1) ; std(f1p,1)] % calculate mean and standard deviation
hammingdiststatp = [mean(hammingp,1) ; std(hammingp,1)]
% errfrobstat = [mean(relfrob,1) ; std(relfrob,1)]
% aveiter = [mean(iterations,1) ; std(iterations,1)]
% avetime = [mean(timing) ; std(timing)]
display(' results using BIC')
f1scorestatbic = [mean(f1bic,1) ; std(f1bic,1)] % calculate mean and standard deviation
hammingdiststatbic = [mean(hammingbic,1) ; std(hammingbic,1)]