% run_cv_only.m
% Cross-validation (5-fold) for selecting lambda in multi-attribute GGM
% Modified on July 8, 2025 to ensure proper reproducibility and randomness
% Initiated on June 30, 2025 to run the Cross validation approach solely

% --- add toolbox source to path ---
root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root,'..','src')));

% --- settings ---
mcrun = 50; % number of Monte Carlo runs
Nv = [200 400 800];
p = 100;
m = 4;
pa = m * p;

model = 4; % 4=BA, 3=ER, 2=Chain (random), 1=Chain (0.2)
fracp = 0.01; % used only for ER (model=3)
print = 1;

adapflag = 1; % 0=lasso, 1=log-sum (LSP reweight), 2=SCAD (if implemented)
niteradap = 1; % number of reweight iterations
rho = 2;
maxiter = 200;
tol = 1e-4;
varrho = 1;
alpha = 0.05;

% --- storage ---
f1cv = zeros(mcrun, length(Nv));
hammingcv = zeros(mcrun, length(Nv));
realfrobcv = zeros(mcrun, length(Nv));

% --- RNG control ---
rng('default');
si = rng;
sl = si;

for ndat = 1:length(Nv)
N = Nv(ndat);
disp(['Sample size: ', num2str(N)]);
fprintf('\n');

for run = 1:mcrun

si = sl;
rng(si);

% Generate graph and data
[Kprec, edgindex, ~] = GenGraphPrec(p, m, fracp, model, print);
K = Kprec;
data = dataGen(Kprec, N, p, m);

sl = rng;

% Full-sample covariance
covmat = (data * data') / N;

% Lambda grid (upper limit via bisection)
[uplim, ~, ~] = bisection_uplim(covmat, N, pa, p, maxiter, rho, alpha, tol, varrho, niteradap, adapflag);
lolim = uplim / 5;
lambdav = logspace(log10(lolim), log10(uplim), 8);

% 5-fold CV over lambda
Kfold = 5;
loss_cv = zeros(length(lambdav), Kfold);
cv = cvpartition(N, "KFold", Kfold);

for i = 1:length(lambdav)
lambda = lambdav(i);

for fold = 1:Kfold
idx_train = training(cv, fold);
idx_test = test(cv, fold);

data_train = data(:, idx_train);
data_test = data(:, idx_test);

cov_train = (data_train * data_train') / sum(idx_train);
cov_test = (data_test * data_test') / sum(idx_test);

[Theta_cv, ~, ~, ~] = optimize_admm_ma(cov_train, pa, p, maxiter, rho, lambda, alpha, tol, varrho, niteradap, adapflag);

% tiny ridge to avoid log(det(.)) issues
Theta_reg = Theta_cv + 1e-6 * eye(pa);

% Negative Gaussian log-likelihood (up to constants)
loss_cv(i, fold) = trace(cov_test * Theta_reg) - log(det(Theta_reg));
end
end

mean_loss = mean(loss_cv, 2);
[~, best_idx_cv] = min(mean_loss);
best_lambda_cv = lambdav(best_idx_cv);

fprintf('Run %d | CV-selected lambda: %.5f\n', run, best_lambda_cv);

% Fit on full data with selected lambda
[Theta_cv_final, ~, edge_cv, ~] = optimize_admm_ma(covmat, pa, p, maxiter, rho, best_lambda_cv, alpha, tol, varrho, niteradap, adapflag);
[f1cv_run, hammingcv_run, errfrob_cv] = performance(edge_cv, edgindex, p, K, Theta_cv_final);

fprintf('CV f1score: %.4f\n\n', f1cv_run);

f1cv(run, ndat) = f1cv_run;
hammingcv(run, ndat) = hammingcv_run;
realfrobcv(run, ndat) = errfrob_cv;
end
end

% --- summary ---
disp('F1score CV:');
f1scorestat_cv = [mean(f1cv, 1); std(f1cv, 1)];
disp(f1scorestat_cv);
fprintf('\n');

disp('Hamming Distance CV:');
hammingdiststat_cv = [mean(hammingcv, 1); std(hammingcv, 1)];
disp(hammingdiststat_cv);