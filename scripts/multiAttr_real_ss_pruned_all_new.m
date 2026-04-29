% multiAttr_real_ss_unpruned_only.m
%
% Stability Selection version for real data
% SHOWS: UNPRUNED weighted adjacency in Figure 520
%

p = 97; % number of nodes (stocks)
m = 4; % number of attributes per node
matr0 = m;
load logMAdata97modv2;
data = logMAdata97modv2;
classind = clas;
pa = p*m;
[pa, N] = size(data); % N is the original sample size
size(data)
attrflag = 1; % if 1 MA, if 0 SA
sampchange = 0; % 1 if you want to use a smaller sample size than N, else 0
if (sampchange == 1)
    % sample size to be used
    N0 = 500;
    Ninit = 100;
    Nfinal = Ninit + N0 -1;
    data = data(:,Ninit:Nfinal);
    size(data)
    N = N0;
end

if (attrflag == 0)
    m=1; % single attribute
    matr = m;
    data1 = zeros(p,N);
    for in =1:1:p
        indd = (in-1)*matr0 + 3;
        data1(in,:) = data(indd,:);
    end
else
    matr = matr0;
    data1 = data;
end

maxiter0 = 500; % maximum number of iterations allowed for optimization
maxiter = maxiter0;
tol0 = 0.0001;
tol = tol0;
rho = 2; % ADMM parameter initial setting
varrho = 1; % varibale ADMM parameter rho

alpha = 0.05; 
alpha0 = alpha;

adapflag = 1; % 0 for lasso, 1 for LSP, 2 for SCAD
niteradap = 1;
print = 1;

% Set method name for titles
if adapflag == 0
    method_name = 'Lasso';
elseif adapflag == 1
    method_name = 'LSP';
elseif adapflag == 2
    method_name = 'SCAD';
end

% Stability selection parameters
nsamp = 20; % number of subsamples
pthresh = 0.6; % pruning threshold
beta = 0.1; % stability parameter
b = floor(10*sqrt(N)); % subsample size

% Compute covariance matrix
covmat = (data1*data1')/N;
pa = m*p;
X = covmat;

% Find upper limit on lambda using bisection
[uplim, count, edmid] = bisection_uplim(X,N,pa,p, ...
    maxiter,rho,alpha,tol,varrho,niteradap,adapflag);
if (print == 1)
    disp([' uplim= ',' ',num2str(uplim),' ',' count= ',' ',num2str(count), ...
        ' ',' edmid= ',' ',num2str(edmid)]);
end

% For stability selection: uplim to uplim/5 (no halving)
lolim = uplim/5;
lambdav = logspace(log10(lolim), log10(uplim), 8);
if (print == 1)
    lambdav
end

% Run stability selection
fprintf('Running stability selection...\n');
rng('default');
si2 = rng;
[edgeopt, seledge0, lambda0, sl2] = stab_selec_modv4(data1,N,lambdav,nsamp,pa,p, ...
    maxiter,rho,alpha,tol,varrho,niteradap,adapflag,pthresh,si2);

% Ensure edgeopt is symmetric for pruning
edgeopt = edgeopt | edgeopt';
edgeopt = double(edgeopt);

disp([' Stability selected lambda = ',' ',num2str(lambda0)]);
disp([' Pruning threshold = ',num2str(pthresh)]);

% Final run with selected lambda on full data
[Theta, Z, edgest, iter] = optimize_admm_ma(covmat,pa,p, ...
    maxiter,rho,lambda0,alpha,tol,varrho,niteradap,adapflag);

nume_unpruned = sum(sum(edgest));
nume_pruned = sum(sum(edgeopt))/2; % Divide by 2 because symmetric
disp([' number of edges (unpruned) = ',' ',num2str(nume_unpruned)]);
disp([' number of edges (pruned) = ',' ',num2str(nume_pruned)]);

% Create weighted adjacency matrix (full unpruned)
Thetastat = zeros(p,p);
mn = matr;
for i=1:1:p
    for j=1:1:p
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
ZZ = Thetastat; % full unpruned matrix

% Create pruned weighted adjacency (for other figures)
Thetastat_pruned = Thetastat .* edgeopt;
Thetastat_pruned = (Thetastat_pruned + Thetastat_pruned')/2;

% Figure 520: UNPRUNED Weighted Adjacency
figure(520)
weights = Thetastat; % Use unpruned weights
imagesc(weights)
title(['weighted adjacency: SS (unpruned), ', method_name]);

% colormap logic
newmap = colormap(linspecer);
minz = double(min(min(weights)));
if minz > 0
    disp('Your range is all above 0, no change');
else
    maxz = double(max(max(weights)));
    if maxz < 0
        disp('Your range is all below 0, no change');
    else
        ncol = size(newmap, 1);
        zratio = (0 - minz) ./ (maxz - minz);
        zpos = max( round(zratio * ncol), 1);
        newmap(zpos,:) = [1 1 1];
        colormap(newmap);
    end
end
hold on;

% Add sector boundaries (same as BIC version)
line([0.5,12.5],[12.5,12.5],'Color','r','LineStyle','--');
line([12.5,12.5],[0.5,12.5],'Color','r','LineStyle','--');
drawblock(12.5,27.5);
drawblock(27.5,44.5);
drawblock(44.5,46.5);
line([46.5,56.5],[46.5,46.5],'Color','r','LineStyle','--');
line([46.5,46.5],[46.5,56.5],'Color','r','LineStyle','--');
line([46.5,56.5],[56.5,56.5],'Color','r','LineStyle','--');
line([56.5,56.5],[46.5,56.5],'Color','r','LineStyle','--');
line([56.5,68.5],[56.5,56.5],'Color','r','LineStyle','--');
line([56.5,56.5],[56.5,68.5],'Color','r','LineStyle','--');
line([56.5,68.5],[68.5,68.5],'Color','r','LineStyle','--');
line([68.5,68.5],[56.5,68.5],'Color','r','LineStyle','--');
line([68.5,76.5],[68.5,68.5],'Color','r','LineStyle','--');
line([68.5,68.5],[68.5,76.5],'Color','r','LineStyle','--');
line([68.5,76.5],[76.5,76.5],'Color','r','LineStyle','--');
line([76.5,76.5],[68.5,76.5],'Color','r','LineStyle','--');
line([76.5,87.5],[76.5,76.5],'Color','r','LineStyle','--');
line([76.5,76.5],[76.5,87.5],'Color','r','LineStyle','--');
line([76.5,87.5],[87.5,87.5],'Color','r','LineStyle','--');
line([87.5,87.5],[76.5,87.5],'Color','r','LineStyle','--');
line([87.5,92.5],[87.5,87.5],'Color','r','LineStyle','--');
line([87.5,87.5],[87.5,92.5],'Color','r','LineStyle','--');
line([87.5,92.5],[92.5,92.5],'Color','r','LineStyle','--');
line([92.5,92.5],[87.5,92.5],'Color','r','LineStyle','--');
line([92.5,93.5],[92.5,92.5],'Color','r','LineStyle','--');
line([92.5,92.5],[92.5,93.5],'Color','r','LineStyle','--');
line([92.5,93.5],[93.5,93.5],'Color','r','LineStyle','--');
line([93.5,93.5],[92.5,93.5],'Color','r','LineStyle','--');
line([93.5,97.5],[93.5,93.5],'Color','r','LineStyle','--');
line([93.5,93.5],[93.5,97.5],'Color','r','LineStyle','--');

colorbar

% SAVE FIGURE 520
filename_fig520 = ['SS_',method_name,'_unpruned_',num2str(nume_unpruned),'edges.png'];
exportgraphics(gcf, filename_fig520, 'Resolution', 300);
fprintf('Saved %s\n', filename_fig520);

% Figure 521: Selection probabilities
figure(521)
imagesc(seledge0)
caxis([0 1])
title('Selection Probabilities');
colormap(flipud(gray));
colorbar

% SAVE FIGURE 521
filename_fig521 = ['SS_',method_name,'_selection_probs.png'];
exportgraphics(gcf, filename_fig521, 'Resolution', 300);
fprintf('Saved %s\n', filename_fig521);

% Figure 199: Force-directed layout (unpruned)
ZZ2 = Thetastat - diag(diag(Thetastat));
ZZ2 = 0.5*(ZZ2 + ZZ2');
G2 = graph(ZZ2);
figure(199)
p0=plot(G2,'NodeColor','r','Layout','force', ...
    'Iterations',50,'UseGravity',true);
p0.MarkerSize =2;
title('estimated graph (unpruned)');
axis off;

% SAVE FIGURE 199
filename_fig199 = ['SS_',method_name,'_force_unpruned.png'];
exportgraphics(gcf, filename_fig199, 'Resolution', 300);
fprintf('Saved %s\n', filename_fig199);

% Figure 198: Circle layout (unpruned)
figure(198)
p0=plot(G2,'NodeColor','r','Layout','circle');
p0.MarkerSize =2;
p0.NodeLabel = {};
title('estimated graph (unpruned)');
axis off;

% SAVE FIGURE 198
filename_fig198 = ['SS_',method_name,'_circle_unpruned.png'];
exportgraphics(gcf, filename_fig198, 'Resolution', 300);
fprintf('Saved %s\n', filename_fig198);

% Figure 197: Pruned binary graph
A = edgeopt - diag(diag(edgeopt));
A = 0.5*(A + A');
figure(197)
Gp = graph(A);
hp = plot(Gp,'NodeColor','r','Layout','force','Iterations',50,'UseGravity',true);
hp.MarkerSize = 2;
hp.NodeLabel = {};
title(['pruned graph (',num2str(nume_pruned),' edges)']);
axis off;

% SAVE FIGURE 197
filename_fig197 = ['SS_',method_name,'_pruned_binary.png'];
exportgraphics(gcf, filename_fig197, 'Resolution', 300);
fprintf('Saved %s\n', filename_fig197);

% Helper function for drawing sector blocks
function drawblock(x1, x2)
    line([x1,x2],[x1,x1],'Color','r','LineStyle','--');
    line([x1,x1],[x1,x2],'Color','r','LineStyle','--');
    line([x1,x2],[x2,x2],'Color','r','LineStyle','--');
    line([x2,x2],[x1,x2],'Color','r','LineStyle','--');
end                                                                                                                                                                                