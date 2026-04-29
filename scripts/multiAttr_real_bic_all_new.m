% multiAttr_real_bic_all_new.m
%
% Updated multiAttr_real_bic_all.m on Feb 6, 2026: use bisection etc
%
% main program for real data to use the ADMM algorithm with BIC for penalty parameter
% selection
%
% Initiated Dec 24, 2024,  based on transferDrive\graphical_new\DTimeSeries
%  real_dtsFlat_varlambdaBIC.m
%
% Functions called: opt_admm_ma.m,
%                   opt_admm_ma_adap.m
%
%
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
%
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
%
maxiter0 = 500; % maximum number of iterations allowed for optimization
maxiter = maxiter0;
tol0 = 0.0001;
tol = tol0;
rho = 2; % ADMM parameter initial setting
varrho = 1; % varibale ADMM parameter rho
%
%
%
alpha = 0.05; 
alpha0 = alpha;
%
adapflag = 0; % 1 if you want to run LSP also, 2 if SCAD , 0 for lasso
niteradap = 1; %2;
print = 1;
printbic = 1;
%
%
covmat = (data1*data1')/N;
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
        ' ',' edmid= ',' ',num2str(edmid)]);
end
% the following is for BIC
uplim = uplim/2;
lolim = uplim/10;
lambdav = logspace(log10(lolim),log10(uplim),8);
if (print == 1)
    lambdav
end
[Theta, Z, edgest, lambda] = bic_selec(X,N,lambdav,pa,p, ...
    maxiter,rho,alpha,tol,varrho,niteradap,adapflag,printbic);
%
nume = sum(sum(edgest));
disp([' number of edges = ',' ',num2str(nume)]);
%
Thetastat = zeros(p,p); % edge detection
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
ZZ = Thetastat; % - diag(diag(Thetastat)); keep the diagonal for plotting
%
figure(520)
weights = ZZ;
imagesc(weights)
title(' weighted adjacency ');
% newmap = hicontrast();       %change as desired, e.g., flag(256)
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
        zpos = max( round(zratio * ncol), 1);  %closest non-zero
        newmap(zpos,:) = [1 1 1];   %set there to white
        colormap(newmap);           %activate it
    end
end
hold on;
% Information technology
% $ drawblock(0.5,12.5);
line([0.5,12.5],[12.5,12.5],'Color','r','LineStyle','--'); % Information technology
line([12.5,12.5],[0.5,12.5],'Color','r','LineStyle','--');
% Health Care
drawblock(12.5,27.5);
% line([12.5,27.5],[12.5,12.5],'Color','r','LineStyle','--'); % Health Care
% line([12.5,12.5],[12.5,27.5],'Color','r','LineStyle','--');
% line([12.5,27.5],[27.5,27.5],'Color','r','LineStyle','--'); % Health Care
% line([27.5,27.5],[12.5,27.5],'Color','r','LineStyle','--');
% Financials
drawblock(27.5,44.5);
% line([27.5,44.5],[27.5,27.5],'Color','r','LineStyle','--'); % Financials
% line([27.5,27.5],[27.5,44.5],'Color','r','LineStyle','--');
% line([27.5,44.5],[44.5,44.5],'Color','r','LineStyle','--'); % Financials
% line([44.5,44.5],[27.5,44.5],'Color','r','LineStyle','--');
% Real Estate
drawblock(44.5,46.5);
% line([44.5,46.5],[44.5,44.5],'Color','r','LineStyle','--'); %
% line([44.5,44.5],[44.5,46.5],'Color','r','LineStyle','--');
% line([44.5,46.5],[46.5,46.5],'Color','r','LineStyle','--'); %
% line([46.5,46.5],[44.5,46.5],'Color','r','LineStyle','--');
% Consumer Discretionary
line([46.5,56.5],[46.5,46.5],'Color','r','LineStyle','--'); %
line([46.5,46.5],[46.5,56.5],'Color','r','LineStyle','--');
line([46.5,56.5],[56.5,56.5],'Color','r','LineStyle','--'); %
line([56.5,56.5],[46.5,56.5],'Color','r','LineStyle','--');
% Industrials
line([56.5,68.5],[56.5,56.5],'Color','r','LineStyle','--'); %
line([56.5,56.5],[56.5,68.5],'Color','r','LineStyle','--');
line([56.5,68.5],[68.5,68.5],'Color','r','LineStyle','--'); %
line([68.5,68.5],[56.5,68.5],'Color','r','LineStyle','--');
% Communications Services
line([68.5,76.5],[68.5,68.5],'Color','r','LineStyle','--'); %
line([68.5,68.5],[68.5,76.5],'Color','r','LineStyle','--');
line([68.5,76.5],[76.5,76.5],'Color','r','LineStyle','--'); %
line([76.5,76.5],[68.5,76.5],'Color','r','LineStyle','--');
% Consumer Staples
line([76.5,87.5],[76.5,76.5],'Color','r','LineStyle','--'); %
line([76.5,76.5],[76.5,87.5],'Color','r','LineStyle','--');
line([76.5,87.5],[87.5,87.5],'Color','r','LineStyle','--'); %
line([87.5,87.5],[76.5,87.5],'Color','r','LineStyle','--');
% Energy
line([87.5,92.5],[87.5,87.5],'Color','r','LineStyle','--'); %
line([87.5,87.5],[87.5,92.5],'Color','r','LineStyle','--');
line([87.5,92.5],[92.5,92.5],'Color','r','LineStyle','--'); %
line([92.5,92.5],[87.5,92.5],'Color','r','LineStyle','--');
% Materials
line([92.5,93.5],[92.5,92.5],'Color','r','LineStyle','--'); %
line([92.5,92.5],[92.5,93.5],'Color','r','LineStyle','--');
line([92.5,93.5],[93.5,93.5],'Color','r','LineStyle','--'); %
line([93.5,93.5],[92.5,93.5],'Color','r','LineStyle','--');
% Utilities
line([93.5,97.5],[93.5,93.5],'Color','r','LineStyle','--'); %
line([93.5,93.5],[93.5,97.5],'Color','r','LineStyle','--');
% line([93.5,97.5],[97.5,97.5],'Color','r','LineStyle','--'); %
% line([97.5,97.5],[93.5,97.5],'Color','r','LineStyle','--');
%
colorbar
filename_fig520 = ['BIC_', 'LSP_', num2str(nume), 'edges.png'];  % Change 'LSP' to 'Lasso' when running for Lasso
exportgraphics(gcf, filename_fig520, 'Resolution', 300);
fprintf('Saved %s\n', filename_fig520);
%
% added Jan 5, 2026
ZZ2 = Thetastat - diag(diag(Thetastat));
ZZ2 = 0.5*(ZZ2 + ZZ2');
G2 = graph(ZZ2);
figure(199)
p0=plot(G2,'NodeColor','r','Layout','force', ...
    'Iterations',50,'UseGravity',true);
p0.MarkerSize =2;
% p0.NodeLabel = {};
title('estimated graph');
axis off;
%
ZZ2 = Thetastat - diag(diag(Thetastat));
ZZ2 = 0.5*(ZZ2 + ZZ2');
% ZZ2 = (ZZ2>0);
G2 = graph(ZZ2);
figure(198)
p0=plot(G2,'NodeColor','r','Layout','circle');
p0.MarkerSize =2;
p0.NodeLabel = {};
title('estimated graph');
axis off;
%