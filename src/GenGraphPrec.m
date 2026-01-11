%
% Generate the precision matrix for the specified model
%
function [Kprec, edgindex, spfac] = ... 
    GenGraphPrec(p,m,fracp,model,print)
%
% modified on Dec 29, 2024 to caompute condition number
%
% Objective: generate mp*mp precision matrix K 
%
%
% Needs functions: BAmodel_mod.m for the BA model
% %
% set parameters
% model =1 -> chain graph: nonzero off-diagonals are 0.2
% model =2 -> chain graph: nonzero off-diagonals are random over +-[0.1,0.4]
% model =3 -> Erdos-renyi with prob=frac
% model =4  -> BA model
%
% option = 1;
mineig = 0.5; %0.5; % 0.05 minimum eigenvalue of precsion matrix K
eigenflag = 1; % if 0, don't manipulate min eigenvalue;
%
% fracp: connection probability in Erdos-Renyi graph
%
%
% adjacency matrix adj(p,p) of Graph X: which nodes are connected
%
adj = zeros(p,p);
if ( (model == 1) | (model == 2) ) % both are chain graphs
    for i=1:1:p
        if (i+1 <= p)
            adj(i,i+1) = 1; % connect to immediate neighbor
            adj(i+1,i) = 1;  % connect to immediate neighbor
        end
    end
elseif (model == 3) % Erdos-Renyi graph
    asign = unifrnd(0,1,p,p);
    adj = (asign <= fracp); % 100*frac percent connected
    %symmetrisize
    Bu = triu(adj);
    Bl = Bu';
    Bl(logical(eye(p))) = zeros(p,1);
    Bu(logical(eye(p))) = ones(p,1);
    adj = Bu+Bl;
elseif (model == 4)
    m0 = 5; % mean degree
    adj = BAModel_mod(p,m0); % diagonals are zeros
    adj = adj +eye(p);
end
pa = m*p;
matr = m;
%
% precision matrix K(pa,pa)
%
K = zeros(pa,pa);
for i=1:1:p
    for j=1:1:p
        for s=1:1:m
            for t=1:1:m
                row = (i-1)*m+s;
                col = (j-1)*m+t;
                if (i == j) % diagonal block
                    K(row,col) = 0.5^(abs(s-t));
                elseif (adj(i,j) == 1)
                    if (m > 1) % multi-attribute
                        if (s ~= t)
                            if (model == 1)
                                K(row,col) = 0.2;
                            elseif ((model == 2) | (model == 3) | (model == 4))
                                abc = rand(1,1);
                                abd = (abc > 0.5) - (abc <= 0.5);
                                K(row,col) = (0.1+0.3*rand(1,1))*abd;
                            end
                        else
                            K(row,col) = 0;
                        end
                    else % single attribute
                        if (model == 1)
                            K(row,col) = 0.2;
                        elseif ((model == 2) | (model == 3) | (model == 4))
                            abc = rand(1,1);
                            abd = (abc > 0.5) - (abc <= 0.5);
                            K(row,col) = (0.1+0.3*rand(1,1))*abd;
                        end
                    end
                end
            end
        end
    end
end
% pick upper triangle and make it symmetric
B = triu(K);
K=B+B' - diag(diag(B));
% "stabilize"
rho = min(eig(K));
if (eigenflag == 1)
    K= K + (mineig-rho)*eye(pa); % minimum eigenvalue is mineig
    condnum = (max(eig(K))+(mineig-rho))/mineig;
else
    condnum = max(eig(K))/rho;
end
Kprec = K;
%
%
%
% determine connected edges: CHECK
for i=1:1:p
    for j=1:1:p
        flag = 0;
        for s=1:1:m
            for t=1:1:m
                if (abs(K(s+(i-1)*m,t+(j-1)*m)) > 0)
                    flag = 1;
                end
            end
        end
        edgindex(i,j) = flag;
    end
end
%
count = sum(sum(edgindex));
count = (count-p)/2; % count only off-diagonal elements
spfac = count/(p*(p-1)/2); % fraction of connected off-diagonal elements
%
%
if (print == 1)
    disp([' spfac= ',' ',num2str(spfac),' ',' condnum= ',' ',num2str(condnum)]);
end
%
%
end
%
