function [f1score, hammingdist, errfrob] = ...
    performance(edgest,edgindex,p,K,Theta)
%
% Objective: calculate performance measure -- F1 score, Hamming distance,
% and the normalized Frobenius norm of inverse covariance estimation error
% 
% Needs functions: none
%
%
% INPUT:
% edgest = estimated edges: 1 at (i,j) if nodes i and j connected, else 0
% edgindex = true edge-set
% p = number of nodes
% K = true precision matrix
% Theta = precision matrix estimate (from ADMM)
% Z = variable splitting of Theta
%
tp =0; % true positive
fp = 0; % false positive
tn = 0; % true negative
fn = 0; % false negative
hammingdist = 0; % hamming distance
for i=1:1:p
    for j=(i+1):1:p
        if (edgest(i,j) == 1)
            ab = 1;
        else
            ab = 0;
        end
        flagtyp1 = (edgindex(i,j) == 0);
        tp = tp + (1-flagtyp1)*ab;
        fp = fp + flagtyp1*ab;
        tn = tn + flagtyp1*(1-ab);
        fn = fn + (1-flagtyp1)*(1-ab);
        hammingdist = hammingdist + ...
            ab*flagtyp1 + (1-flagtyp1)*(1-ab);
    end
end
prec = tp/(tp+fp +eps);
recall = tp/(tp+fn + eps);
f1score = 2*prec*recall/(prec+recall + eps); %F1-score
errfrob = norm(Theta-K,'fro')/norm(K,'fro'); % normalized error Frobenius norm
end

