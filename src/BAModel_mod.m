function adj = BAModel_mod(p, m0)
% Barabási-Albert model graph generation
% p = number of nodes
% m0 = number of initial nodes to connect

adj = zeros(p);
% Start with a complete graph of m0 nodes
for i = 1:m0
    for j = i+1:m0
        adj(i,j) = 1;
        adj(j,i) = 1;
    end
end

% Add remaining nodes with preferential attachment
for i = m0+1:p
    degrees = sum(adj(1:i-1, 1:i-1), 2);
    prob = degrees / sum(degrees);
    targets = randsample(i-1, m0, true, prob);
    adj(i, targets) = 1;
    adj(targets, i) = 1;
end
end