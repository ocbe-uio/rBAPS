function [cliqcounts, sepcounts] = computeCounts(cliques, separators, npops)

global PARTITION;
ncliq = size(cliques,1);
nsep = size(separators,1);

cliqPartition = zeros(ncliq, size(cliques,2));
sepPartition = zeros(nsep, size(separators, 2));

apuCliq = find(cliques > 0);
apuSep = find(separators > 0);

cliqPartition(apuCliq) = PARTITION(cliques(apuCliq));
sepPartition(apuSep) = PARTITION(separators(apuSep));


cliqcounts = zeros(ncliq, npops);
for i = 1:npops
    cliqcounts(:,i) = sum(cliqPartition == i, 2);
end
    
sepcounts = zeros(nsep, npops);
for i = 1:npops
    sepcounts(:,i) = sum(sepPartition == i, 2);
end

%-------------------------------------------------------------------------