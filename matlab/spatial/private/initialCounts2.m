function [sumcounts, counts] = initialCounts2(partition, data, npops, nLetters)
% initialize counts and sumcounts for the initial partition
%    npops: number of populations in the partition
% nLetters: the maximum number of different symbols over all loci
% Lu Cheng, 25.05.2011

[nSeq nLoci] = size(data);

counts = zeros(nLetters,nLoci,npops);
sumcounts = zeros(nLoci,npops);

for i=1:npops
    inds = (partition==i);
    counts(:,:,i) = histc(data(inds,:),1:nLetters,1);
    sumcounts(:,i) = sum(counts(:,:,i),1);
end

