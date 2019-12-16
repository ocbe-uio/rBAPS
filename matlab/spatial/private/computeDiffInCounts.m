function [counts sumcounts] = computeDiffInCounts(rows, data, nLetters)
% calculate the counts of the given rows of the data (ninds*nLoci)
% nLetters is the maximum number of different symbols over all loci
% Lu Cheng, 25.05.2011

counts = histc(data(rows,:),1:nLetters,1);
sumcounts = sum(counts,1)';
