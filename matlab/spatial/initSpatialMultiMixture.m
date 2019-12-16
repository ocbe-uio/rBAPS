function [partition, counts, sumcounts] = initSpatialMultiMixture(initData, ... 
    npops, Z, rows, noalle, dist, adjprior, priorTerm, fixedK);
% Etsii spatial multimixturelle alkutilan baps 3.1:n ahneella algoritmilla.
% toimii!

global PARTITION; global COUNTS;
global SUMCOUNTS; global POP_LOGML;

c.data = initData; c.Z = Z; c.rows=rows; c.rowsFromInd=0; c.noalle=noalle;
c.dist = dist; c.adjprior = adjprior; c.priorTerm = priorTerm;

indMix_fixK(c,npops,1,0);

partition = PARTITION; counts = COUNTS; sumcounts = SUMCOUNTS;

  


