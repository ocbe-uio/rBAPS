if ~exist('BAPS_package','dir')
    mkdir('BAPS_package');
end

mcc -m ./general/baps6.m -a ./admixture -a ./general -a ./graph -a ./independent -a ./linkage -a ./parallel -a ./spatial -d ./BAPS_package

