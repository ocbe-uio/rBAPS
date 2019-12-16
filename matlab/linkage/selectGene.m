function [isOK, genename] = selectGene(genename)
ngenes = size(genename,1);
for i=1:ngenes
    items(i).name = genename{i};
    items(i).default = 1;
    items(i).values = [];
end

title = 'STEP 3';
nGenestr = sprintf('%d',ngenes);
msg = sprintf(['The allelic profile contains ' nGenestr ' genes named below.\n'...
               'Select the individual sequence data that you want to load.\n'...
               'It is recommended that all the genes are selected.']);
out = CSEFlagDialog(items, title, msg);
if ~(isempty(out)),
    for i=1:ngenes
        if ~out(i).answer
        genename{i}=[];
        end
    end
    isOK = 1;
else
    isOK = 0;
    genename = {[]};
end   