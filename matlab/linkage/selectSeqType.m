function nextstep = selectSeqType()

items(1).name = 'SEQUENCE TYPE:';
items(1).default = 1;
items(1).linked = [2 3 4];
items(1).values = [];

items(2).name = '1 = Non-coding nucleotide';
items(2).default = 0;
items(2).exclusive = [3 4];
items(2).indent = 1;

items(3).name = '2 = Coding nucleotide';
items(3).default = 1;
items(3).exclusive = [2 4];
items(3).indent = 1;

items(4).name = '3 = Protein';
items(4).default = 0;
items(4).exclusive = [2 3];
items(4).indent = 1;

items(5).name = 'GENETIC CODE:  ';
items(5).default = 1;
items(5).indent = 0;
items(5).values = {'1 - Standard';'2 - Vertebrate Mithchondrial'; ...
                   '3 - Yeast Mithchondrial';'4 - Mold Mithchondrial'; ...
		   '5 - Invertebrate Mithchondrial';'6 - Protozoan Mithchondrial'; ...
                   '7 - Coelenterate Mithchondrial';'8 - Mycoplasma'};
items(5).help = 'Please select a genetic code';

title = 'Sequence type and genetic code';
% msg = sprintf(['Please select sequence type and genetic code']);
out = CSEFlagDialog(items, title);
if ~(isempty(out)),
    if(out(2).answer==1)
        seqtype=1;
    elseif(out(3).answer==1)
        seqtype=2;
    elseif(out(4).answer==1)
        seqtype=3;
    end
    geneticcode=out(5).answer;
else
    seqtype=[];
    geneticcode=[];
end   