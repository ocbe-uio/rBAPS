function [isOK, datatype] = selectDataType()

items(1).name = 'MLST DATA';
items(1).default = 1;
items(1).linked = [2 3];
items(1).exclusive = [4 5 6];
items(1).values = {1};

items(2).name = '1 = Allelic profile';
items(2).default = 1;
items(2).exclusive = [3 4 5 6];
items(2).indent = 1;

items(3).name = '2 = FASTA-format';
items(3).default = 0;
items(3).exclusive = [2 4 5 6];
items(3).indent = 1;


title = 'STEP 1';
msg = sprintf(['Please specify the data format:']);
out = CSEFlagDialog(items, title, msg);
if ~(isempty(out)),
    if(out(2).answer==1)
        datatype = 1;
    elseif(out(3).answer==1)
        datatype = 2;
    elseif(out(5).answer==1)
        datatype = 3;
    elseif(out(6).answer==1)
        datatype = 4;
    end
    isOK = 1;
else
    isOK = 0;
    datatype = 0;
end   