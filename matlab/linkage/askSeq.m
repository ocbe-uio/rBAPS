function [isOK,returnvalue] = askSeq()

items(1).name = 'Yes. Continue loading the sequence profile.';
items(1).default = 1;
items(1).exclusive = [2];
items(1).values = [];

items(2).name = 'No. Stop loading data.';
items(2).default = 0;
items(2).exclusive = [1];
items(2).values = [];

title = 'Load Sequence?';
msg = sprintf(['The allelic profile has been loaded.\nWould you like to continue loading the corresponding gene sequence data?']);
out = CSEFlagDialog(items, title, msg);
if ~(isempty(out)),
    if(out(1).answer==1)
        returnvalue = 1;
    elseif(out(2).answer==1)
        returnvalue = 2;
    end
    isOK = 1;
else
    isOK = 0;
    returnvalue = 0;
end   