function popnames = initPopNames(nameFile, indexFile)
%Palauttaa tyhjän, mikäli nimitiedosto ja indeksitiedosto
% eivät olleet yhtä pitkiä.

popnames = [];
indices = load(indexFile);

fid = fopen(nameFile);
if fid == -1
    %File didn't exist
    msgbox('Loading of the population names was unsuccessful', ...
        'Error', 'error');
    return;
end;
line = fgetl(fid);
counter = 1;
while (line ~= -1) && ~isempty(line)
    names{counter} = line;
    line = fgetl(fid);
    counter = counter + 1;
end;
fclose(fid);

if length(names) ~= length(indices)
    disp('The number of population names must be equal to the number of ');
    disp('entries in the file specifying indices of the first individuals of ');
    disp('each population.');
    return;
end

popnames = cell(length(names), 2);
for i = 1:length(names)
    popnames{i,1} = names(i);
    popnames{i,2} = indices(i);
end
