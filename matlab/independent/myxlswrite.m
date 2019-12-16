function myxlswrite(file, A)
% A is a cell matrix, each element is a string
% Lu Cheng, 25.11.2010

h = fopen(file,'w+');
[nRow nCol] = size(A);

for i=1:nRow
    %tmpLine = '';
    for j=1:nCol-1
        if isnumeric(A{i,j})
            A{i,j} = num2str(A{i,j});
        end
        fprintf(h,'%s\t',A{i,j});
    end
    fprintf(h,'%s\n',A{i,nCol});
end

fclose(h);