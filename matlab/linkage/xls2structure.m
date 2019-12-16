function  xls2structure(source, target)
% XLS2STRUCTURE converts a MLST data in xls format into a Structure data

[data, component_mat, popnames] = processxls(source);

data = data(:,[1:end-1]);
[ninds, nloci] = size(data);

fid = fopen(target,'w');
if (fid ~= -1)
    for i = 1:ninds
     fprintf(fid,'%s\t', popnames{i,1}{1});
       for j = 1:nloci
           fprintf(fid,'%d\t',data(i,j));
       end
       fprintf(fid,'\n');
    end
end
fclose(fid);
