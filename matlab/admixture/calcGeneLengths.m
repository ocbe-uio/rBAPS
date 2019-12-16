function gene_lengths = calcGeneLengths(component_mat)
[ngenes, y] = size(component_mat);
gene_lengths = zeros(ngenes,1);
for i = 1:ngenes
    gene_lengths(i) = length(find(component_mat(i,:)>0));
end
