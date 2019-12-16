function summat = sumCell(counts)
nbatches = size(counts,2);
sumcell = cell(1,nbatches);
[nalleles, nloci, ninds] = size(counts{1});
summat = zeros(nalleles, nloci,'uint16');
for i = 1:nbatches
    sumcell{i} = uint16(sum(counts{i},3)); % sum as double format.
    summat = summat + sumcell{i};
end
clear counts;
