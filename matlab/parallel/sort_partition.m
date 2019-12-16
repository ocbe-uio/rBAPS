function srt_partition = sort_partition(partition)
% SORT_PARTITION sorts a given partition (row vector) into the canonical order, where every
% new class has the smallest possible index.

% input:
%    partition is a column vector.
% output:
%    srt_partition is a row vector.

n_classes=max(partition);
srt_partition=zeros(1,n_classes);
for i=1:n_classes
    nonz=find(partition);
    here=find(partition==partition(nonz(1)));
    srt_partition(here)=i;
    partition(here)=0;
end