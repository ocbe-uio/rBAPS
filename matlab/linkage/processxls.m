function [data, component_mat, popnames] = processxls(filename)
%
% a bug in line 64-68 was fixed
data = [];
component_mat = [];
popnames = [];
try
    [A,B] = xlsread(filename);
catch
    display('*** ERROR: Wrong Excel format');
    return
end

if size(A,2)~=1   % more than one columns containing numeric ST values
    display('*** ERROR: multiple columns of numeric values');
    data = []; component_mat = []; popnames = [];
    return
end

if size(A,1)~=size(B,1)-1
    display('*** ERROR: Wrong format');
    data = []; component_mat = []; popnames = [];   
    return
end

B = deblank(B); % remove any trailing blanks
nstrains = size(B,1)-1;
nheader = size(B,2);
for i = 1:nheader
    if strcmpi('ST',B{1,i}) ix_ST = i; end
    if strcmpi('Strain', B{1,i}) || strcmpi('Isolate',B{1,i})
        ix_Strain = i; 
    end
end
if ~exist('ix_ST') 
    display('*** ERROR: ST column needed');
    data = []; component_mat = []; popnames = [];
    return
end

if ~exist('ix_Strain')
    ix_gene = setdiff([1:nheader],ix_ST);
else
    ix_gene = setdiff([1:nheader],[ix_ST ix_Strain]);
end
    
ngenes = length(ix_gene);

C = cell(nstrains,ngenes);
if ~isempty(A)
    for i=1:nstrains
        B{i+1,ix_ST}=num2str(A(i));
        for j=1:ngenes
            C{i,j}=uint16(i_encode_n(B{i+1,ix_gene(j)})); % save the memory.
        end
    end
end
genesize=cellfun('size',C(1,:),2);
data=cell2mat(C);
data=[data uint16([1:nstrains]')]; 
component_mat = zeros(ngenes,max(genesize));
cum = cumsum(genesize);
component_mat(1,[1:genesize(1)]) = [1:cum(1)];
for i=2:ngenes
    component_mat(i,[1:genesize(i)]) = [(cum(i-1)+1):cum(i)];
end

if ~exist('ix_Strain')
    popnames = num2cell(B([2:end],ix_ST));
else  % store the strain names only
    popnames = num2cell(B([2:end],ix_Strain)); 
end
popnames(:,2)=num2cell([1:nstrains]');

display('---------------------------------------------------');
display(['# of strains: ', num2str(nstrains)]);
display(['# of genes: ', num2str(ngenes)]);