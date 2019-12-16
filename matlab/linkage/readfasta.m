function data = readfasta(varargin)
MAXNAME = 200;

if nargin == 0
    [filename, pathname] = uigetfile( ...
        {'*.fasta;*.fas;*.txt', 'FASTA Format Files (*.fasta, *.fas, *.txt)';
        '*.*',  'All Files (*.*)'}, ...
        'Pick a FASTA file');
    if ~(filename), aln=[]; return; end
    filename=[pathname,filename];
end

if nargin == 1 
    [filename,pathname] = uigetfile( ...
        {'*.fasta;*.fas;*.txt', 'FASTA Format Files (*.fasta, *.fas, *.txt)';
        '*.*',  'All Files (*.*)'}, ...
        sprintf('Pick the FASTA file for gene %s',varargin{1}) );
    if ~(filename), aln=[];
        data = [];
        return;
    end
    filename=[pathname,filename];
end

if nargin == 2
    [filename,pathname] = uigetfile( ...
        {'*.fasta;*.fas;*.txt', 'FASTA Format Files (*.fasta, *.fas, *.txt)';
        '*.*',  'All Files (*.*)'}, ...
        sprintf('Pick the FASTA file for gene %s',varargin{1}) );
    if ~(filename), aln=[];
        data = [];
        return;
    end
    filename=[pathname,filename];
    chosen_index = varargin{2};
end

if nargin < 3
    % [seqtype, geneticcode]=selectSeqTypeAndGeneticCode;
    seqtype = 2;
    geneticcode = 1;
    if (isempty(seqtype)|isempty(geneticcode)), aln=[],data=[]; return; end
end


pause(0.0001);
if ~ischar(filename)
    error('BAPS:InvalidInput','Input must be a character array')
    data = [];
    return;
end

if ~(exist(filename,'file') | exist(fullfile(cd,filename),'file')),
    %  is a valid filename ?
    error('BAPS:InvalidInput','Input must be a valid file')
    data = [];
    return;
end

file = fopen(filename, 'r');
display('---------------------------------------------------');
display(['Reading fasta sequence from: ',filename,'...']);
display('---------------------------------------------------');
% Now we are looking for the maximum length of the sequence
n=0;    % the number of sequences
m=0;    % the maximum length
cm = 0; % current sequence length

while 1
    [x,nr] = fscanf(file,'%c',1);
    if nr == 0 break; end;
    if x =='>'  % new sequence started
        if cm ~=m & m >0
            fprintf(['*** ERROR: Different sequence length found in allelic ','%d','.\n'],n+1);
            data = [];
            return;
        end
        if cm > m m=cm; end;
        cm = 0;
        fgets(file);
        n=n+1;
    else
        if isletter(x) | x=='-' | x == '?'
            cm=cm+1;
        end;
    end;
end

if cm > m m=cm; end;

% go throught the file
if  (m==0 | n==0)
    display(['*** ERROR: Unmatched data for gene ' varargin{1}]);
    data = [];
    return;
end

Ss = char(m); S = [];
str = zeros(1,MAXNAME);
sizes = zeros(1,n);
frewind(file);
% names=[];
names={};
i=0;j=1;
id = 0;
while 1
    [x,nr] = fscanf(file,'%c',1);
    if nr == 0 
        break; 
    end;
    if x =='>'  % new sequence started
        if i~= 0 % save the sequence
            [x, sizes(i)]=size(Ss);
            S=strvcat(S,Ss);
            Ss = []; Ss = char(m);
        end;
        str=fgetl(file); % read the name, we remove the '>' symbol
        % names=strvcat(names,str);
        % pos=find(str==' ');
        % if ~(isempty(pos))
        %    str=str(1:pos(1,1));
        % end
        i=i+1;
        names{i}=str;
        
        if nargin == 1
            if isempty(findstr(str, varargin{1}))
                display(['*** ERROR: Unmatched data for gene ' varargin{1}]);
                data = [];
                return
            end
        end
        % disp(['Processing in: ' str]);
        id = id + 1;
        j=1;
    else
        if isletter(x) | x== '-' | x=='?'
            % processing the sequence symbol
            Ss(j) = upper(x);
            j=j+1;
        end;
    end;
end
S=strvcat(S,Ss);
[x, sizes(i)]=size(Ss);
if exist('chosen_index','var')
   S = S(chosen_index,:);
   names = names(chosen_index);
end
aln.seqtype = seqtype;
aln.geneticcode = geneticcode;
aln.seqnames = names;
aln.seq = S;
aln = encodealn(aln);
data = aln.seq;
% if nargin == 1
%     if isempty(findstr(names{1},varargin{1}))
%         disp(['*** ERROR: The file does not contain the required gene ' varargin{1}]);
%         data = [];
%         return;
%     end
% end
display(['# of allelic types: ' num2str(size(aln.seq,1))]);
display(['# of nucleotides: ' num2str(size(aln.seq,2))]);
% ----------------------
% order_index = [1:size(aln.seq,1)]';
% data = [aln.seq order_index]; % Append the index column in the end.
% ----------------------

fclose(file);

