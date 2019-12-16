function [data, filename] = silentReadBaps(varargin)
% This function is modified from readbaps.m to avoid invoke a dialog in Linux
% The input must be a vaild BAPS sequence file
% modified by Lu Cheng, 29.06.2010

MAXNAME = 200;

if nargin == 0
    [filename, pathname] = uigetfile( ...
        {'*.txt', 'BAPS Sequence Files (*.txt)';
        '*.*',  'All Files (*.*)'}, ...
        'Load BAPS sequence data');
    if ~(filename), data=[]; return; end
    filename=[pathname,filename];
end

if nargin == 1 
%    [filename,pathname] = uigetfile( ...
%        {'*.txt', 'BAPS Sequence Files (*.txt)';
%        '*.*',  'All Files (*.*)'}, ...
%        sprintf('Load the BAPS sequence data for gene %s',varargin{1}) );
%    if ~(filename), data = [];
%        return;
%    end
%    filename = [pathname,filename];
    filename=varargin{1};  % added by Lu Cheng, 29.06.2010
end

if nargin == 2
    [filename,pathname] = uigetfile( ...
        {'*.txt', 'BAPS Sequence Files (*.txt)';
        '*.*',  'All Files (*.*)'}, ...
        sprintf('Load the BAPS sequence file for gene %s',varargin{1}) );
    if ~(filename),data = [];
        return;
    end
    filename=[pathname,filename];
    chosen_index = varargin{2};
end

if nargin < 3
    % [seqtype, geneticcode]=selectSeqTypeAndGeneticCode;
    seqtype = 2;
    geneticcode = 1;
    if (isempty(seqtype)|isempty(geneticcode)), data=[]; return; end
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
display(['Reading BAPS sequence data from: ',filename,'...']);
display('---------------------------------------------------');
% Now we are looking for the maximum length of the sequence
n=0;    % the number of sequences
m=0;    % the maximum length
cm = 0; % current sequence length

while 1
    [x,nr] = fscanf(file,'%c',1);
    if nr == 0 break; end;
    if x==' '  % new sequence started
        if cm ~=m & m >0
            % fprintf(['*** ERROR: Different sequence length found in allelic ','%d','.\n'],n+1);
            disp('***ERROR: Incorrect BAPS sequence data.');
            data = [];
            fclose(file);
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
    % display(['*** ERROR: Unmatched data for gene ' varargin{1}]);
    disp('***ERROR: Incorrect BAPS sequence data.');
    data = [];
    fclose(file);
    return;
end

Ss = char(m); S = [];
str = zeros(1,MAXNAME);
sizes = zeros(1,n);
frewind(file);
% names=[];
names={};
i=1;j=1;
id = 0;
while 1
    [x,nr] = fscanf(file,'%c',1);
    if nr == 0 
        break; 
    end;
    if x==' '  % new sequence started
        if i~= 0 % save the sequence
            % str=x;
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
        names{i}=str2num(str);
        i=i+1;
        %if nargin == 1
        if nargin == 1000   % modified by Lu Cheng, 29.06.2010
            if isempty(findstr(str, varargin{1}))
                display(['*** ERROR: Unmatched data for gene ' varargin{1}]);
                data = [];
                fclose(file);
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
% S=strvcat(S,Ss);
[x, sizes(i)]=size(Ss);
if ~isempty(find(S==' '))
    disp('***ERROR: unequal sequence length.');
    data = [];
    fclose(file);
    return
end

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
% display(['# of allelic types: ' num2str(size(aln.seq,1))]);
% display(['# of nucleotides: ' num2str(size(aln.seq,2))]);
% ----------------------
% order_index = [1:size(aln.seq,1)]';
% data = [aln.seq order_index]; % Append the index column in the end.
% ----------------------
try,
data = [data cell2mat(names)'];
catch,
    disp('*** ERROR: Failed in loading the BAPS data.');
    disp('*** ERROR: Inidividual indices are not numerically sequential.');
    data = [];
end
fclose(file);

