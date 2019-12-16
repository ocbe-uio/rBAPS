function paras = semiReadScript(script_file)
% This function extracts parameter information from the script file
% Script Command Table
% datafile('train|test','c:\BAPS5\DATA.xls');  Here only .xls and .mat file
% input file is supported.
% savePreproFile('train|test','c:\BAPS5\predata.mat'); 
% setK('16 17 18');
% outputmat('c:\BAPS5\output.mat')
% Lu Cheng, 11.03.2010

paras.train_file_format = [];
paras.train_file_name = [];

paras.save_prepro_train_data = []; paras.save_prepro_train_data = 'No';
paras.train_prepro_file = [];

paras.test_file_format = [];
paras.test_file_name = [];

paras.save_prepro_test_data = []; paras.save_prepro_test_data = 'No';
paras.test_prepro_file = [];

paras.cluster_num_upperbounds = [];

paras.save_results = []; paras.save_results = 'No';
paras.result_file = [];

T = readfile(script_file);

n = length(T);
for i=1:n
    %line = regexprep(T{i},'\s+','');
    line = T{i};
    [res toks] = regexp(line,'(.+)\((.+)\)','once','match','tokens');
    
    if isempty(res)
        continue;
    else
        %toks
        paras = parseCmd(toks{1}, toks{2}, paras);
    end
end

% -------------------------------------------------------------------------
function prog_paras = parseCmd(cmd, paras, prog_paras)
% cmd is the script command
% paras are the parameters of the script command
% prog_paras is a stucture of the global parameters

switch cmd
    case 'datafile'
        paras = regexprep(paras,'\s+','');
        toks = regexp(paras,'''([^,]+)''','tokens');
        option = toks{1}{:};
        filename = toks{2}{:};
        if exist(filename,'file')~=2
            error(cat(2,'File not exist! File: ',filename));
        end
        filetype = getFileType(filename);
        if isequal(option,'train')
            prog_paras.train_file_format = filetype;
            prog_paras.train_file_name = filename;
        elseif isequal(option,'test')
            prog_paras.test_file_format = filetype;
            prog_paras.test_file_name = filename;
        else
            error(cat(2,'Unkown option: ',option,'! Expect train or test.'));
        end
        
    case 'savePreprocFile'
        paras = regexprep(paras,'\s+','');
        toks = regexp(paras,'''([^,]+)''','tokens');
        option = toks{1}{:};
        filename = toks{2}{:};
        
        filetype = getFileType(filename);
        if ~isequal(filetype,'.mat')
            error(cat(2,'The saved file should end with .mat! ',filename));
        end
        
        if isequal(option,'train')
            prog_paras.save_prepro_train_data = 'Yes';
            prog_paras.train_prepro_file = filename;
        elseif isequal(option,'test')
            prog_paras.save_prepro_test_data = 'Yes';
            prog_paras.test_prepro_file = filename;
        else
            error(cat(2,'Unkown option: ',option,'! Expect train or test.'));
        end
    case 'setK'
        prog_paras.cluster_num_upperbounds = paras(2:end-1);
    case 'outputmat'
        filename = paras(2:end-1);
        filetype = getFileType(filename);
        if ~isequal(filetype,'.mat')
            error(cat(2,'The saved file should end with .mat! ',filename));
        end
        prog_paras.save_results = 'Yes';
        prog_paras.result_file = filename;
    otherwise
        error('Can not parse the cmd: %s in the script!', cmd);
end

% -------------------------------------------------------------------------
function filetype = getFileType(filename)
filetype = filename(end-3:end);
if ~isequal(filetype,'.xls') && ~isequal(filetype,'.mat')
    error(cat(2,'Unknown option: ', filename, '! Expect .xls or .mat file'));
end

% -------------------------------------------------------------------------
function T = readfile(filename)
f = fopen(filename,'r');
if f == -1 
    error(cat(2,'*** ERROR: invalid input file: ',filename));
    T = [];
    return
end

i = 1;
while 1
  clear line;
  line = fgetl(f);
  if ~ischar(line), break, end  
  T{i} = line;
  i = i+1;
end
fclose(f);