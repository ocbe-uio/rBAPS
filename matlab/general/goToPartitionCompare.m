function goToPartitionCompare
% GOTOPARTITIONCOMPARE goes to the partition comaparing mode


% Load the partition result
[filename, pathname] = uigetfile('*.*', 'Load the priori specified partitions');
if (sum(filename)==0) || (sum(pathname)==0)
    return;
end
disp('---------------------------------------------------');
disp('In loading the partition result...');
try
    c = load([pathname filename]);
catch
    fprintf(1,'***ERROR: incorrect partition result.\n');
    return
end
if sum(c(1,:))~=1
    fprintf(1,'***ERROR: invalid prior density.\n');
    return
else
    prior = c(1,:);
    c = c([2:end],:);
end
[ninds, npartitions] = size(c);
fprintf(1,'# of sampling units: %d\n', ninds);
fprintf(1,'# of partitions in comparision: %d\n', npartitions);



h1 = findobj('Tag','partitioncompare_menu');

% Choose data type and model type
items(1).name = 'Model:';
items(1).default = 1;
items(1).indent = 1;
items(1).values = {'Independent';'Spatial';'Linkage'};
items(1).linked = [2 3 4];

items(2).name = 'Data type';
items(2).indent = 1;
items(2).values = {1};
items(3).name = 'individual level';
items(3).default = 1;
items(3).exclusive = 4;
items(3).indent = 2;
items(4).name = 'group level';
items(4).default = 0;
items(4).exclusive = 3;
items(4).indent = 2;

title = 'Specify data and model types';
out = CSEFlagDialog1(items, title);
if isempty(out)
    disp(['cancelled.']);
    return
end

userdata.partitions = c;
userdata.logmls = [];
set(h1,'UserData',userdata);

try
    if out(1).answer == 1
        if out(3).answer == 1 % independent model with individual level
            disp(['Model: Independent clustering - Individual level']);
            greedyMix(-1);
        else % independent model with group level
            disp(['Model: Independent clustering - Group level']);
            greedyPopMix();
        end

    elseif out(1).answer == 2
        if out(3).answer == 1 % spatial model with individual level
            disp(['Model: Spatial clustering - Individual level']);
            spatialMixture;
        else
            disp(['Model: Spatial clustering - Group level']);
            spatialPopMixture;
        end
        
    elseif out(1).answer == 3
            disp(['Model: Linkage clustering']);
            linkageMixture_speed;

    end

    userdata = get(h1,'UserData');
    logmls = userdata.logmls;
    if isempty(logmls)
        disp(['*** ERROR: program stopped.']);
        set(findobj('Tag','filename1_text'),'String',[]);
    else
        diary('baps5_partitioncompare.out');
        disp('---------------------------------------------------');
        disp(['Partition  Prior  LogLikelihood  Posterior']);
        posterior = zeros(1, npartitions);
        sum_posterior = exp(logmls)*prior';
        if sum_posterior == 0 % meaning that one partition dominates
            dominate_partition = find(logmls==max(logmls));
            posterior(dominate_partition) = 1;
        else
            for i = 1:npartitions
            posterior(i) = exp(userdata.logmls(i))*prior(i)/sum_posterior;
            end
        end
    
        for i = 1:npartitions
            disp(['        ' ownNum2Str(i) '   ' ownNum2Str(prior(i)) '  ' ...
                   ownNum2Str(logmls(i)) '  ' ownNum2Str(posterior(i))]);
        end
        diary off
        
        save_preproc = questdlg('Do you wish to save the partition compare result?',...
            'Save result?',...
            'Yes','No','Yes');
        if isequal(save_preproc,'Yes');
            waitALittle;
            [filename, pathname] = uiputfile('*.txt','Save result as');
            if (sum(filename)==0) || (sum(pathname)==0)
                % Cancel was pressed
                return;
            else
                % copy 'baps4_output.baps' into the text file with the same name.
                if exist('baps5_partitioncompare.out','file')
                    copyfile('baps5_partitioncompare.out',[pathname filename])
                    delete('baps5_partitioncompare.out')
                    disp('result saved.');
                else
                    disp('*** ERROR: result cannot be saved.');
                end
            end;
        end
    end
catch
   disp('*** ERROR: incorrect format. Check model and data types');
end
set(h1,'UserData',[]);
