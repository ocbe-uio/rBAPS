function linkageMixture_speed
base = findobj('Tag','base_figure');


% check whether fixed k mode is selected
h0 = findobj('Tag','fixk_menu');
fixedK = get(h0, 'userdata');

if fixedK
    if ~(fixKWarning == 1) % call function fixKWarning
        return
    end
end

% check whether partition compare mode is selected
h1 = findobj('Tag','partitioncompare_menu');
partitionCompare = get(h1, 'userdata');

% Data handling
input_type = questdlg('Specify the format of your data: ',...
    'Specify Data Format', ...
    'MLST-format', 'BAPS-format','Pre-processed data', 'MLST-format');

switch input_type

    case 'MLST-format'

        %waitALittle;
        mlst_type = questdlg('Choose data type: ',...
            'Specify MLST format', ...
            'Separate allelic profiles(TXT)', ...
            'Concatenate allelic sequences(EXCEL)','Separate allelic profiles(TXT)');
        switch mlst_type
            case 'Concatenate allelic sequences(EXCEL)'
                %waitALittle
                setWindowOnTop(base,'false')
                [filename,pathname] = uigetfile('*.xls','Load new concatenate sequence profile(*.xls)');
                if (sum(filename)==0) || (sum(pathname)==0)
                    return;
                end
                display('---------------------------------------------------');
                display(['Reading sequence profile from: ',[pathname filename],'...']);

                h0 = findobj('Tag','filename1_text');
                set(h0,'String',filename);
                [data, component_mat, popnames] = processxls([pathname filename]);

                if isempty(data)
                    display('*** ERROR: Failed in loading the data');
                    return;
                end
     
            case 'Separate allelic profiles(TXT)'
                % Ask the allelic profile
                setWindowOnTop(base,'false')
                [filename,pathname] = uigetfile('*.pl;*.txt','Load new allelic profile(*.pl, *.txt)');
                if (sum(filename)==0) || (sum(pathname)==0)
                    return;
                end
                display('---------------------------------------------------');
                display(['Reading allelic profile from: ',[pathname filename],'...']);

                h0 = findobj('Tag','filename1_text');
                set(h0,'String',filename);

                % Preprocess the profile
                output = processprofile([pathname filename]);
                headercount = size(output,2);
                flag = zeros(1,headercount);
                for i = 1:headercount
                    if strcmpi('ST',output{1,i}) 
                        flag(i) = 1;
                    else
                        if strcmpi('Isolate', output{1,i}) || strcmpi('Strain',output{1,i})
                            flag(i) = 2;
                        else
                            if strcmpi('Species', output{1,i}) 
                                flag(i) = 3;
                            end
                        end
                    end
                end

                if ~any(flag)
                    h = errordlg(['Loading of the specified file was unsuccessful. ' ...
                        'Please see the tutorial to find out the correct file ' ...
                        'format.'] ,'Error','modal');
                    handle = [h,base];
                    setWindowOnTop(handle,{'true','true'})

                    uiwait(h);
                    fprintf(1,'\n*** ERROR: Failed in loading the allelic profile.\n');
                    return;
                end

                index = (1:size(output,1)-1)';


                % Species selection
                species_loc = find(flag==3);
                if ~isempty(species_loc)
                    species_col = output((2:end), species_loc);
                    [species_str, m, n] = unique(species_col);

                    [s1,v1] = listdlg('PromptString','Select species:',...
                        'SelectionMode','multiple',...
                        'Name','Select Species',...
                        'ListString',species_str');
                    removed = setdiff(n,s1);

                    if ~v1 || isempty(s1) 
                        dispCancel;return
                    end
                else
                    n = index;
                    removed = 0;
                end

                % Isolate/Strain selection
                isolate_loc = find(flag==2);
                if ~isempty(isolate_loc)
                    isolate_col = output((2:end), isolate_loc);
                    isolate_str = isolate_col(logical(~ismember(n,removed)));
                    index(ismember(n,removed)) = [];

                    [s2,v2] = choosebox('Name','Select Isolates','PromptString',...
                        'Isolates in the sample(Ctrl+A to select all):','SelectString', ...
                        'Isolates you have selected','ListString', isolate_str', ...
                        'InitialValue',1);
                    if isempty(s2) || ~v2 
                        dispCancel; return
                    end
                else
                    % ST selection
                    isolate_loc = find(flag==1);
                    isolate_col = output((2:end), isolate_loc);
                    % isolate_str = isolate_col(find(~ismember(n,removed)));
                    isolate_str = isolate_col(logical(~ismember(n,removed)));
                    index(ismember(n,removed)) = [];

                    [s2,v2] = choosebox('Name','Select STs','PromptString',...
                        'STs in the sample (Ctrl+A to select all):','SelectString', ...
                        'STs you have selected','ListString', isolate_str');
                    if isempty(s2) || ~v2 
                        dispCancel; return
                    end
                end

                % Read the data
                isolate_index = index(s2);
                if ~isempty(isolate_loc)
                    popnames(:,1) = isolate_str(s2);
                else
                    popnames(:,1) = num2cell(index(s2));
                end
                popnames(:,2) = num2cell([1:length(isolate_index)]');

                % remove empty elements
                strmat = output([isolate_index+1]',find(flag==0));
                [i,j] = find(cellfun('isempty',strmat));
                ij = [i j];
                for k=1:length(i)
                    strmat(ij(k,1),ij(k,2))={''};
                end
                [i,j] = find(strcmp(strmat,''));
                ij = [i j];
                for k=1:length(i)
                    strmat(ij(k,1),ij(k,2))={'0'};
                end

                data_allele = zeros(size(strmat));
                for i = 1:size(strmat,1)
                    for j = 1:size(strmat,2)
                        data_allele(i,j) =  str2num(char(strmat(i,j)));
                    end
                end

                % remove columns containing empty values
                realgene = find(all(data_allele));
                data_allele = data_allele(:,realgene);
                partition_index = 1:size(strmat,1);
                data_allele = [data_allele partition_index'];
                genename = output(1,find(flag==0));
                genename = genename(realgene);


                if isempty(genename) || isempty(data_allele)
                    msgbox(['Loading of the specified file was unsuccessful. ' ...
                        'Please see the tutorial to find out the correct file ' ...
                        'format.'] ,'Error', ...
                        'error');
                    fprintf(1,'\n*** ERROR: Failed in loading the allelic profile.\n');
                    return;
                else
                    display('---------------------------------------------------');
                    display(['# of strains: ', num2str(size(data_allele,1))]);
                    display(['# of genes: ', num2str(size(data_allele,2)-1)]);
                end

                %waitALittle;
                % Ask the individual gene sequence
                %             [isOK,returnvalue] = askSeq;
                %             if (isOK & returnvalue~= 1)  % allelic profile loaded only
                %                 data = data_allele;
                %                 component_mat = [1:size(data,2)-1]'; % assume independence
                %             elseif
                %               isOK & returnvalue == 1,
                [s3,v3] = listdlg('PromptString','Select genes:',...
                    'SelectionMode','multiple',...
                    'Name','Select Genes',...
                    'ListString',genename);

                if isempty(s3) || ~v3
                    dispCancel;
                    return
                else
                    m = size(s3,2); % number of genes
                    % data_seq = cell(m,1);
                    data = [];
                    genesize = zeros(1,m);
                    for i=1:m
                        %waitALittle;
                        data_gene = readfasta(genename{s3(i)}); % read fasta
                        if (isempty(data_gene))
                            display('*** ERROR: Failed in loading the sequence data');
                            return;
                        else
                            data_gene(:,find(sum(data_gene)==0)) = []; % NB! remove all the gaps
                            % data_seq{i} = data_gene;
                            selected_data = data_gene(data_allele(:,s3(i))',:); % Store only those selected strains

                            emptyloci = find(all(selected_data(:,[1:end])<0));
                            if ~isempty(emptyloci)
                                disp('Removing empty loci...');
                            end
                            selected_data(:,emptyloci) = []; % remove empty loci
                            data = [data selected_data];
                        end
                        genesize(i) = size(selected_data,2); % NB! could be different than the original gene length
                    end
                    data = [data data_allele(:,end)]; % add the index
                    % determine the component matrix
                    component_mat = zeros(m,max(genesize));
                    cum = cumsum(genesize);
                    component_mat(1,1:genesize(1)) = 1:cum(1);
                    for i=2:m
                        component_mat(i,1:genesize(i)) = (cum(i-1)+1):cum(i);
                    end
                end
                %             elseif isOK ==0
                %                 return
                %            end
            otherwise
                return
        end

        %waitALittle;
        display('---------------------------------------------------');
        fprintf(1,'Preprocessing the data ...');

        % Make the missing data complete
        % missing values are denoted as -999
        data = uint16(data);
        % data = uint8(data);
        data = makecomplete(data);
        if isempty(data)
            display('*** ERROR: Failed in completing the missing data');
            return;
        end

        isRational = isTheLoadedNewDataRational(data);
        if isRational == 0
            return;
        end

        [data, rowsFromInd, alleleCodes, noalle, adjprior, priorTerm] = ...
            handleData(data);


        % Distance between individuals is computed as if the loci are
        % independent.
        [Z,dist] = newGetDistances(data,rowsFromInd);
        fprintf(1,'Finished.\n');
        ninds = max(data(:,end));
        popnames = fixPopnames(popnames, ninds);

        c.data = uint16(data); c.rowsFromInd = rowsFromInd; c.alleleCodes = alleleCodes;
        c.noalle = noalle; c.adjprior = adjprior; c.priorTerm = priorTerm;
        c.popnames = popnames; c.component_mat = component_mat;
        c.dist = dist; c.Z = Z;     

        %waitALittle;
        save_preproc = questdlg('Do you wish to save the pre-processed data?',...
            'Save pre-processed data?',...
            'Yes','No','Yes');
        if isequal(save_preproc,'Yes');
            %waitALittle;
            [filename, pathname] = uiputfile('*.mat','Save pre-processed data as');
            if (sum(filename)==0) || (sum(pathname)==0)
                return;
            else
                kokonimi = [pathname filename];
%                 save(kokonimi,'c');
                save(kokonimi,'c','-v7.3'); % added by Lu Cheng, 08.06.2012
            end
        end;

        %waitALittle;
        linkage_model = questdlg('Specify the linkage model',...
            'Specify the linkage model?',...
            'Linear','Codon', 'Independent', 'Linear');
        if isequal(linkage_model,'Linear')
            linkage_model = 'linear';
            display('Linear model was selected.');
        elseif isequal(linkage_model,'Codon')
            linkage_model = 'codon';
            display('Codon model was selected.');
        elseif isequal(linkage_model,'Independent')
            display('Independent model was selected.');
            c.data = double(c.data);
            greedyMix(c);
            return;
        else
            dispCancel;
            return;
        end;

        % Data transformation
        fprintf(1,'Transforming the data ...');
        index = data(:,end);
        [data_clique, data_separator, noalle_clique, noalle_separator] = ...
            transform4(data, component_mat, linkage_model);
        data_clique = [data_clique index];
        data_separator = [data_separator index];

        % Count the data
        [counts_cq, nalleles_cq, prior_cq, adjprior_cq, genotypes_cq]...
            = allfreqsnew2(data_clique, double(noalle_clique));
        [counts_sp, nalleles_sp, prior_sp, adjprior_sp, genotypes_sp]...
            = allfreqsnew2(data_separator, double(noalle_separator));

        counts_cq = uint8(counts_cq);
        counts_sp = uint8(counts_sp);
        fprintf(1,'Finished.\n');

        clear data_clique data_separator

        save_preproc = questdlg('Do you wish to save the fully pre-processed data?',...
            'Save pre-processed data?',...
            'Yes','No','Yes');
        if isequal(save_preproc,'Yes');
            %waitALittle;
            [filename, pathname] = uiputfile('*.mat','Save fully pre-processed data as');
            if (sum(filename)==0) || (sum(pathname)==0)
                return;
            else
                kokonimi = [pathname filename];
                c.counts_cq = counts_cq; c.adjprior_cq = adjprior_cq;
                c.counts_sp = counts_sp; c.adjprior_sp = adjprior_sp;
                c.linkage_model = linkage_model;
%                 save(kokonimi,'c');
                save(kokonimi,'c','-v7.3'); % added by Lu Cheng, 08.06.2012
            end
        end;
        clear c;

    case 'BAPS-format'

        input_type = questdlg('Specify the format of your data: ',...
            'Specify BAPS Format', ...
            'BAPS sequence data', 'BAPS numeric data', 'BAPS sequence data');
        switch input_type
            case 'BAPS numeric data'
                %waitALittle;
                setWindowOnTop(base,'false')
                [filename,pathname] = uigetfile('*.txt', 'Load BAPS numeric data');
                if (sum(filename)==0) || (sum(pathname)==0)
                    %cancel was pressed; do nothing.
                    return;
                end;

                display('---------------------------------------------------');
                display(['Reading BAPS numeric data from: ',[pathname filename],'...']);

                try
                    data = load([pathname filename]);
                catch
                    disp('*** ERROR: Incorrect BAPS numerical data.');
                    return
                end
            case 'BAPS sequence data'
                waitALittle;%waitALittle;
                [data, filename] = readbaps;
                if isempty(data) 
                    return
                end
            otherwise
                return;
        end

        h0 = findobj('Tag','filename1_text');
        set(h0,'String',filename); clear h0;
        waitALittle;%waitALittle;
        input_pops = questdlg(['When using data which are in BAPS-format, '...
            'you can specify the sampling populations of the individuals by '...
            'giving two additional files: one containing the names of the '...
            'populations, the other containing the indices of the first '...
            'individuals of the populations. Do you wish to specify the '...
            'sampling populations?'], ...
            'Specify sampling populations?',...
            'Yes', 'No', 'No');
        if isequal(input_pops,'Yes')
            %waitALittle;
            display('Reading name and index files...');
            setWindowOnTop(base,'false')
            [namefile, namepath] = uigetfile('*.txt', 'Load population names');
            if namefile==0
                kysyToinen = 0;
            else
                kysyToinen = 1;
            end
            if kysyToinen==1
                %waitALittle;
                setWindowOnTop(base,'false')
                [indicesfile, indicespath] = uigetfile('*.txt', 'Load population indices');
                if indicesfile==0
                    % popnames = [];
                    dispCancel;
                    return
                else
                    popnames = initPopNames([namepath namefile],[indicespath indicesfile]);
                end
            else
                % popnames = [];
                dispCancel;
                return
            end
        else
            % popnames = [];
            popnames(:,1) = num2cell(unique(data(:,end)));
            popnames(:,2) = popnames(:,1);
            ninds = max(data(:,end));
            popnames = fixPopnames(popnames, ninds);
        end

        % check that popnames is correct
        if isempty(popnames)
            display('*** ERROR: error in reading popnames.')
            return
        end
        
        data = uint16(data);
        % Check that the data is rational:
        isRational = isTheLoadedNewDataRational(data);
        if isRational == 0
            msgbox(['Loaded file contained incorrect data. The last column of the ' ...
                'data file must contain sampling unit identifiers. Identifier specifies the ' ...
                'unit from which the genetic data on that particular row was collected. ' ...
                'Identifiers must be positive integers. If the biggest unit identifier is i.e. 27, there must ' ...
                'be at least one row for each unit 1-27'] ,'Error', ...
                'error');
            disp('*** ERROR: Failed in loading the BAPS data.');
            return;
        else
            display(['# of haplotypes: ', num2str(size(data,1))]);
            display(['# of loci: ', num2str(size(data,2)-1)]);
            % display('Finished.');
        end;

        % Check if the data is discrete or continuous
        if any(any(fix(data)~=data))
            disp('Found decimal numbers. Continuous model will be used.');
%             input_type = questdlg('Choose the method for the continuous data: ',...
%                 'Specify the method', ...
%                 'BEC', 'Gibbs sampling with WINBUGS','BEC');
%             switch input_type
%                 case 'BEC'
%                     disp('Using the BEC mixture model...');
%                     becMixture(data, popnames);
%                 case 'Gibbs sampling with WINBUGS'
%                     disp('Preparing the WINBUGS code...');
%                     isok = makeBUGS(data);
%                     if isok disp('Finished.');
% 
%                     end
%                 otherwise
%                     return
%             end
            disp('** CANCELLED: continuous model is under construction.');
            return;
        end
        
        display('---------------------------------------------------');
        fprintf(1,'Preprocessing the data ...');
        % Make the missing data complete
        data = makecomplete(data);
        if isempty(data)
            display('*** ERROR: Failed in completing the missing data');
            return;
        end

        [data, rowsFromInd, alleleCodes, noalle, adjprior, priorTerm] = ...
            handleData(data);

        % Distance between individuals is computed as if the loci are
        % independent.
        [Z,dist] = newGetDistances(data,rowsFromInd);
        fprintf(1,'Finished.\n');

        c.data = uint16(data); c.rowsFromInd = rowsFromInd; c.alleleCodes = alleleCodes;
        c.noalle = noalle; c.adjprior = adjprior; c.priorTerm = priorTerm;
        c.dist = dist; c.popnames = popnames; 
        c.Z = Z;

        input_linkage = questdlg('Do you wish to load the linkage map?',...
            'Load Linkage Map',...
            'Yes','No','Yes');
        if isequal(input_linkage,'Yes');
            display('---------------------------------------------------');
            %%waitALittle;
            setWindowOnTop(base,'false')
            [linkage_filename, linkage_pathname] = uigetfile('*.txt', 'Load Linkage Map');

            if isempty(linkage_filename) && isempty(linkage_pathname)
                return;
            else
                display(['Reading linkage map from: ',[linkage_pathname linkage_filename],'...']);
            end;

            try
                component_mat = load([linkage_pathname linkage_filename]);
            catch
                disp('*** ERROR: Incorrect linkage map.');
                return;
            end

            % Check if the linkage map matches the data
            if (size(data,2)-1) ~= max(component_mat(:))
                msgbox(['Loading of the specified file was unsuccessful. ' ...
                    'The linkage map dose not match with the data.'] ,'Error', ...
                    'error');
                disp('*** ERROR: Failed in loading the linkage map.');
                return;
            else
                display(['# of linkage groups: ', num2str(size(component_mat,1))]);
            end;
            display('---------------------------------------------------');
            h0 = findobj('Tag','filename1_text');
            set(h0,'String',[filename '/' linkage_filename]); clear h0;
            c.component_mat = component_mat;
        else
            display('Independent model was selected.');
            c.data = double(c.data);
            greedyMix(c);
            return;
        end


        save_preproc = questdlg('Do you wish to save the pre-processed data?',...
            'Save pre-processed data?',...
            'Yes','No','Yes');
        if isequal(save_preproc,'Yes');
            %%waitALittle;
            [filename, pathname] = uiputfile('*.mat','Save pre-processed data as');
            if (sum(filename)==0) || (sum(pathname)==0)
                return;
            end
            kokonimi = [pathname filename];
%             save(kokonimi,'c');
            save(kokonimi,'c','-v7.3'); % added by Lu Cheng, 08.06.2012
        end;

        linkage_model = questdlg('Specify the linkage model',...
            'Specify the linkage model?',...
            'Linear','Codon', 'Independent', 'Linear');
        if isequal(linkage_model,'Linear')
            linkage_model = 'linear';
            display('Linear model was selected.');
        elseif isequal(linkage_model,'Codon')
            linkage_model = 'codon';
            display('Codon model was selected.');
        elseif isequal(linkage_model,'Independent')
            display('Independent model was selected.');
            c.data = double(c.data);
            greedyMix(c);
            return;
        else
            dispCancel;
            return;
        end;

        % Data transformation
        % display('---------------------------------------------------');
        fprintf(1,'Transforming the data ...');
        index = data(:,end);
        [data_clique, data_separator, noalle_clique, noalle_separator] = ...
            transform4(data, component_mat, linkage_model);
        data_clique = [data_clique index];
        data_separator = [data_separator index];

        [counts_cq, nalleles_cq, prior_cq, adjprior_cq, genotypes_cq]...
            = allfreqsnew2(data_clique, double(noalle_clique));
        clear data_clique;
        [counts_sp, nalleles_sp, prior_sp, adjprior_sp, genotypes_sp]...
            = allfreqsnew2(data_separator, double(noalle_separator));
        clear data_separator;
        counts_cq = uint8(counts_cq);
        counts_sp = uint8(counts_sp);
        fprintf(1,'Finished.\n');

        save_preproc = questdlg('Do you wish to save the fully pre-processed data?',...
            'Save pre-processed data?',...
            'Yes','No','Yes');
        if isequal(save_preproc,'Yes');
            %waitALittle;
            [filename, pathname] = uiputfile('*.mat','Save fully pre-processed data as');
            if (sum(filename)==0) || (sum(pathname)==0)
                return;
            end
            kokonimi = [pathname filename];
            c.counts_cq = counts_cq; c.adjprior_cq = adjprior_cq;
            c.counts_sp = counts_sp; c.adjprior_sp = adjprior_sp;
            c.linkage_model = linkage_model;
%             save(kokonimi,'c');
            save(kokonimi,'c','-v7.3'); % added by Lu Cheng, 08.06.2012
        end;
        clear c;

    case 'Pre-processed data'
        % This is basically the same format as the "Pre-processed data" in
        % the basic clustering. The only difference is that the file
        % includes also the component_mat
        % %waitALittle;
        setWindowOnTop(base,'false')
        [filename, pathname] = uigetfile('*.mat', 'Load pre-processed data');
        if filename==0
            return;
        end
        display('---------------------------------------------------');
        display(['Reading preprocessed data from: ',[pathname filename],'...']);
        h0 = findobj('Tag','filename1_text');
        set(h0,'String',filename); clear h0;

        struct_array = load([pathname filename]);
        if isfield(struct_array,'c')  %Matlab versio
            c = struct_array.c;
            if ~isfield(c,'dist')
                display('*** ERROR: Incorrect file format');
                return
            end
            clear struct_array;
        elseif isfield(struct_array,'dist')  %Mideva versio
            c = struct_array;
            clear struct_array;
        else
            display('*** ERROR: Incorrect file format');
            return;
        end

        % The following are the same as in the basic clustering
        data = c.data;  popnames = c.popnames; Z = c.Z;
        noalle = c.noalle; adjprior = c.adjprior;
        rowsFromInd = c.rowsFromInd; alleleCodes = c.alleleCodes;
        dist = c.dist; priorTerm = c.priorTerm;

        if ~isfield(c,'component_mat')
            display('*** ERROR: Incorrect file format');
            return
        end
        
        % This is new
        component_mat = c.component_mat;
        data = uint16(data);
        
        display(['# of haplotypes: ', num2str(size(data,1))]);
        display(['# of loci: ', num2str(size(data,2)-1)]);
        display(['# of linkage groups: ', num2str(size(component_mat,1))]);

        if ~isfield(c, 'linkage_model')
            %%%waitALittle;
            % Independent is not an option, since it can be computed with the
            % basic clustering which is much faster
            linkage_model = questdlg('Specify the linkage model',...
                'Specify the linkage model?',...
                'Linear','Codon','Linear');
            if isequal(linkage_model,'Linear')
                linkage_model = 'linear';
                display('Linear model was selected.');
            elseif isequal(linkage_model,'Codon')
                linkage_model = 'codon';
                display('Codon model was selected.');
            else
                dispCancel;
                return;
            end;
            
            clear c; % save the memory usage

            % Data transformation
            fprintf(1,'Transforming the data ...');
            index = data(:,end);
%             [data_clique, data_separator, noalle_clique, noalle_separator] = ...
%                 transform2(data, component_mat, linkage_model);
            [data_clique, data_separator, noalle_clique, noalle_separator] = ...
                transform4(data, component_mat, linkage_model);
            data_clique = [data_clique index];
            data_separator = [data_separator index];

            [counts_cq, nalleles_cq, prior_cq, adjprior_cq, genotypes_cq]...
                = allfreqsnew2(data_clique, double(noalle_clique));
            clear data_clique;
            [counts_sp, nalleles_sp, prior_sp, adjprior_sp, genotypes_sp]...
                = allfreqsnew2(data_separator, double(noalle_separator));
            clear data_separator;
            counts_cq = uint8(counts_cq);
            counts_sp = uint8(counts_sp);
            fprintf(1,'Finished.\n');

            save_preproc = questdlg('Do you wish to save the fully pre-processed data?',...
                'Save pre-processed data?',...
                'Yes','No','Yes');
            if isequal(save_preproc,'Yes');
                %%%waitALittle;
                [filename, pathname] = uiputfile('*.mat','Save fully pre-processed data as');
                if (sum(filename)==0) || (sum(pathname)==0)
                    return;
                end

                kokonimi = [pathname filename]; c.data = data;
                c.counts_cq = counts_cq; c.adjprior_cq = adjprior_cq;
                c.counts_sp = counts_sp; c.adjprior_sp = adjprior_sp;
                c.linkage_model = linkage_model;
                c.rowsFromInd = rowsFromInd;
                c.alleleCodes = alleleCodes;
                c.noalle = noalle;
                c.adjprior = adjprior;
                c.priorTerm = priorTerm;
                c.popnames = popnames;
                c.component_mat = component_mat;
                c.dist = dist; c.Z = Z;
%                 save(kokonimi,'c');
                save(kokonimi,'c','-v7.3'); % added by Lu Cheng, 08.06.2012
            end;
            clear c;
        else
            %Linkage model is specified in the preprocessed file.
            counts_cq = c.counts_cq; adjprior_cq = c.adjprior_cq;
            counts_sp = c.counts_sp; adjprior_sp = c.adjprior_sp;
            linkage_model = c.linkage_model;
            counts_cq = uint8(counts_cq);
            counts_sp = uint8(counts_sp);
            clear c;
            display(['linkage model: ', linkage_model]);
        end

    otherwise
        return;
end


global POP_LOGML;       global PARTITION;
global CQ_COUNTS;       global SP_COUNTS;       %These counts are for populations
global CQ_SUMCOUNTS;    global SP_SUMCOUNTS;    %not for individuals
clearGlobalVars;

c.noalle = noalle; 
c.adjprior = adjprior; %priorTerm = c.priorTerm;
c.rowsFromInd = rowsFromInd;
c.counts_cq = counts_cq; 
c.adjprior_cq = adjprior_cq;
c.counts_sp = counts_sp; 
c.adjprior_sp = adjprior_sp;
c.dist = dist; c.Z = Z;

% partition compare mode
if ~isempty(partitionCompare)
    partitions = partitionCompare.partitions;
    npartitions = size(partitions,2);
    partitionLogml = zeros(1,npartitions);
    for i = 1:npartitions
        % number of unique partition lables
        try
            [cq_counts, cq_sumcounts] = ...
                initialCounts(counts_cq, partitions(:,i));
            [sp_counts, sp_sumcounts] = initialCounts(counts_sp, partitions(:,i));
            partitionLogml(i) = computeLogml(adjprior_cq, adjprior_sp, ...
                                 cq_counts, cq_sumcounts, ...
                                 sp_counts, sp_sumcounts);
        catch
            disp('*** ERROR: unmatched data.');
            return
        end
    end
    % return the logml result
    partitionCompare.logmls = partitionLogml;
    set(h1, 'userdata', partitionCompare);
    return
end

if fixedK
    [logml, npops, partitionSummary] = linkageMix_fixK(c);
else
    [logml, npops, partitionSummary] = linkageMix(c);
end

if logml==1
    return;
end

h0 = findobj('Tag','filename1_text');  inp = get(h0,'String');
h0 = findobj('Tag','filename2_text');
outp = get(h0,'String');

clear c; % save the memory
%This is basically the same as in BAPS 3.
changesInLogml = writeMixtureInfo(logml, counts_cq, counts_sp, adjprior_cq, ...
    adjprior_sp, outp, inp, partitionSummary, popnames, linkage_model, ...
    fixedK);

viewMixPartition(PARTITION, popnames);

% ---------------------------------------------------------------------
% Save the result.
% Jing - 26.12.2005
talle = questdlg(['Do you want to save the mixture populations ' ...
    'so that you can use them later in admixture analysis?'], ...
    'Save results?','Yes','No','Yes');
if isequal(talle,'Yes')
    %%%waitALittle;
    [filename, pathname] = uiputfile('*.mat','Save results as');

    if (sum(filename)==0) || (sum(pathname)==0)
        % Cancel was pressed
        return
    else % copy 'baps4_output.baps' into the text file with the same name.
        if exist('baps4_output.baps','file')
            copyfile('baps4_output.baps',[pathname filename '.txt'])
            delete('baps4_output.baps')
        end
    end

    [sumcounts, counts] = indLociCounts(PARTITION, data, npops, noalle);
    % NB! Index column is removed in data matrix.
    c.PARTITION = PARTITION; c.CQ_COUNTS = CQ_COUNTS; c.CQ_SUMCOUNTS = CQ_SUMCOUNTS;
    c.SP_COUNTS = SP_COUNTS; c.SP_SUMCOUNTS = SP_SUMCOUNTS;
    c.alleleCodes = alleleCodes; c.adjprior_cq = adjprior_cq; c.adjprior_sp = adjprior_sp; c.popnames = popnames;
    c.rowsFromInd = rowsFromInd; c.data = uint16(data); c.npops = npops;
    % c.nalleles_cq = nalleles_cq; c.nalleles_sp = nalleles_sp;
    if strcmp(linkage_model,'linear') % Added on 03.11.06
        c.mixtureType = 'linear_mix';
    elseif strcmp(linkage_model,'codon')
        c.mixtureType = 'codon_mix';
    end
    c.changesInLogml = changesInLogml; % this variable stores the change of likelihoods.
                                       % [ncluster ninds]
                                       % -Added on 02.11.2006

    % The next ones are for the admixture input
    c.COUNTS = counts; c.SUMCOUNTS = sumcounts;
    c.adjprior = adjprior; c.rowsFromInd = rowsFromInd; c.noalle = noalle; c.alleleCodes = alleleCodes;
    
    % The two variables below are for the new linkage admixture model
    c.gene_lengths = calcGeneLengths(component_mat);
    
    % The logml is saved for parallel computing
    c.logml = logml;
    
    fprintf(1,'Saving the result...')
    try
%         save([pathname filename], 'c');
        save([pathname filename], 'c', '-v7.3'); % added by Lu Cheng, 08.06.2012
        fprintf(1,'Finished.\n');
    catch
        display('*** ERROR in saving the result.');
    end
else
    if exist('baps4_output.baps','file')
        delete('baps4_output.baps')
    end
end
% -----------------------------------------------------------------------


%--------------------------------------------------------------------------
% The next three functions are for computing the initial partition
% according to the distance between the individuals

function initial_partition=admixture_initialization(nclusters,Z)
T=cluster_own(Z,nclusters);
initial_partition=T;

%--------------------------------------------------------------------------
function T = cluster_own(Z,nclust)
% true=logical(1);
% false=logical(0);

maxclust = nclust;
% Start of algorithm
m = size(Z,1)+1;
T = zeros(m,1);
% maximum number of clusters based on inconsistency
if m <= maxclust
    T = (1:m)';
elseif maxclust==1
    T = ones(m,1);
else
    clsnum = 1;
    for k = (m-maxclust+1):(m-1)
        i = Z(k,1); % left tree
        if i <= m % original node, no leafs
            T(i) = clsnum;
            clsnum = clsnum + 1;
        elseif i < (2*m-maxclust+1) % created before cutoff, search down the tree
            T = clusternum(Z, T, i-m, clsnum);
            clsnum = clsnum + 1;
        end
        i = Z(k,2); % right tree
        if i <= m  % original node, no leafs
            T(i) = clsnum;
            clsnum = clsnum + 1;
        elseif i < (2*m-maxclust+1) % created before cutoff, search down the tree
            T = clusternum(Z, T, i-m, clsnum);
            clsnum = clsnum + 1;
        end
    end
end

function T = clusternum(X, T, k, c)
m = size(X,1)+1;
while(~isempty(k))
    % Get the children of nodes at this level
    children = X(k,1:2);
    children = children(:);

    % Assign this node number to leaf children
    t = (children<=m);
    T(children(t)) = c;

    % Move to next level
    k = children(~t) - m;
end

%--------------------------------------------------------------------------

function dist2 = laskeOsaDist(inds2, dist, ninds)
% Muodostaa dist vektorista osavektorin, joka sis‰lt‰‰ yksilˆiden inds2
% v‰liset et‰isyydet. ninds=kaikkien yksilˆiden lukum‰‰r?

ninds2 = length(inds2);
apu = zeros(nchoosek(ninds2,2),2);
rivi = 1;
for i=1:ninds2-1
    for j=i+1:ninds2
        apu(rivi, 1) = inds2(i);
        apu(rivi, 2) = inds2(j);
        rivi = rivi+1;
    end
end
apu = (apu(:,1)-1).*ninds - apu(:,1) ./ 2 .* (apu(:,1)-1) + (apu(:,2)-apu(:,1));
dist2 = dist(apu);

%--------------------------------------------------------------------------

function Z = computeLinkage(Y, method)
[k, n] = size(Y);
m = (1+sqrt(1+8*n))/2;
if k ~= 1 || m ~= fix(m)
    error('The first input has to match the output of the PDIST function in size.');
end
if nargin == 1 % set default switch to be 'co'
    method = 'co';
end
method = lower(method(1:2)); % simplify the switch string.
% monotonic = 1;
Z = zeros(m-1,3); % allocate the output matrix.
N = zeros(1,2*m-1);
N(1:m) = 1;
n = m; % since m is changing, we need to save m in n.
R = 1:n;
for s = 1:(n-1)
    X = Y;
    [v, k] = min(X);
    i = floor(m+1/2-sqrt(m^2-m+1/4-2*(k-1)));
    j = k - (i-1)*(m-i/2)+i;
    Z(s,:) = [R(i) R(j) v]; % update one more row to the output matrix A
    I1 = 1:(i-1); I2 = (i+1):(j-1); I3 = (j+1):m; % these are temp variables.
    U = [I1 I2 I3];
    I = [I1.*(m-(I1+1)/2)-m+i i*(m-(i+1)/2)-m+I2 i*(m-(i+1)/2)-m+I3];
    J = [I1.*(m-(I1+1)/2)-m+j I2.*(m-(I2+1)/2)-m+j j*(m-(j+1)/2)-m+I3];

    switch method
        case 'si' %single linkage
            Y(I) = min(Y(I),Y(J));
        case 'av' % average linkage
            Y(I) = Y(I) + Y(J);
        case 'co' %complete linkage
            Y(I) = max(Y(I),Y(J));
        case 'ce' % centroid linkage
            K = N(R(i))+N(R(j));
            Y(I) = (N(R(i)).*Y(I)+N(R(j)).*Y(J)-(N(R(i)).*N(R(j))*v^2)./K)./K;
        case 'wa'
            Y(I) = ((N(R(U))+N(R(i))).*Y(I) + (N(R(U))+N(R(j))).*Y(J) - ...
                N(R(U))*v)./(N(R(i))+N(R(j))+N(R(U)));
    end
    J = [J i*(m-(i+1)/2)-m+j];
    Y(J) = []; % no need for the cluster information about j.

    % update m, N, R
    m = m-1;
    N(n+s) = N(R(i)) + N(R(j));
    R(i) = n+s;
    R(j:(n-1))=R((j+1):n);
end

%--------------------------------------------------------------------------

function changes = computeChanges(ind, adjprior_cq, adjprior_sp, ...
    indCqCounts, indSpCounts)
% Computes changes in log-marginal likelihood if individual ind is
% moved to another population
%
% Input:
% ind - the individual to be moved
% adjprior_cq & _sp  - adjpriors for cliques and separators
% indCqCounts, indSpCounts - counts for individual ind
%
% Output:
% changes - table of size 1*npops. changes(i) = difference in logml if
% ind is move to population i.

global CQ_COUNTS;   global CQ_SUMCOUNTS;
global SP_COUNTS;   global SP_SUMCOUNTS;
global PARTITION;   global POP_LOGML;
npops = size(CQ_COUNTS,3);
changes = zeros(npops,1);

i1 = PARTITION(ind);
i1_logml = POP_LOGML(i1);
sumCq = uint16(sum(indCqCounts,1));
sumSp = uint16(sum(indSpCounts,1));

CQ_COUNTS(:,:,i1) = CQ_COUNTS(:,:,i1)-indCqCounts;
CQ_SUMCOUNTS(i1,:) = CQ_SUMCOUNTS(i1,:)-sumCq;
SP_COUNTS(:,:,i1) = SP_COUNTS(:,:,i1)-indSpCounts;
SP_SUMCOUNTS(i1,:) = SP_SUMCOUNTS(i1,:)-sumSp;

new_i1_logml = computePopulationLogml(i1, adjprior_cq, adjprior_sp);

CQ_COUNTS(:,:,i1) = CQ_COUNTS(:,:,i1)+indCqCounts;
CQ_SUMCOUNTS(i1,:) = CQ_SUMCOUNTS(i1,:)+sumCq;
SP_COUNTS(:,:,i1) = SP_COUNTS(:,:,i1)+indSpCounts;
SP_SUMCOUNTS(i1,:) = SP_SUMCOUNTS(i1,:)+sumSp;


i2 = [1:i1-1 , i1+1:npops];
i2_logml = POP_LOGML(i2);

CQ_COUNTS(:,:,i2) = CQ_COUNTS(:,:,i2)+repmat(indCqCounts, [1 1 npops-1]);
CQ_SUMCOUNTS(i2,:) = CQ_SUMCOUNTS(i2,:)+repmat(sumCq,[npops-1 1]);
SP_COUNTS(:,:,i2) = SP_COUNTS(:,:,i2)+repmat(indSpCounts, [1 1 npops-1]);
SP_SUMCOUNTS(i2,:) = SP_SUMCOUNTS(i2,:) + repmat(sumSp,[npops-1 1]);

new_i2_logml = computePopulationLogml(i2, adjprior_cq, adjprior_sp);

CQ_COUNTS(:,:,i2) = CQ_COUNTS(:,:,i2)-repmat(indCqCounts, [1 1 npops-1]);
CQ_SUMCOUNTS(i2,:) = CQ_SUMCOUNTS(i2,:)-repmat(sumCq,[npops-1 1]);
SP_COUNTS(:,:,i2) = SP_COUNTS(:,:,i2)-repmat(indSpCounts, [1 1 npops-1]);
SP_SUMCOUNTS(i2,:) = SP_SUMCOUNTS(i2,:) - repmat(sumSp,[npops-1 1]);
% a = repmat(sumSp,[npops-1 1]);

changes(i2) = new_i1_logml - i1_logml ...
    + new_i2_logml - i2_logml;

%------------------------------------------------------------------------------------

function changes = computeChanges2(i1, adjprior_cq, adjprior_sp)
% Computes changes in log marginal likelihood if population i1 is combined
% with another population
%
% Input:
% i1 - the population to be combined
% adjprior_cq & _sp  - adjpriors for cliques and separators
%
% Output:
% changes - table of size 1*npops. changes(i) = difference in logml if
% i1 is combined with population i.

global CQ_COUNTS;   global CQ_SUMCOUNTS;
global SP_COUNTS;   global SP_SUMCOUNTS;
global POP_LOGML;
npops = size(CQ_COUNTS,3);
changes = zeros(npops,1);

i1_logml = POP_LOGML(i1);
indCqCounts = CQ_COUNTS(:,:,i1);
indSpCounts = SP_COUNTS(:,:,i1);
sumCq = uint16(sum(indCqCounts,1));
sumSp = uint16(sum(indSpCounts,1));

new_i1_logml = 0;

i2 = [1:i1-1 , i1+1:npops];
i2_logml = POP_LOGML(i2);

CQ_COUNTS(:,:,i2) = CQ_COUNTS(:,:,i2)+repmat(indCqCounts, [1 1 npops-1]);
CQ_SUMCOUNTS(i2,:) = CQ_SUMCOUNTS(i2,:)+repmat(sumCq,[npops-1 1]);
SP_COUNTS(:,:,i2) = SP_COUNTS(:,:,i2)+repmat(indSpCounts, [1 1 npops-1]);
SP_SUMCOUNTS(i2,:) = SP_SUMCOUNTS(i2,:)+ repmat(sumSp,[npops-1 1]);
% a = repmat(sumSp,[npops-1 1]);
% if ~any(sumSp)
%     a(:,[1:size(a,2)])=[];
% end
% SP_SUMCOUNTS(i2,:) = SP_SUMCOUNTS(i2,:)+ a ;


new_i2_logml = computePopulationLogml(i2, adjprior_cq, adjprior_sp);

CQ_COUNTS(:,:,i2) = CQ_COUNTS(:,:,i2)-repmat(indCqCounts, [1 1 npops-1]);
CQ_SUMCOUNTS(i2,:) = CQ_SUMCOUNTS(i2,:)-repmat(sumCq,[npops-1 1]);
SP_COUNTS(:,:,i2) = SP_COUNTS(:,:,i2)-repmat(indSpCounts, [1 1 npops-1]);
SP_SUMCOUNTS(i2,:) = SP_SUMCOUNTS(i2,:)- repmat(sumSp,[npops-1 1]);

changes(i2) = new_i1_logml - i1_logml ...
    + new_i2_logml - i2_logml;




%------------------------------------------------------------------------------------


function changes = computeChanges3(T2, inds2, i1, counts_cq, counts_sp, ...
    adjprior_cq, adjprior_sp)
% Computes changes in log marginal likelihood if subpopulation of i2 is
% moved to another population
%
% Input:
% T2 - partition of inds2 to subpopulations
% inds2 - individuals in population i1
% i2
% counts_cq, counts_sp - counts for individuals
%
% Output:
% changes - table of size length(unique(T2))*npops.
% changes(i,j) = difference in logml if subpopulation inds2(find(T2==i)) of
% i2 is moved to population j

global CQ_COUNTS;   global CQ_SUMCOUNTS;
global SP_COUNTS;   global SP_SUMCOUNTS;
global POP_LOGML;
npops = size(CQ_COUNTS,3);
npops2 = length(unique(T2));
changes = zeros(npops2,npops);

%cq_counts = CQ_COUNTS;
%sp_counts = SP_COUNTS;
%cq_sumcounts = CQ_SUMCOUNTS;
%sp_sumcounts = SP_SUMCOUNTS;


i1_logml = POP_LOGML(i1);

for pop2 = 1:npops2
    % inds = inds2(find(T2==pop2));
    inds = inds2(logical(T2==pop2));
    ninds = length(inds);
    if ninds>0
        indCqCounts = uint16(sum(counts_cq(:,:,inds),3));
        indSpCounts = uint16(sum(counts_sp(:,:,inds),3));
        sumCq = uint16(sum(indCqCounts,1));
        sumSp = uint16(sum(indSpCounts,1));

        CQ_COUNTS(:,:,i1) = CQ_COUNTS(:,:,i1)-indCqCounts;
        CQ_SUMCOUNTS(i1,:) = CQ_SUMCOUNTS(i1,:)-sumCq;
        SP_COUNTS(:,:,i1) = SP_COUNTS(:,:,i1)-indSpCounts;
        SP_SUMCOUNTS(i1,:) = SP_SUMCOUNTS(i1,:)-sumSp;

        new_i1_logml = computePopulationLogml(i1, adjprior_cq, adjprior_sp);

        CQ_COUNTS(:,:,i1) = CQ_COUNTS(:,:,i1)+indCqCounts;
        CQ_SUMCOUNTS(i1,:) = CQ_SUMCOUNTS(i1,:)+sumCq;
        SP_COUNTS(:,:,i1) = SP_COUNTS(:,:,i1)+indSpCounts;
        SP_SUMCOUNTS(i1,:) = SP_SUMCOUNTS(i1,:)+sumSp;

        i2 = [1:i1-1 , i1+1:npops];
        i2_logml = POP_LOGML(i2)';

        CQ_COUNTS(:,:,i2) = CQ_COUNTS(:,:,i2)+repmat(indCqCounts, [1 1 npops-1]);
        CQ_SUMCOUNTS(i2,:) = CQ_SUMCOUNTS(i2,:)+repmat(sumCq,[npops-1 1]);
        SP_COUNTS(:,:,i2) = SP_COUNTS(:,:,i2)+repmat(indSpCounts, [1 1 npops-1]);
        SP_SUMCOUNTS(i2,:) = SP_SUMCOUNTS(i2,:)+ repmat(sumSp,[npops-1 1]);
    
        new_i2_logml = computePopulationLogml(i2, adjprior_cq, adjprior_sp)';

        CQ_COUNTS(:,:,i2) = CQ_COUNTS(:,:,i2)-repmat(indCqCounts, [1 1 npops-1]);
        CQ_SUMCOUNTS(i2,:) = CQ_SUMCOUNTS(i2,:)-repmat(sumCq,[npops-1 1]);
        SP_COUNTS(:,:,i2) = SP_COUNTS(:,:,i2)-repmat(indSpCounts, [1 1 npops-1]);
        SP_SUMCOUNTS(i2,:) = SP_SUMCOUNTS(i2,:)- repmat(sumSp,[npops-1 1]);

        changes(pop2,i2) = new_i1_logml - i1_logml ...
            + new_i2_logml - i2_logml;
    end
end

%--------------------------------------------------------------------------

function changes = computeChanges5(inds, i1, i2, counts_cq, counts_sp, ...
    adjprior_cq, adjprior_sp)
% Computes change in logml if individual of inds is moved between
% populations i1 and i2

global CQ_COUNTS;   global CQ_SUMCOUNTS;
global SP_COUNTS;   global SP_SUMCOUNTS;
global POP_LOGML;   global PARTITION;

ninds = length(inds);
changes = zeros(ninds,1);

i1_logml = POP_LOGML(i1);
i2_logml = POP_LOGML(i2);

for i = 1:ninds
    ind = inds(i);
    if PARTITION(ind)==i1
        pop1 = i1;  %from
        pop2 = i2;  %to
    else
        pop1 = i2;
        pop2 = i1;
    end
    indCqCounts = uint16(counts_cq(:,:,ind));
    indSpCounts = uint16(counts_sp(:,:,ind));
    sumCq = uint16(sum(indCqCounts,1));
    sumSp = uint16(sum(indSpCounts,1));

    CQ_COUNTS(:,:,pop1) = CQ_COUNTS(:,:,pop1)-indCqCounts;
    CQ_SUMCOUNTS(pop1,:) = CQ_SUMCOUNTS(pop1,:)-sumCq;
    SP_COUNTS(:,:,pop1) = SP_COUNTS(:,:,pop1)-indSpCounts;
    SP_SUMCOUNTS(pop1,:) = SP_SUMCOUNTS(pop1,:) - sumSp;

    CQ_COUNTS(:,:,pop2) = CQ_COUNTS(:,:,pop2)+indCqCounts;
    CQ_SUMCOUNTS(pop2,:) = CQ_SUMCOUNTS(pop2,:)+sumCq;
    SP_COUNTS(:,:,pop2) = SP_COUNTS(:,:,pop2)+indSpCounts;
    SP_SUMCOUNTS(pop2,:) = SP_SUMCOUNTS(pop2,:) + sumSp;

    new_logmls = computePopulationLogml([i1 i2], adjprior_cq, adjprior_sp);
    changes(i) = sum(new_logmls);

    CQ_COUNTS(:,:,pop1) = CQ_COUNTS(:,:,pop1)+indCqCounts;
    CQ_SUMCOUNTS(pop1,:) = CQ_SUMCOUNTS(pop1,:)+sumCq;
    SP_COUNTS(:,:,pop1) = SP_COUNTS(:,:,pop1)+indSpCounts;
    SP_SUMCOUNTS(pop1,:) = SP_SUMCOUNTS(pop1,:)+sumSp;
    CQ_COUNTS(:,:,pop2) = CQ_COUNTS(:,:,pop2)-indCqCounts;
    CQ_SUMCOUNTS(pop2,:) = CQ_SUMCOUNTS(pop2,:)-sumCq;
    SP_COUNTS(:,:,pop2) = SP_COUNTS(:,:,pop2)-indSpCounts;
    SP_SUMCOUNTS(pop2,:) = SP_SUMCOUNTS(pop2,:)-sumSp;
end

changes = changes - i1_logml - i2_logml;


%-------------------------------------------------------------------------------------


function updateGlobalVariables(ind, i2, indCqCounts, indSpCounts, ...
    adjprior_cq, adjprior_sp)
% Updates global variables when individual ind is moved to population i2

global CQ_COUNTS;   global CQ_SUMCOUNTS;
global SP_COUNTS;   global SP_SUMCOUNTS;
global PARTITION;   global POP_LOGML;

i1 = PARTITION(ind);
PARTITION(ind)=i2;

sumCq = uint16(sum(indCqCounts,1));
sumSp = uint16(sum(indSpCounts,1));

CQ_COUNTS(:,:,i1) = CQ_COUNTS(:,:,i1)-indCqCounts;
CQ_SUMCOUNTS(i1,:) = CQ_SUMCOUNTS(i1,:)-sumCq;
SP_COUNTS(:,:,i1) = SP_COUNTS(:,:,i1)-indSpCounts;
SP_SUMCOUNTS(i1,:) = SP_SUMCOUNTS(i1,:)-sumSp;

CQ_COUNTS(:,:,i2) = CQ_COUNTS(:,:,i2)+indCqCounts;
CQ_SUMCOUNTS(i2,:) = CQ_SUMCOUNTS(i2,:)+sumCq;
SP_COUNTS(:,:,i2) = SP_COUNTS(:,:,i2)+indSpCounts;
SP_SUMCOUNTS(i2,:) = SP_SUMCOUNTS(i2,:)+sumSp;


POP_LOGML([i1 i2]) = computePopulationLogml([i1 i2], adjprior_cq, adjprior_sp);


%---------------------------------------------------------------------------------


function updateGlobalVariables2(i1, i2, adjprior_cq, adjprior_sp)
% Updates global variables when all individuals from population i1 are moved
% to population i2

global CQ_COUNTS;   global CQ_SUMCOUNTS;
global SP_COUNTS;   global SP_SUMCOUNTS;
global PARTITION;   global POP_LOGML;

% inds = find(PARTITION==i1);
% PARTITION(inds) = i2;
PARTITION(logical(PARTITION==i1)) = i2;

CQ_COUNTS(:,:,i2) = CQ_COUNTS(:,:,i2)+CQ_COUNTS(:,:,i1);
CQ_SUMCOUNTS(i2,:) = CQ_SUMCOUNTS(i2,:)+CQ_SUMCOUNTS(i1,:);
SP_COUNTS(:,:,i2) = SP_COUNTS(:,:,i2)+SP_COUNTS(:,:,i1);
SP_SUMCOUNTS(i2,:) = SP_SUMCOUNTS(i2,:)+SP_SUMCOUNTS(i1,:);

CQ_COUNTS(:,:,i1) = 0;
CQ_SUMCOUNTS(i1,:) = 0;
SP_COUNTS(:,:,i1) = 0;
SP_SUMCOUNTS(i1,:) = 0;

POP_LOGML(i1) = 0;
POP_LOGML(i2) = computePopulationLogml(i2, adjprior_cq, adjprior_sp);


%------------------------------------------------------------------------------------


function updateGlobalVariables3(muuttuvat, i2, indCqCounts, indSpCounts, ...
    adjprior_cq, adjprior_sp)
% Updates global variables when individuals muuttuvat are moved to
% population i2

global CQ_COUNTS;   global CQ_SUMCOUNTS;
global SP_COUNTS;   global SP_SUMCOUNTS;
global PARTITION;   global POP_LOGML;

i1 = PARTITION(muuttuvat(1));
PARTITION(muuttuvat) = i2;

sumCq = uint16(sum(indCqCounts,1));
sumSp = uint16(sum(indSpCounts,1));

CQ_COUNTS(:,:,i1) = CQ_COUNTS(:,:,i1)-indCqCounts;
CQ_SUMCOUNTS(i1,:) = CQ_SUMCOUNTS(i1,:)-sumCq;
SP_COUNTS(:,:,i1) = SP_COUNTS(:,:,i1)-indSpCounts;
SP_SUMCOUNTS(i1,:) = SP_SUMCOUNTS(i1,:)-sumSp;

CQ_COUNTS(:,:,i2) = CQ_COUNTS(:,:,i2)+indCqCounts;
CQ_SUMCOUNTS(i2,:) = CQ_SUMCOUNTS(i2,:)+sumCq;
SP_COUNTS(:,:,i2) = SP_COUNTS(:,:,i2)+indSpCounts;
SP_SUMCOUNTS(i2,:) = SP_SUMCOUNTS(i2,:)+sumSp;


POP_LOGML([i1 i2]) = computePopulationLogml([i1 i2], adjprior_cq, adjprior_sp);

%----------------------------------------------------------------------


function inds = returnInOrder(inds, pop, counts_cq, counts_sp, ...
    adjprior_cq, adjprior_sp)
% Returns individuals inds in order according to the change in the logml if
% they are moved out of the population pop

global CQ_COUNTS;   global CQ_SUMCOUNTS;
global SP_COUNTS;   global SP_SUMCOUNTS;

ninds = length(inds);
apuTaulu = [inds, zeros(ninds,1)];

for i=1:ninds
    ind = inds(i);
    indCqCounts = uint16(counts_cq(:,:,ind));
    indSpCounts = uint16(counts_sp(:,:,ind));
    sumCq = uint16(sum(indCqCounts,1));
    sumSp = uint16(sum(indSpCounts,1));

    CQ_COUNTS(:,:,pop) = CQ_COUNTS(:,:,pop)-indCqCounts;
    CQ_SUMCOUNTS(pop,:) = CQ_SUMCOUNTS(pop,:)-sumCq;
    SP_COUNTS(:,:,pop) = SP_COUNTS(:,:,pop)-indSpCounts;
    SP_SUMCOUNTS(pop,:) = SP_SUMCOUNTS(pop,:)-sumSp;

    apuTaulu(i, 2) = computePopulationLogml(pop, adjprior_cq, adjprior_sp);

    CQ_COUNTS(:,:,pop) = CQ_COUNTS(:,:,pop)+indCqCounts;
    CQ_SUMCOUNTS(pop,:) = CQ_SUMCOUNTS(pop,:)+sumCq;
    SP_COUNTS(:,:,pop) = SP_COUNTS(:,:,pop)+indSpCounts;
    SP_SUMCOUNTS(pop,:) = SP_SUMCOUNTS(pop,:)+sumSp;
end
apuTaulu = sortrows(apuTaulu,2);
inds = apuTaulu(ninds:-1:1,1);


%-------------------------------------------------------------------------

function [Z, dist] = newGetDistances(data, rowsFromInd)

ninds = max(data(:,end));
nloci = size(data,2)-1;
riviLkm = nchoosek(double(ninds),2);

% empties = find(data<0);
% data(empties)=0;
data(logical(data<0)) = 0;
data = uint16(data);

pariTaulu = zeros(riviLkm,2);
aPointer=1;

for a=1:ninds-1
    pariTaulu(aPointer:aPointer+double(ninds-1-a),1) = ones(ninds-a,1,'uint16')*a;
    pariTaulu(aPointer:aPointer+double(ninds-1-a),2) = uint16((a+1:ninds)');
    aPointer = aPointer+double(ninds-a);
end

eka = pariTaulu(:,ones(1,rowsFromInd));
eka = eka * rowsFromInd;
miinus = repmat(rowsFromInd-1 : -1 : 0, [riviLkm 1]);
eka = eka - miinus;

toka = pariTaulu(:,ones(1,rowsFromInd)*2);
toka = toka * rowsFromInd;
toka = toka - miinus;

eka = uint16(eka);
toka = uint16(toka);

clear pariTaulu; clear miinus;

summa = uint16(zeros(riviLkm,1));
vertailuja = uint16(zeros(riviLkm,1));

x = zeros(size(eka));    x = uint16(x);
y = zeros(size(toka));   y = uint16(y);
% fprintf(1,'%%10');
for j=1:nloci;
    
    for k=1:rowsFromInd
        x(:,k) = data(eka(:,k),j);
        y(:,k) = data(toka(:,k),j);
    end

    for a=1:rowsFromInd
        for b=1:rowsFromInd
            vertailutNyt = uint16(x(:,a)>0 & y(:,b)>0);
            vertailuja = vertailuja + vertailutNyt;
            lisays = (x(:,a)~=y(:,b) & vertailutNyt);
            summa = summa + uint16(lisays);
        end
    end
    % fprintf(1,'\b\b');
    % fprintf(1,'%d',floor(10+80*j/nloci));
end

clear x;    clear y;   clear vertailutNyt;
clear eka; clear toka; clear data; clear lisays;
dist = zeros(length(vertailuja),1);
% nollat = find(vertailuja==0);
% dist(nollat) = 1;
dist(logical(vertailuja==0)) = 1;
muut = find(vertailuja>0);
dist(muut) = double(summa(muut))./double(vertailuja(muut));
clear summa; clear vertailuja; clear muut;

Z = computeLinkage(dist');
% fprintf(1,'\b\b');
% fprintf(1,'%d\n',100);
%--------------------------------------------------------------------------

function clearGlobalVars

global CQ_COUNTS; CQ_COUNTS = [];
global CQ_SUMCOUNTS; CQ_SUMCOUNTS = [];
global SP_COUNTS; SP_COUNTS = [];
global SP_SUMCOUNTS; SP_SUMCOUNTS = [];
global PARTITION; PARTITION = [];
global POP_LOGML; POP_LOGML = [];

%--------------------------------------------------------------------------

function npops = removeEmptyPops
% Removes empty pops from all global COUNTS variables.
% Updates PARTITION and npops

global CQ_COUNTS;
global CQ_SUMCOUNTS;
global SP_COUNTS;
global SP_SUMCOUNTS;
global PARTITION;

notEmpty = find(any(CQ_SUMCOUNTS,2));
CQ_COUNTS = CQ_COUNTS(:,:,notEmpty);
CQ_SUMCOUNTS = CQ_SUMCOUNTS(notEmpty,:);
SP_COUNTS = SP_COUNTS(:,:,notEmpty);
SP_SUMCOUNTS = SP_SUMCOUNTS(notEmpty,:);

for n=1:length(notEmpty)
%     apu = find(PARTITION==notEmpty(n));
%     PARTITION(apu)=n;
PARTITION(logical(PARTITION==notEmpty(n))) = n;
end
npops = length(notEmpty);

%--------------------------------------------------------------------------

function [partitionSummary, added] = addToSummary(logml, partitionSummary, worstIndex)
% Tiedet‰‰n, ett?annettu logml on isompi kuin huonoin arvo
% partitionSummary taulukossa. Jos partitionSummary:ss?ei viel?ole
% annettua logml arvoa, niin lis‰t‰‰n worstIndex:in kohtaan uusi logml ja
% nykyist?partitiota vastaava nclusters:in arvo. Muutoin ei tehd?mit‰‰n.
global PARTITION;
apu = isempty(find(abs(partitionSummary(:,2)-logml)<1e-5,1));
if apu
    % Nyt lˆydetty partitio ei ole viel?kirjattuna summaryyn.
    npops = length(unique(PARTITION));
    partitionSummary(worstIndex,1) = npops;
    partitionSummary(worstIndex,2) = logml;
    added = 1;
else
    added = 0;
end

%--------------------------------------------------------------------------
function [counts, sumcounts] = initialCounts(ind_counts, PARTITION)
% This version is for partition compare mode

pops = unique(PARTITION);
npops = max(pops);

counts = zeros(size(ind_counts,1), size(ind_counts,2), npops,'uint16');
sumcounts = zeros(npops, size(ind_counts,2),'uint16');

for i = 1:npops
    inds = find(PARTITION == i);
    counts(:,:,i) = sum(ind_counts(:,:,inds), 3);
    sumcounts(i,:) = sum(counts(:,:,i),1);
end
%--------------------------------------------------------------------------


function logml = computeLogml(adjprior_cq, adjprior_sp, ...
                              CQ_COUNTS, CQ_SUMCOUNTS, ...
                              SP_COUNTS, SP_SUMCOUNTS)

                          % for partition compare mode
cq_counts = double(CQ_COUNTS);
cq_sumcounts = double(CQ_SUMCOUNTS);
sp_counts = double(SP_COUNTS);
sp_sumcounts = double(SP_SUMCOUNTS);

npops = size(CQ_COUNTS, 3);

cq_logml = sum(sum(sum(gammaln(cq_counts+repmat(adjprior_cq,[1 1 npops]))))) ...
    - npops*sum(sum(gammaln(adjprior_cq))) - ...
    sum(sum(gammaln(1+cq_sumcounts)));

sp_logml = sum(sum(sum(gammaln(sp_counts+repmat(adjprior_sp,[1 1 npops]))))) ...
    - npops*sum(sum(gammaln(adjprior_sp))) - ...
    sum(sum(gammaln(1+sp_sumcounts)));

logml = cq_logml - sp_logml;
clear cq_counts cq_sumcounts sp_counts sp_sumcounts;

%--------------------------------------------------------------------------

function popLogml = computePopulationLogml(pops, adjprior_cq, adjprior_sp)
% Palauttaa length(pops)*1 taulukon, jossa on laskettu korikohtaiset
% logml:t koreille, jotka on m‰‰ritelty pops-muuttujalla.

global CQ_COUNTS;   global CQ_SUMCOUNTS;
global SP_COUNTS;   global SP_SUMCOUNTS;

cq_counts = double(CQ_COUNTS);
cq_sumcounts = double(CQ_SUMCOUNTS);
sp_counts = double(SP_COUNTS);
sp_sumcounts = double(SP_SUMCOUNTS);

nall_cq = size(CQ_COUNTS,1);
nall_sp = size(SP_COUNTS, 1);
ncliq = size(CQ_COUNTS,2);
nsep = size(SP_COUNTS, 2);

z = length(pops);

popLogml_cq = ...
    squeeze(sum(sum(reshape(...
    gammaln(repmat(adjprior_cq,[1 1 z]) + cq_counts(:,:,pops)) ...
    ,[nall_cq ncliq z]),1),2)) - sum(gammaln(1+cq_sumcounts(pops,:)),2) - ...
    sum(sum(gammaln(adjprior_cq)));

popLogml_sp = ...
    squeeze(sum(sum(reshape(...
    gammaln(repmat(adjprior_sp,[1 1 z]) + sp_counts(:,:,pops)) ...
    ,[nall_sp nsep z]),1),2)) - sum(gammaln(1+sp_sumcounts(pops,:)),2) - ...
    sum(sum(gammaln(adjprior_sp)));

popLogml = popLogml_cq - popLogml_sp;
clear cq_counts cq_sumcounts sp_counts sp_sumcounts;

%-------------------------------------------------------------------


function changesInLogml = writeMixtureInfo(logml, counts_cq, counts_sp, adjprior_cq, ...
    adjprior_sp, outPutFile, inputFile, partitionSummary, popnames, linkage_model,...
    fixedK)

global PARTITION;
global CQ_COUNTS; 

%global CQ_SUMCOUNTS;
%global SP_COUNTS;  global SP_SUMCOUNTS;
ninds = length(PARTITION);
npops =  size(CQ_COUNTS,3);
names = (size(popnames,1) == ninds);    %Tarkistetaan ett?nimet viittaavat yksilˆihin

if length(outPutFile)>0
    fid = fopen(outPutFile,'a');
else
    fid = -1;
    diary('baps4_output.baps'); % save in text anyway.
end

dispLine;
disp('RESULTS OF INDIVIDUAL LEVEL MIXTURE ANALYSIS:');
disp(['Data file/ Linkage map: ' inputFile]);
disp(['Model: Codon']);
disp(['Number of clustered individuals: ' ownNum2Str(ninds)]);
disp(['Number of groups in optimal partition: ' ownNum2Str(npops)]);
disp(['Log(marginal likelihood) of optimal partition: ' ownNum2Str(logml)]);
disp(' ');
if (fid ~= -1)
    fprintf(fid,'%s \n', ['RESULTS OF INDIVIDUAL LEVEL MIXTURE ANALYSIS:']); fprintf(fid,'\n');
    fprintf(fid,'%s \n', ['Data file: ' inputFile]); fprintf(fid,'\n');
    fprintf(fid,'%s \n', ['Number of clustered individuals: ' ownNum2Str(ninds)]); fprintf(fid,'\n');
    fprintf(fid,'%s \n', ['Number of groups in optimal partition: ' ownNum2Str(npops)]); fprintf(fid,'\n');
    fprintf(fid,'%s \n', ['Log(marginal likelihood) of optimal partition: ' ownNum2Str(logml)]); fprintf(fid,'\n');
end

cluster_count = length(unique(PARTITION));
disp('Best Partition: ');
if (fid ~= -1)
    fprintf(fid,'%s \n','Best Partition: '); fprintf(fid,'\n');
end
for m=1:cluster_count
    indsInM = find(PARTITION==m);
    length_of_beginning = 11 + floor(log10(m));
    cluster_size = length(indsInM);

    if names
        text = ['Cluster ' num2str(m) ': {' char(popnames{indsInM(1)})];
        for k = 2:cluster_size
            text = [text ', ' char(popnames{indsInM(k)})];
        end;
    else
        text = ['Cluster ' num2str(m) ': {' num2str(indsInM(1))];
        for k = 2:cluster_size
            text = [text ', ' num2str(indsInM(k))];
        end;
    end
    text = [text '}'];
    while length(text)>58
        %Take one line and display it.
        new_line = takeLine(text,58);
        text = text(length(new_line)+1:end);
        disp(new_line);
        if (fid ~= -1)
            fprintf(fid,'%s \n',new_line);
            fprintf(fid,'\n');
        end
        if length(text)>0
            text = [blanks(length_of_beginning) text];
        else
            text = [];
        end;
    end;
    if ~isempty(text)
        disp(text);
        if (fid ~= -1)
            fprintf(fid,'%s \n',text);
            fprintf(fid,'\n');
        end
    end;
end

if npops == 1
    changesInLogml = [];
else
    disp(' ');
    disp(' ');
    disp('Changes in log(marginal likelihood) if indvidual i is moved to group j:');
    if (fid ~= -1)
        fprintf(fid, '%s \n', ' '); fprintf(fid, '\n');
        fprintf(fid, '%s \n', ' '); fprintf(fid, '\n');
        fprintf(fid, '%s \n', 'Changes in log(marginal likelihood) if indvidual i is moved to group j:'); fprintf(fid, '\n');
    end

    if names
        nameSizes = zeros(ninds,1);
        for i = 1:ninds
            nimi = char(popnames{i});
            nameSizes(i) = length(nimi);
        end
        maxSize = max(nameSizes);
        maxSize = max(maxSize, 5);
        erotus = maxSize - 5;
        alku = blanks(erotus);
        ekarivi = [alku '  ind' blanks(6+erotus)];
    else
        ekarivi = '  ind      ';
    end

    for i = 1:cluster_count
        ekarivi = [ekarivi ownNum2Str(i) blanks(8-floor(log10(i)))];
    end
    disp(ekarivi);
    if (fid ~= -1)
        fprintf(fid, '%s \n', ekarivi); fprintf(fid, '\n');
    end

    %ninds = size(data,1)/rowsFromInd;
    changesInLogml = zeros(npops,ninds);
    for ind = 1:ninds
        indCqCounts = uint16(counts_cq(:,:,ind));
        indSpCounts = uint16(counts_sp(:,:,ind));
        changesInLogml(:,ind) = computeChanges(ind, adjprior_cq, ...
            adjprior_sp, indCqCounts, indSpCounts);
        
        % transform the logml change to conditional posterior probabilities
        % Added by Lu Cheng, 28.03.2010
        tmp_changelogml = exp(changesInLogml(:,ind));
        changesInLogml(:,ind) = tmp_changelogml ./ sum(tmp_changelogml);
        clear tmp_changelogml;
        %----------------------
        
        if names
            nimi = char(popnames{ind});
            rivi = [blanks(maxSize - length(nimi)) nimi ':'];
        else
            rivi = [blanks(4-floor(log10(ind))) ownNum2Str(ind) ':'];
        end
        for j = 1:npops
            rivi = [rivi '  ' logml2String(omaRound(changesInLogml(j,ind)))];
        end
        disp(rivi);
        if (fid ~= -1)
            fprintf(fid, '%s \n', rivi); fprintf(fid, '\n');
        end
    end


    %     % KL-divergence has to be calculated otherwise...
    %     % {
    %     disp(' '); disp(' ');
    %     disp('KL-divergence matrix:');
    %
    %     if (fid ~= -1)
    %         fprintf(fid, '%s \n', [' ']); fprintf(fid, '\n');
    %         fprintf(fid, '%s \n', [' ']); fprintf(fid, '\n');
    %         fprintf(fid, '%s \n', ['KL-divergence matrix:']); fprintf(fid, '\n');
    %     end
    %
    %     maxnoalle = size(COUNTS,1);
    %     nloci = size(COUNTS,2);
    %     d = zeros(maxnoalle, nloci, npops);
    %     prior = adjprior;
    %     prior(find(prior==1))=0;
    %     nollia = find(all(prior==0));  %Lokukset, joissa oli havaittu vain yht?alleelia.
    %     prior(1,nollia)=1;
    %     for pop1 = 1:npops
    %         d(:,:,pop1) = (squeeze(COUNTS(:,:,pop1))+prior) ./ repmat(sum(squeeze(COUNTS(:,:,pop1))+prior),maxnoalle,1);
    %         dist1(pop1) = (squeeze(COUNTS(:,:,pop1))+adjprior) ./ repmat((SUMCOUNTS(pop1,:)+adjprior), maxnoalle, 1);
    %     end
    %     ekarivi = blanks(7);
    %     for pop = 1:npops
    %         ekarivi = [ekarivi num2str(pop) blanks(7-floor(log10(pop)))];
    %     end
    %     disp(ekarivi);
    %     if (fid ~= -1)
    %         fprintf(fid, '%s \n', [ekarivi]); fprintf(fid, '\n');
    %     end
    %
    %     for pop1 = 1:npops
    %         rivi = [blanks(2-floor(log10(pop1))) num2str(pop1) '  '];
    %         for pop2 = 1:pop1-1
    %             dist1 = d(:,:,pop1); dist2 = d(:,:,pop2);
    %             div12 = sum(sum(dist1.*log2((dist1+10^-10) ./ (dist2+10^-10))))/nloci;
    %             div21 = sum(sum(dist2.*log2((dist2+10^-10) ./ (dist1+10^-10))))/nloci;
    %             div = (div12+div21)/2;
    %             rivi = [rivi kldiv2str(div) '  '];
    %         end
    %         disp(rivi);
    %         if (fid ~= -1)
    %             fprintf(fid, '%s \n', [rivi]); fprintf(fid, '\n');
    %         end
    %     end
    %     % }
end

disp(' ');
disp(' ');
disp('List of sizes of 10 best visited partitions and corresponding log(ml) values');

if (fid ~= -1)
    fprintf(fid, '%s \n', ' '); fprintf(fid, '\n');
    fprintf(fid, '%s \n', ' '); fprintf(fid, '\n');
    fprintf(fid, '%s \n', 'List of sizes of 10 best visited partitions and corresponding log(ml) values'); fprintf(fid, '\n');
end

partitionSummary = sortrows(partitionSummary,2);
partitionSummary = partitionSummary(size(partitionSummary,1):-1:1 , :);
% partitionSummary = partitionSummary(find(partitionSummary(:,2)>-1e49),:);
partitionSummary = partitionSummary(logical(partitionSummary(:,2)>-1e49),:);
if size(partitionSummary,1)>10
    vikaPartitio = 10;
else
    vikaPartitio = size(partitionSummary,1);
end
for part = 1:vikaPartitio
    line = [num2str(partitionSummary(part,1)) '    ' num2str(partitionSummary(part,2))];
    disp(line);
    if (fid ~= -1)
        fprintf(fid, '%s \n', line); fprintf(fid, '\n');
    end
end

if ~fixedK

    disp(' ');
    disp(' ');
    disp('Probabilities for number of clusters');

    if (fid ~= -1)
        fprintf(fid, '%s \n', ' '); fprintf(fid, '\n');
        fprintf(fid, '%s \n', ' '); fprintf(fid, '\n');
        fprintf(fid, '%s \n', 'Probabilities for number of clusters'); fprintf(fid, '\n');
    end

    npopsTaulu = unique(partitionSummary(:,1));
    len = length(npopsTaulu);
    probs = zeros(len,1);
    partitionSummary(:,2) = partitionSummary(:,2)-max(partitionSummary(:,2));
    sumtn = sum(exp(partitionSummary(:,2)));
    for i=1:len
        % npopstn = sum(exp(partitionSummary(find(partitionSummary(:,1)==npopsTaulu(i)),2)));
        npopstn = sum(exp(partitionSummary(logical(partitionSummary(:,1)==npopsTaulu(i)),2)));
        probs(i) = npopstn / sumtn;
    end
    for i=1:len
        if probs(i)>1e-5
            line = [num2str(npopsTaulu(i)) '   ' num2str(probs(i))];
            disp(line);
            if (fid ~= -1)
                fprintf(fid, '%s \n', line); fprintf(fid, '\n');
            end
        end
    end
    
end

if (fid ~= -1)
    fclose(fid);
else
    diary off
end

    
%--------------------------------------------------------------
function newline = takeLine(description,width)
%Returns one line from the description: line ends to the first
%space after width:th mark.
% newLine = description(1:width);
n = width+1;
while ~isspace(description(n)) && n<length(description)
    n = n+1;
end;
newline = description(1:n);


function dispLine
disp('---------------------------------------------------');

function dispCancel
disp('** CANCELLED');

function num2 = omaRound(num)
% Pyˆrist‰‰ luvun num 1 desimaalin tarkkuuteen
num = num*10;
num = round(num);
num2 = num/10;

%---------------------------------------------------------


function digit = palautaYks(num,yks)
% palauttaa luvun num 10^yks termin kertoimen
% string:in?
% yks t‰ytyy olla kokonaisluku, joka on
% v‰hint‰‰n -1:n suuruinen. Pienemmill?
% luvuilla tapahtuu jokin pyˆristysvirhe.

if yks>=0
    digit = rem(num, 10^(yks+1));
    digit = floor(digit/(10^yks));
else
    digit = num*10;
    digit = floor(rem(digit,10));
end
digit = num2str(digit);

%-------------------------------------------------------------------------

function [newData, rowsFromInd, alleleCodes, noalle, adjprior, priorTerm] = ...
    handleData(raw_data)
% Alkuper‰isen datan viimeinen sarake kertoo, milt?yksilˆlt?
% kyseinen rivi on per‰isin. Funktio tutkii ensin, ett?montako
% rivi?maksimissaan on per‰isin yhdelt?yksilˆlt? jolloin saadaan
% tiet‰‰ onko kyseess?haploidi, diploidi jne... T‰m‰n j‰lkeen funktio
% lis‰‰ tyhji?rivej?niille yksilˆille, joilta on per‰isin v‰hemm‰n
% rivej?kuin maksimim‰‰r?
%   Mik‰li jonkin alleelin koodi on =0, funktio muuttaa t‰m‰n alleelin
% koodi pienimm‰ksi koodiksi, joka isompi kuin mik‰‰n k‰ytˆss?oleva koodi.
% T‰m‰n j‰lkeen funktio muuttaa alleelikoodit siten, ett?yhden lokuksen j
% koodit saavat arvoja v‰lill?1,...,noalle(j).

data = raw_data;
nloci=size(raw_data,2)-1;

dataApu = data(:,1:nloci);
nollat = find(dataApu==0);
if ~isempty(nollat)
    isoinAlleeli = max(max(dataApu));
    dataApu(nollat) = isoinAlleeli+1;
    data(:,1:nloci) = dataApu;
end
% dataApu = []; 
% nollat = []; 
% isoinAlleeli = [];

noalle=zeros(1,nloci);
alleelitLokuksessa = cell(nloci,1);
for i=1:nloci
    alleelitLokuksessaI = unique(data(:,i));
   %alleelitLokuksessa{i,1} = alleelitLokuksessaI(find(alleelitLokuksessaI>=0));
    alleelitLokuksessa{i,1} = alleelitLokuksessaI(logical(alleelitLokuksessaI>=0));
    noalle(i) = length(alleelitLokuksessa{i,1});
end
alleleCodes = zeros(max(noalle),nloci);
for i=1:nloci
    alleelitLokuksessaI = alleelitLokuksessa{i,1};
    puuttuvia = max(noalle)-length(alleelitLokuksessaI);
    alleleCodes(:,i) = [alleelitLokuksessaI; zeros(puuttuvia,1)];
end

for loc = 1:nloci
    for all = 1:noalle(loc)
        % data(find(data(:,loc)==alleleCodes(all,loc)), loc)=all;
        data(logical(data(:,loc)==alleleCodes(all,loc)), loc)=all;
    end;
end;

nind = max(data(:,end));
nrows = size(data,1);
ncols = size(data,2);
rowsFromInd = zeros(nind,1);
for i=1:nind
    rowsFromInd(i) = length(find(data(:,end)==i));
end
maxRowsFromInd = max(rowsFromInd);
a = -999;
emptyRow = repmat(a, 1, ncols);
lessThanMax = find(rowsFromInd < maxRowsFromInd);
missingRows = maxRowsFromInd*nind - nrows;
data = [data; zeros(missingRows, ncols)];
pointer = 1;
for ind=lessThanMax'    %K‰y l‰pi ne yksilˆt, joilta puuttuu rivej?
    miss = maxRowsFromInd-rowsFromInd(ind);  % T‰lt?yksilˆlt?puuttuvien lkm.
    for j=1:miss
        rowToBeAdded = emptyRow;
        rowToBeAdded(end) = ind;
        data(nrows+pointer, :) = rowToBeAdded;
        pointer = pointer+1;
    end
end
data = sortrows(data, ncols);   % Sorttaa yksilˆiden mukaisesti
newData = data;
rowsFromInd = maxRowsFromInd;

adjprior = zeros(max(noalle),nloci);
priorTerm = 0;
for j=1:nloci
    adjprior(:,j) = [repmat(1/noalle(j), [noalle(j),1]) ; ones(max(noalle)-noalle(j),1)];
    priorTerm = priorTerm + noalle(j)*gammaln(1/noalle(j));
end

%--------------------------------------------------------------------------

function [emptyPop, pops] = findEmptyPop(npops)
% Palauttaa ensimm‰isen tyhj‰n populaation indeksin. Jos tyhji?
% populaatioita ei ole, palauttaa -1:n.

global PARTITION;
pops = unique(PARTITION)';
if (length(pops) ==npops)
    emptyPop = -1;
else
    popDiff = diff([0 pops npops+1]);
    emptyPop = min(find(popDiff > 1));
end

%--------------------------------------------------------------------------

function popnames = fixPopnames(popnames, ninds)

if length(popnames) == ninds
    for i=1:ninds
        if isnumeric(popnames{i})
            popnames{i} = num2str(popnames{i});
            % popnames(i) = num2str(popnames{i});
        end
        popnames{i} = cellstr(popnames{i});
        % popnames(i) = cellstr(popnames{i});
    end
end

%--------------------------------------------------------------------------
function isRational = isTheLoadedNewDataRational(data)
% The last column of the data must include numbers 1-npops
% If so, isRational = 1, otherwise isRational = 0.
% The row numbers must be larger than 1.
if size(data,1) == 1 
    isRational = 0;
    display('*** ERROR: Sample size must be larger than one');
    return;
end
last_column = data(:,end);
last_column = sort(last_column);
current = 1;
if last_column(1) ~= current
    isRational = 0;
    display('*** ERROR: Wrong Indexes in the data');
    return;
end;
lengthcol = length(last_column);
for n = 2:lengthcol
    if ~(last_column(n) == current || last_column(n) == current + 1)
        %Some population is missing from the last column
        isRational = 0;
        display('*** ERROR: Missing indexes in the data');
        return;
    end;
    current = last_column(n);
end;
isRational = 1;


% %-------------------------------------------------------------------------
% function isRational = isTheLoadedNewLinkageRational(linkage_data)
% % Each positive element must be unique.
% % If so, isRational = 1, otherwise isRational = 0;
% nonzero = find(linkage_data~=0);
% dif = diff(linkage_data(nonzero));
% if ~all(dif)
%     isRational = 0; return;
% end;
% isRational = 1;

%--------------------------------------------------------------------------

function [sumcounts, counts] = ...
    indLociCounts(partition, data, npops, noalle)

nloci=size(data,2)-1;
% ninds = size(data,1);

counts = zeros(max(noalle),nloci,npops);
sumcounts = zeros(npops,nloci);
for i=1:npops
    for j=1:nloci
        % havainnotLokuksessa = find(partition==i & data(:,j)>=0);
        havainnotLokuksessa = find(ismember(data(:,end),find(partition==i)));
        sumcounts(i,j) = length(havainnotLokuksessa);
        for k=1:noalle(j)
            alleleCode = k;
            N_ijk = length(find(data(havainnotLokuksessa,j)==alleleCode));
            counts(k,j,i) = N_ijk;
        end
    end
end

%-----------------------------------------------------------------------------------


function popnames = initPopNames(nameFile, indexFile)
%Palauttaa tyhj‰n, mik‰li nimitiedosto ja indeksitiedosto
% eiv‰t olleet yht?pitki?

popnames = [];
try
indices = load(indexFile);
catch
    msgbox('Loading of the index file was unsuccessful', ...
        'Error', 'error');
    return
end
fid = fopen(nameFile);
if fid == -1
    % File does not exist
    msgbox('Loading of the name file was unsuccessful', ...
        'Error', 'error');
    return;
end

line = fgetl(fid);
counter = 1;

while sum(line~=-1) && ~isempty(line)
    names{counter} = line;
    line = fgetl(fid);
    counter = counter + 1;
end;
fclose(fid);

if length(names) ~= length(indices)
    disp('The number of population names must be equal to the number of ');
    disp('entries in the file specifying indices of the first individuals of ');
    disp('each population.');
    return;
end

popnames = cell(length(names), 2);
for i = 1:length(names)
    popnames{i,1} = names(i);
    popnames{i,2} = indices(i);
end
