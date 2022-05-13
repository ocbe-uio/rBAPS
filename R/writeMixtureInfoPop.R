function changesInLogml = writeMixtureInfoPop(logml, rows, data, adjprior, ...
  priorTerm, outPutFile, inputFile, partitionSummary, popnames, fixedK)

global PARTITION;
global COUNTS;
global SUMCOUNTS;
global LOGDIFF;
ninds = size(rows,1);
npops =  size(COUNTS,3);
names = (size(popnames,1) == ninds);    %Tarkistetaan ett?nimet viittaavat yksilöihin
changesInLogml = [];
if length(outPutFile)>0
  fid = fopen(outPutFile,'a');
else
  fid = -1;
  diary('baps4_output.baps'); % save in text anyway.
end

dispLine;
disp('RESULTS OF GROUP LEVEL MIXTURE ANALYSIS:');
disp(['Data file: ' inputFile]);
disp(['Number of clustered groups: ' ownNum2Str(ninds)]);
disp(['Number of clusters in optimal partition: ' ownNum2Str(npops)]);
disp(['Log(marginal likelihood) of optimal partition: ' ownNum2Str(logml)]);
disp(' ');
if (fid ~= -1)
  fprintf(fid,'%s \n', ['RESULTS OF GROUP LEVEL MIXTURE ANALYSIS:']); fprintf(fid,'\n');
  fprintf(fid,'%s \n', ['Data file: ' inputFile]); fprintf(fid,'\n');
  fprintf(fid,'%s \n', ['Number of clustered groups: ' ownNum2Str(ninds)]); fprintf(fid,'\n');
  fprintf(fid,'%s \n', ['Number of clusters in optimal partition: ' ownNum2Str(npops)]); fprintf(fid,'\n');
  fprintf(fid,'%s \n', ['Log(marginal likelihood) of optimal partition: ' ownNum2Str(logml)]); fprintf(fid,'\n');
  fprintf(fid,'\n');
end

cluster_count = length(unique(PARTITION));
disp(['Best Partition: ']);
if (fid ~= -1)
  fprintf(fid,'%s \n',['Best Partition: ']); fprintf(fid,'\n');
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
          fprintf(fid,'%s \n',[new_line]);
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
          fprintf(fid,'%s \n',[text]);
          fprintf(fid,'\n');
      end
  end;
end

if npops > 1
  disp(' ');
  disp(' ');
  disp('Changes in log(marginal likelihood) if group i is moved to cluster j:');
  if (fid ~= -1)
      fprintf(fid, '%s \n', [' ']); %fprintf(fid, '\n');
      fprintf(fid, '%s \n', [' ']); %fprintf(fid, '\n');
      fprintf(fid, '%s \n', ['Changes in log(marginal likelihood) if group i is moved to cluster j:']); %fprintf(fid, '\n');
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
      ekarivi = [alku 'group' blanks(6+erotus)];
  else
      ekarivi = 'group      ';
  end
  for i = 1:cluster_count
      ekarivi = [ekarivi ownNum2Str(i) blanks(8-floor(log10(i)))];
  end
  disp(ekarivi);
  if (fid ~= -1)
      fprintf(fid, '%s \n', [ekarivi]); %fprintf(fid, '\n');
  end

  changesInLogml = LOGDIFF';
  for ind = 1:ninds
      %[muutokset, diffInCounts] = laskeMuutokset(ind, rows, data, ...
      %    adjprior, priorTerm);
      %changesInLogml(:,ind) = muutokset;
      muutokset = changesInLogml(:,ind);
      if names
          nimi = char(popnames{ind});
          rivi = [blanks(maxSize - length(nimi)) nimi ':'];
      else
          rivi = [blanks(4-floor(log10(ind))) ownNum2Str(ind) ':'];
      end
      for j = 1:npops
          rivi = [rivi '  ' logml2String(omaRound(muutokset(j)))];
      end
      disp(rivi);
      if (fid ~= -1)
          fprintf(fid, '%s \n', [rivi]); fprintf(fid, '\n');
      end
  end

  disp(' '); disp(' ');
  disp('KL-divergence matrix in PHYLIP format:');
  dist_mat = zeros(npops, npops);
  if (fid ~= -1)
      fprintf(fid, '%s \n', [' ']); %fprintf(fid, '\n');
      fprintf(fid, '%s \n', [' ']); %fprintf(fid, '\n');
      fprintf(fid, '%s \n', ['KL-divergence matrix in PHYLIP format:']); %fprintf(fid, '\n');
  end

  maxnoalle = size(COUNTS,1);
  nloci = size(COUNTS,2);
  d = zeros(maxnoalle, nloci, npops);
  prior = adjprior;
  prior(find(prior==1))=0;
  nollia = find(all(prior==0));  %Lokukset, joissa oli havaittu vain yht?alleelia.
  prior(1,nollia)=1;
  for pop1 = 1:npops
      d(:,:,pop1) = (squeeze(COUNTS(:,:,pop1))+prior) ./ repmat(sum(squeeze(COUNTS(:,:,pop1))+prior),maxnoalle,1);
      %dist1(pop1) = (squeeze(COUNTS(:,:,pop1))+adjprior) ./ repmat((SUMCOUNTS(pop1,:)+adjprior), maxnoalle, 1);
  end
%     ekarivi = blanks(7);
%     for pop = 1:npops
%         ekarivi = [ekarivi num2str(pop) blanks(7-floor(log10(pop)))];
%     end
  ekarivi = num2str(npops);
  disp(ekarivi);

  if (fid ~= -1)
      fprintf(fid, '%s \n', [ekarivi]); %fprintf(fid, '\n');
  end

  for pop1 = 1:npops
      rivi = [blanks(2-floor(log10(pop1))) num2str(pop1) '  '];
      for pop2 = 1:pop1-1
          dist1 = d(:,:,pop1); dist2 = d(:,:,pop2);
          div12 = sum(sum(dist1.*log2((dist1+10^-10) ./ (dist2+10^-10))))/nloci;
          div21 = sum(sum(dist2.*log2((dist2+10^-10) ./ (dist1+10^-10))))/nloci;
          div = (div12+div21)/2;
          % rivi = [rivi kldiv2str(div) '  '];
          dist_mat(pop1,pop2) = div;
      end
%         disp(rivi);
%         if (fid ~= -1)
%             fprintf(fid, '%s \n', [rivi]); fprintf(fid, '\n');
%         end
  end



  dist_mat = dist_mat + dist_mat'; % make it symmetric
  for pop1 = 1:npops
      rivi = ['Cluster_' num2str(pop1) ' '];
      for pop2 = 1:npops
          rivi = [rivi kldiv2str(dist_mat(pop1,pop2)) ' '];
      end
      disp(rivi);
      if (fid ~= -1)
          fprintf(fid, '%s \n', [rivi]); %fprintf(fid, '\n');
      end
  end

end

disp(' ');
disp(' ');
disp('List of sizes of 10 best visited partitions and corresponding log(ml) values');

if (fid ~= -1)
  fprintf(fid, '%s \n', [' ']); fprintf(fid, '\n');
  fprintf(fid, '%s \n', [' ']); fprintf(fid, '\n');
  fprintf(fid, '%s \n', ['List of sizes of 10 best visited partitions and corresponding log(ml) values']); fprintf(fid, '\n');
end

partitionSummary = sortrows(partitionSummary,2);
partitionSummary = partitionSummary(size(partitionSummary,1):-1:1 , :);
partitionSummary = partitionSummary(find(partitionSummary(:,2)>-1e49),:);
if size(partitionSummary,1)>10
  vikaPartitio = 10;
else
  vikaPartitio = size(partitionSummary,1);
end
for part = 1:vikaPartitio
  line = [num2str(partitionSummary(part,1)) '    ' num2str(partitionSummary(part,2))];
  disp(line);
  if (fid ~= -1)
      fprintf(fid, '%s \n', [line]); fprintf(fid, '\n');
  end
end

if ~fixedK
  disp(' ');
  disp(' ');
  disp('Probabilities for number of clusters');

  if (fid ~= -1)
      fprintf(fid, '%s \n', [' ']); fprintf(fid, '\n');
      fprintf(fid, '%s \n', [' ']); fprintf(fid, '\n');
      fprintf(fid, '%s \n', ['Probabilities for number of clusters']); fprintf(fid, '\n');
  end

  npopsTaulu = unique(partitionSummary(:,1));
  len = length(npopsTaulu);
  probs = zeros(len,1);
  partitionSummary(:,2) = partitionSummary(:,2)-max(partitionSummary(:,2));
  sumtn = sum(exp(partitionSummary(:,2)));
  for i=1:len
      npopstn = sum(exp(partitionSummary(find(partitionSummary(:,1)==npopsTaulu(i)),2)));
      probs(i) = npopstn / sumtn;
  end
  for i=1:len
      if probs(i)>1e-5
          line = [num2str(npopsTaulu(i)) '   ' num2str(probs(i))];
          disp(line);
          if (fid ~= -1)
              fprintf(fid, '%s \n', [line]); fprintf(fid, '\n');
          end
      end
  end
end

if (fid ~= -1)
  fclose(fid);
else
  diary off
end
