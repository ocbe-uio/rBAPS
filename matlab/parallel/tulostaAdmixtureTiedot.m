%--------------------------------------------------------------------------


function tulostaAdmixtureTiedot(proportions, uskottavuus, alaRaja, niter)
h0 = findobj('Tag','filename1_text');
inputf = get(h0,'String');
h0 = findobj('Tag','filename2_text');
outf = get(h0,'String'); clear h0;

if length(outf)>0
    fid = fopen(outf,'a');
else
    fid = -1;
    diary('baps4_output.baps'); % save in text anyway.
end

ninds = length(uskottavuus);
npops = size(proportions,2);
disp(' ');
dispLine;
disp('RESULTS OF ADMIXTURE ANALYSIS BASED');
disp('ON MIXTURE CLUSTERING OF INDIVIDUALS');
disp(['Data file: ' inputf]);
disp(['Number of individuals: ' num2str(ninds)]);
disp(['Results based on ' num2str(niter) ' simulations from posterior allele frequencies.']);
disp(' ');
if fid ~= -1
    fprintf(fid, '\n');
    fprintf(fid,'%s \n', ['--------------------------------------------']); fprintf(fid, '\n');
    fprintf(fid,'%s \n', ['RESULTS OF ADMIXTURE ANALYSIS BASED']); fprintf(fid, '\n');
    fprintf(fid,'%s \n', ['ON MIXTURE CLUSTERING OF INDIVIDUALS']); fprintf(fid, '\n');
    fprintf(fid,'%s \n', ['Data file: ' inputf]); fprintf(fid, '\n');
    fprintf(fid,'%s \n', ['Number of individuals: ' num2str(ninds)]); fprintf(fid, '\n');
    fprintf(fid,'%s \n', ['Results based on ' num2str(niter) ' simulations from posterior allele frequencies.']); fprintf(fid, '\n');
    fprintf(fid, '\n');
end

ekaRivi = blanks(6);
for pop = 1:npops
    ekaRivi = [ekaRivi blanks(3-floor(log10(pop))) num2str(pop) blanks(2)];
end
ekaRivi = [ekaRivi blanks(1) 'p']; % Added on 29.08.06
disp(ekaRivi);
for ind = 1:ninds
    rivi = [num2str(ind) ':' blanks(4-floor(log10(ind)))];
    if any(proportions(ind,:)>0)
        for pop = 1:npops-1
            rivi = [rivi proportion2str(proportions(ind,pop)) blanks(2)];
        end
        rivi = [rivi proportion2str(proportions(ind,npops)) ':  '];
        rivi = [rivi ownNum2Str(uskottavuus(ind))];
    end
    disp(rivi);
    if fid ~= -1
        fprintf(fid,'%s \n',[rivi]); fprintf(fid,'\n');
    end
end
if fid ~= -1
    fclose(fid);
else
    diary off
end

