function viewPartition4(osuudet, popnames, npops, admixnpops, ...
                        inliers, partition, filename)
% Kuin viewParition.m, mutta toimii itsenäisesti. Ei siis tarvitse
% globaalia COUNTS-muuttujaa.

nind = size(osuudet,1);
% inliers = unique(partition);
% inliers = find(any(osuudet));
if admixnpops>36
    disp(' ');
    disp('Figure can be drawn only if the number of populations');
    disp('is less or equal to 36.');
    disp(' ');
    return;
end

varit = giveColors(npops);
korkeinviiva = 1.05;
pieninarvo = -korkeinviiva;

h0 = figure;
set(h0, 'NumberTitle', 'off'); %image_figure;   %Muutettu
set(h0,'Tag','admixture_result');
set(h0,'Name',['Admixture Result - ' filename]); 
tiedot.popnames = popnames;
tiedot.info = osuudet;
set(h0,'UserData',tiedot);

set(gca, 'Xlim', [-.5 ,nind+.5], 'YLim', [pieninarvo ,korkeinviiva], ...
    'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', []);

for i=1:nind
    
    if any(osuudet(i,:)>0)
        cumOsuudet = cumsum(osuudet(i,:));
    
        % Pylvään piirtäminen
        for j = 1:admixnpops
            if j == 1
                if cumOsuudet(1)>0
                    h0 =patch([i-1, i, i, i-1], [0, 0, cumOsuudet(1), cumOsuudet(1)], varit(inliers(j),:));
                    set(h0,'EdgeColor','none'); % Midevaa varten kommentoitava!
                end
            else
                if (cumOsuudet(j)>cumOsuudet(j-1))
                    h0 = patch([i-1, i, i, i-1], [cumOsuudet(j-1), cumOsuudet(j-1), ...
                    cumOsuudet(j), cumOsuudet(j)], varit(inliers(j),:));
                    set(h0,'EdgeColor','none'); % Midevaa varten kommentoitava!
                end
            end
        end
    end
end



if ~isempty(popnames)
     npops = size(popnames,1); % NB! npops now stands for the number of individuals.
     for i = 1:npops
        firstInd = popnames{i,2};
        if size(popnames,1)~=nind
            line([firstInd-1, firstInd-1], [0,1], 'Color', 'k');  %Populaatioiden rajat
        end
        
        % The determination of x_paikka is changed, since popnames could be
        % non-sequential. - Jing
        if i<npops
            % x_paikka = popnames{i,2}-1+(popnames{i+1,2}-popnames{i,2})/2;
             x_paikka = i-1;
        else
            % x_paikka = popnames{i,2}-1+(nind+1-popnames{i,2})/2;
            x_paikka = i-1+(nind+1-i)/2;
        end
               
        % korkeuskerroin = pieninarvo / -0.2;
        suhdekerroin = npops/6;
        max_length = length(popnames{i,1}{1});
        for letter_num = 1: max_length 
            letter= popnames{i,1}{1}(letter_num); % alter .004|
%             text(x_paikka+korjaus(letter)*suhdekerroin, ...
%                 0.0005*korkeuskerroin-0.02*letter_num*korkeuskerroin, ...
%                 letter, 'Interpreter','none');
            text(x_paikka+korjaus(letter)*suhdekerroin, ...
                pieninarvo*letter_num/max_length, ...
                letter, 'Interpreter','none');
        end
    end
    line([nind,nind],[0,1],'Color','k');
end

%-------------------------------------------------------------------------------------

function extra = korjaus(letter)
if any(letter == 'ijlI')
    extra = 0.022;
elseif any(letter == 'r')
    extra = 0.016;
elseif any(letter == 'k')
    extra = 0.009;
elseif any(letter == 'f')
    extra = 0.013;
elseif any(letter == 't')
    extra = 0.014;
elseif any(letter == 'w')
    extra = -0.003;
else
    extra = 0;
end;