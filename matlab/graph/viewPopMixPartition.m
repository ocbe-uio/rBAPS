function viewPopMixPartition(partition, rows, popnames)
%"Yksilön" leveys on verrannollinen yksilön rivien lukumäärään.

notEmptyPops = length(unique(partition));
disp(['Number of populations: ' num2str(notEmptyPops)]);
if notEmptyPops>30
    disp(' ');
    disp('Figure can be drawn only if the number of populations');
    disp('is less or equal to 30.');
    disp(' ');
    return;
end

nind = length(partition);
totalNumRows = 0;
for ind = 1:nind
    totalNumRows = totalNumRows+rows(ind,2)-rows(ind,1)+1;
    %totalNumRows = totalNumRows+length(rows{ind});
end
%npops = max(partition);
npops = notEmptyPops;

varit = giveColors(npops);
korkeinviiva = 1.05;
pieninarvo = -korkeinviiva;

h0 = figure('NumberTitle', 'off'); %image_figure;   %Muutettu
tiedot.rows = rows;
tiedot.info = partition;
tiedot.popnames = popnames;
set(h0,'UserData',tiedot);

set(gca, 'Xlim', [-.5 ,totalNumRows+.5], 'YLim', [pieninarvo ,korkeinviiva], ...
    'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', []);

eiTyhjatPopulaatiot = unique(partition);

for i=1:nind
    pop = partition(i);
    pop = find(eiTyhjatPopulaatiot==pop);
           
    % Pylväiden piirtäminen
    for rivi = rows(i,1):rows(i,2)
        h0 =patch([rivi-1, rivi, rivi, rivi-1], [0, 0, 1, 1], varit(pop,:));
        set(h0,'EdgeColor','none'); % Midevaa varten kommentoitava!
    end
end


npops = size(rows,1);
for i=1:npops
    firstRow = rows(i,1);
    line([firstRow-1, firstRow-1], [0,1], 'Color', 'k');  %Populaatioiden rajat
end

if ~isempty(popnames)
    for i=1:npops
        %rivi = rows(i};
        %x_paikka = (rivi(1)-1+rivi(end))/2;    
        x_paikka = (rows(i,1) - 1 + rows(i,2)) / 2;
        
        korkeuskerroin = pieninarvo / -0.2;
        suhdekerroin = npops/6;
        for letter_num = 1:length(popnames{i,1}{1})
            letter= popnames{i,1}{1}(letter_num);%alter .004|
            text(x_paikka+korjaus(letter)*suhdekerroin, ...
                0.0005*korkeuskerroin-0.02*letter_num*korkeuskerroin, ...
                letter, 'Interpreter','none');
        end
    end
end

line([totalNumRows,totalNumRows],[0,1],'Color','k');

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
