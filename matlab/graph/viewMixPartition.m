function viewMixPartition(partition,  popnames)

notEmptyPops = length(unique(partition));
if notEmptyPops>30
    disp(['Number of populations: ' num2str(notEmptyPops)]);
    disp(' ');
    disp('Figure can be drawn only if the number of populations');
    disp('is less or equal to 30.');
    disp(' ');
    return;
end

nind = length(partition);
%npops = max(partition);
npops = notEmptyPops;

varit = giveColors(npops);
korkeinviiva = 1.05;
pieninarvo = -korkeinviiva;

h0 = figure;
set(h0, 'NumberTitle', 'off'); %image_figure;   %Muutettu
tiedot.popnames = popnames;
tiedot.info = partition;
set(h0,'UserData',tiedot);

set(gca, 'Xlim', [-.5 ,nind+.5], 'YLim', [pieninarvo ,korkeinviiva], ...
    'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', []);

eiTyhjatPopulaatiot = unique(partition);

for i=1:nind
    % Suhteellisten osuuksien laskeminen
    pop = partition(i);
    pop = find(eiTyhjatPopulaatiot==pop);
           
    % Pylv‰‰n piirt‰minen
    h0 =patch([i-1, i, i, i-1], [0, 0, 1, 1], varit(pop,:));
    set(h0,'EdgeColor','none'); % Midevaa varten kommentoitava!   

end



if ~isempty(popnames)
    npops = size(popnames,1);
    for i=1:npops
        firstInd = popnames{i,2};
        if size(popnames,1) ~=nind
            line([firstInd-1, firstInd-1], [0,1], 'Color', 'k');  %Populaatioiden rajat
        end            
        if i<npops
            x_paikka = popnames{i,2}-1+(popnames{i+1,2}-popnames{i,2})/2;
        else
            x_paikka = popnames{i,2}-1+(nind+1-popnames{i,2})/2;
        end
               
        korkeuskerroin = pieninarvo / -0.2;
        suhdekerroin = npops/6;
        for letter_num = 1:length(popnames{i,1}{1})
            letter= popnames{i,1}{1}(letter_num);%alter .004|
            text(x_paikka+korjaus(letter)*suhdekerroin, ...
                0.0005*korkeuskerroin-0.02*letter_num*korkeuskerroin, ...
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
