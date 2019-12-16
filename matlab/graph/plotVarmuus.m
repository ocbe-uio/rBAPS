function plotVarmuus(V, C, pointers, varmuus, coordinates, partition, tekstit)

if nargin < 7
    tekstit = pointers;
end

notEmptyPops = length(unique(partition));
if notEmptyPops>30
    disp(['Number of populations: ' num2str(notEmptyPops)]);
    disp(' ');
    disp('Figure can be drawn only if the number of populations');
    disp('is less or equal to 30.');
    disp(' ');
    return;
end


h1 = figure;
hold on

colors=giveColors(notEmptyPops);

[I, J] = find(coordinates>0 | coordinates<0);
I=unique(I);
xmin = min(coordinates(I,1));
xmax = max(coordinates(I,1));
xdiff = (xmax-xmin);
xmean = xmin + xdiff/2;

ymin = min(coordinates(I,2));
ymax = max(coordinates(I,2));
ydiff = (ymax-ymin);
ymean = ymin + ydiff/2;

pituus = max(ydiff,xdiff)*1.1/2;

zmax = 0.8*max(varmuus) + 0.2;

axis([xmean-pituus xmean+pituus ymean-pituus ymean+pituus 0 zmax]);
grid(gca);
d = [1 2 3 1];



for i=1:length(C)
    koko = length(C{i});
    soluPisteet = V(C{i},:);
    center = [mean(soluPisteet(:,1)) mean(soluPisteet(:,2))];
    center = repmat(center, [koko 1]);
    soluPisteet = soluPisteet + (center - soluPisteet)./1000;
    
    
    apu = zeros(2*koko, 3);
    apu(1:koko, 1:2) = soluPisteet;
    apu(koko+1:end, 1:2) = soluPisteet;
    apu(koko+1:end, 3) = varmuus(i);
    
    taulu = pointers{i};
    if length(taulu)>0
        color = colors(partition(taulu(1)),:);
        pisteet =[1:koko 1];
        patch('XData', apu(pisteet,1), 'YData', apu(pisteet,2), ...
            'ZData', apu(pisteet,3), 'FaceColor',color, 'Clipping', ...
            'on', 'EdgeColor','k', 'LineWidth', 1);
        
        pisteet = pisteet+koko;
        patch('XData', apu(pisteet,1), 'YData', apu(pisteet,2), ...
            'ZData', apu(pisteet,3), 'FaceColor',color, 'Clipping', ...
            'on', 'EdgeColor','k', 'LineWidth', 1);
        
        for j = 1:koko-1
            pisteet = [j j+1 j+koko+1 j+koko j];
            patch('XData', apu(pisteet,1), 'YData', apu(pisteet,2), ...
                'ZData', apu(pisteet,3), 'FaceColor',color, 'Clipping', ...
                'on', 'EdgeColor','k', 'LineWidth', 1);
        end
        
        pisteet = [koko 1 koko+1 2*koko koko];
        patch('XData', apu(pisteet,1), 'YData', apu(pisteet,2), ...
            'ZData', apu(pisteet,3), 'FaceColor',color, 'Clipping', ...
            'on', 'EdgeColor','k', 'LineWidth', 1);
    end
end

if ~isequal(tekstit, -1)
    for i=1:length(pointers)
        taulu = pointers{i};
        teksti = tekstit{i};
        if isnumeric(teksti)
            teksti = num2str(teksti);
        end
        if length(taulu)>0
            text(coordinates(taulu(1),1),coordinates(taulu(1),2), ...
                varmuus(i) + zmax/100, teksti, 'FontSize', 10);
        end
    end
end

view(3);
hold off    


