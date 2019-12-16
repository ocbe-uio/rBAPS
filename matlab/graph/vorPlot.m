function vorPlot(V,C,partition, pointers, coordinates, tekstit)

if nargin < 6
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

npops = length(unique(partition));

if npops > 30
    return
end

colors=giveColors(npops);

h1 = figure('NumberTitle', 'off', 'Name', 'Colored Voronoi tessellation');
hold on


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

axis([xmean-pituus xmean+pituus ymean-pituus ymean+pituus]);


for i=1:length(C)
    X=V(C{i},:);
    %k=convhull(X(:,1),X(:,2),{'QJ', 'Pp'});
    k=convhull(X(:,1),X(:,2));
    taulu = pointers{i};
    if length(taulu)>0
        color=colors(partition(taulu(1)),:);
        patch(X(k,1),X(k,2),color);
        plot(coordinates(taulu(1),1),coordinates(taulu(1),2),'Color',[1 1 1], 'MarkerSize', 50);
        %text(coordinates(taulu(1),1),coordinates(taulu(1),2),num2str(taulu), 'FontSize', 8);
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
            text(coordinates(taulu(1),1),coordinates(taulu(1),2),teksti, ...
                'Interpreter', 'none', 'FontSize', 8);
        end
    end
end
%[I,J] = find(coordinates(:,1) > 0);

%plot(coordinates(I,1), coordinates(I,2), 'k.');
    
%hold off
%{
%h0 = image_figure;
%hold on;
for i=1:length(partition)
    if coordinates(i,1)>=0
        %plot(coordinates(i,1),coordinates(i,2),'k.');
        plot(coordinates(i,1),coordinates(i,2),'Color',[1 1 1], 'MarkerSize', 40);
        text(coordinates(i,1),coordinates(i,2),num2str(i), 'FontSize', 8);
    end
end
%}
hold off;