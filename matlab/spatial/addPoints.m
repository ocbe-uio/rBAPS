function uusiData=addPoints(data)
%Lis‰‰ koordinaatipisteiden joukkoon pisteit‰, jotta jokainen datapiste
%kuuluisi ‰‰relliseen voronoi soluun voronoi tessellaatiota
%muodostettaessa. Apupisteet lis‰t‰‰n muodostamalla hila
%koordinaattipisteiden p‰‰lle ja ottamalla voronoi tessellaatio hilasta. Ne
%hilan pisteet, joita vastaavien solujen sis‰ll‰ ei ole yht‰‰n 
%koordinaattipistett‰, j‰‰v‰t apupisteiksi

x = data(:,1);
y = data(:,2);

xmax = max(x);
xmin = min(x);
ymax = max(y);
ymin = min(y);

npist = size(unique(data, 'rows'),1);
nstep = ceil(npist^0.4) + 7;
xstep = (xmax-xmin)/(nstep-7);
ystep = (ymax-ymin)/(nstep-7);

apuPisteet = zeros(nstep^2,2);

for i=1:nstep
    apuPisteet((i-1)*nstep+1 : i*nstep,1) = xmin + (i-4)*xstep;
    apuPisteet((i-1)*nstep+1 : i*nstep,2) = ymin + ((1:nstep)-4)*ystep;
end



[V,C] = voronoin(apuPisteet,{'Qt','Qbb','Qc','Qz'});

if 0
    figure
    hold on
    for i=1:length(C)
        if isempty(find(C{i} == 1))
            X = V(C{i},:);
            hull = convhull(X(:,1),X(:,2));
            plot(X(hull,1), X(hull,2));
        end
    end
    axis([-2 7 -2 8]);
    plot(data(:,1), data(:,2), 'r*');
    plot(apuPisteet(:,1), apuPisteet(:,2), 'b+');

    hold off
end
empty = zeros(nstep^2,1);

for i = 1:length(C)
    if isempty(find(C{i} == 1))  %Tutkitaan vain rajoitetut solut
        vx = V(C{i},1);
        vy = V(C{i},2);
        IN = any(inpolygon(x,y,vx,vy));
        if IN == 0
            empty(i) = 1;
        end
               
    end
end

empty = find(empty == 1);
C = C(empty);

apuPisteet = apuPisteet(empty, :);

if 0
    figure
    hold on
    for i=1:length(C)
        if isempty(find(C{i} == 1))
            X = V(C{i},:);
            hull = convhull(X(:,1),X(:,2));
            plot(X(hull,1), X(hull,2));
        end
    end
    plot(data(:,1), data(:,2), 'r*');
    plot(apuPisteet(:,1), apuPisteet(:,2), 'b+');
    axis([-2 7 -2 8]);
    hold off
end

uusiData = [data; apuPisteet];
    

