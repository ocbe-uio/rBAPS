function diffInCliqCounts = computeDiffInCliqCounts(cliques, inds)
% Laskee muutoksen CLIQCOUNTS:ssa (tai SEPCOUNTS:ssa, jos syötteen?
% separators) kun yksilöt inds siirretään.
% diffInCliqcounts on ncliq*1 taulu, joka on CLIQCOUNTS:n sarakkeesta josta
% yksilöt inds siirretään ja lisättäv?sarakkeeseen, johon yksilöt
% siirretään.

% taken from spatial model of Jukka Siren's code
% Lu Cheng
% 15.12.2012

ncliq = size(cliques,1);
diffInCliqCounts = zeros(ncliq,1);
ninds = length(inds);
for i = 1:ninds
    ind = inds(i);
    rivit = sum((cliques == ind),2);
    diffInCliqCounts = diffInCliqCounts + rivit;
end