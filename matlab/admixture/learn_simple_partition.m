function part = learn_simple_partition(ordered_points, fii)
% Goes through all the ways to divide the points into two or three groups.
% Chooses the partition which obtains highest logml.

npoints = length(ordered_points);

% One cluster:
val = calculatePopLogml(ordered_points,fii);
bestValue = val;
best_type = 'single';

% Two clusters:
for i=1:npoints-1
    % The right endpoint of the first cluster.
    val_1 = calculatePopLogml(ordered_points(1:i),fii);
    val_2 = calculatePopLogml(ordered_points(i+1:end),fii);
    total = val_1 + val_2;
    if total>bestValue
        bestValue = total;
        best_type = 'double';
        best_i = i;
    end
end

% Three clusters:
for i=1:npoints-2
    for j=i+1:npoints-1
        val_1 = calculatePopLogml(ordered_points(1:i),fii);
        val_2 = calculatePopLogml(ordered_points(i+1:j),fii);
        val_3 = calculatePopLogml(ordered_points(j+1:end),fii);
        total = val_1 + val_2 + val_3;
        if total>bestValue
            bestValue = total;
            best_type = 'triple';
            best_i = i;
            best_j = j;
        end
    end
end

part = zeros(npoints,1);

switch best_type
    case 'single'
        part = ones(npoints,1);    
    case 'double'
        part(1:best_i) = 1;
        part(best_i+1:end) = 2;
    case 'triple'
        part(1:best_i) = 1;
        part(best_i+1:best_j) = 2;
        part(best_j+1:end) = 3;
end


%------------------------------------------


function val = calculatePopLogml(points,fii)
% Calculates fuzzy (log) marginal likelihood for a population of real
% values using estimate "fii" for the dispersion value, and Jeffreys prior
% for the mean parameter.

n = length(points);
fuzzy_ones = sum(points);
fuzzy_zeros = n-fuzzy_ones;

val = gammaln(1) - gammaln(1 + n/fii) ...
    + gammaln(0.5 + fuzzy_ones/fii) + gammaln(0.5 + fuzzy_zeros/fii) ...
    - gammaln(0.5) - gammaln(0.5);