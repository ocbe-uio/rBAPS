function showColors(n_col)

if n_col>36
    error('Maximum number of colors 36');
end

figure('NumberTitle','off','Name','Colors');

set(gca, 'Xlim', [-.5 , n_col+.5], 'YLim', [0,1], ...
    'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', []);

varit = giveColors(n_col);
for k=1:length(varit)
    h0=patch([k-1 k k k-1], [0 0 1 1], varit(k,:));
    set(h0,'EdgeColor','none');   
end