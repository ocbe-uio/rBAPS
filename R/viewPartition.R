viewPartition <- function(osuudet, popnames, COUNTS = matrix(0, 0, 0)) {

    npops <- size(COUNTS, 3)
    nind <- size(osuudet,1)

# TODO: translate if necessary. Remove if this function won't be used 
#     disp(['Number of populations: ' num2str(npops)]);
#     if npops>30
#         disp(' ');
#         disp('Figure can be drawn only if the number of populations');
#         disp('is less or equal to 30.');
#         disp(' ');
#         return;
#     end
    
    
#     varit = givecolors(npops);
#     korkeinviiva = 1.05;
#     pieninarvo = -korkeinviiva;
    
    
#     h0 = figure;
#     set(h0, 'NumberTitle', 'off'); %image_figure;   %Muutettu
#     tiedot.popnames = popnames;
#     tiedot.info = osuudet;
#     set(h0,'UserData',tiedot);
    
#     set(gca, 'Xlim', [-.5 ,nind+.5], 'YLim', [pieninarvo ,korkeinviiva], ...
#         'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', []);
    
#     for i=1:nind
        
#         if any(osuudet(i,:)>0)
#             cumOsuudet = cumsum(osuudet(i,:));
        
#             % Pylv��n piirt�minen
#             for j=1:npops
#                 if j==1
#                     if cumOsuudet(1)>0
#                         h0 =patch([i-1, i, i, i-1], [0, 0, cumOsuudet(1), cumOsuudet(1)], varit(j,:));
#                         set(h0,'EdgeColor','none'); % Midevaa varten kommentoitava!
#                     end
#                 else
#                     if (cumOsuudet(j)>cumOsuudet(j-1))
#                         h0 = patch([i-1, i, i, i-1], [cumOsuudet(j-1), cumOsuudet(j-1), ...
#                         cumOsuudet(j), cumOsuudet(j)], varit(j,:));
#                         set(h0,'EdgeColor','none'); % Midevaa varten kommentoitava!
#                     end
#                 end
#             end
#         end
#     end
    
    
    
#     if ~isempty(popnames)
#         npops = size(popnames,1);
#         for i=1:npops
#             firstInd = popnames{i,2};
#             line([firstInd-1, firstInd-1], [0,1], 'Color', 'k');  %Populaatioiden rajat
                
#             if i<npops
#                 x_paikka = popnames{i,2}-1+(popnames{i+1,2}-popnames{i,2})/2;
#             else
#                 x_paikka = popnames{i,2}-1+(nind+1-popnames{i,2})/2;
#             end
                   
#             korkeuskerroin = pieninarvo / -0.2;
#             suhdekerroin = npops/6;
#             for letter_num = 1:length(popnames{i,1}{1})
#                 letter= popnames{i,1}{1}(letter_num);%alter .004|
#                 text(x_paikka+korjaus(letter)*suhdekerroin, ...
#                     0.0005*korkeuskerroin-0.02*letter_num*korkeuskerroin, ...
#                     letter, 'Interpreter','none');
#             end
#         end
#         line([nind,nind],[0,1],'Color','k');
#     end
}

korjaus <- function(letter) {
    if (any(letter %in% c('i', 'j', 'l', 'I'))) {
            extra <- 0.022
    } else if (any(letter == 'r')) {
            extra <- 0.016
    } else if (any(letter == 'k')) {
            extra <- 0.009
    } else if (any(letter == 'f')) {
            extra <- 0.013
    } else if (any(letter == 't')) {
            extra <- 0.014
    } else if (any(letter == 'w')) {
            extra <- -0.003
    } else {
            extra <- 0
    }
    return(extra)
}

giveColors <- function(n) {
        if (n > 36) stop('Maximum number of colors 36')
        colors <- matrix(
                data = c(
                        1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1,
                        0.4, 0, 0, 0, 0.4, 0, 0, 0, 0.4, 0.4, 0.4, 0, 0.4, 0,
                        0.4, 0, 0.4, 0.4, 0.2, 0, 0, 0, 0.2, 0, 0, 0, 0.2, 0.2,
                        0.2, 0, 0.2, 0, 0.2, 0, 0.2, 0.2, 0.8, 0, 0, 0, 0.8, 0,
                        0, 0, 0.8, 0.8, 0.8, 0, 0.8, 0, 0.8, 0, 0.8, 0.8, 
                        0.6, 0, 0, 0, 0.6, 0, 0, 0, 0.6, 0.6, 0.6, 0, 0.6, 0,
                        0.6, 0, 0.6, 0.6, 0.6, 0.2, 0.4, 0.2, 0.4, 0.8, 0.8,
                        0.4, 0.2, 0, 0.6, 0.2, 0.2, 0.8, 0.6, 0.5, 0.2, 0.1,
                        0.6, 0.3, 0.1
                ),
                ncol = 3,
                byrow = TRUE
        )
        colors = colors[1:n, ]
        # red; green; blue; yellow
        # RGB format: [red green blue]
        return(colors)
}