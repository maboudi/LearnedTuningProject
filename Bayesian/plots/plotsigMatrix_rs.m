function plotsigMatrix_rs(sigMatrix, currTitle, okColorbar)


rsthresholds = 0:0.1:0.9;
clthresholds = 0:0.1:0.9; % covered length of the track

cl = [-2 -0.6]; % color limit
cm = redblue; % red blue colormap 

imagesc(clthresholds, rsthresholds, log10(sigMatrix), cl);  % 

hold on 

 
% gcaxrange = xlim;
% gcayrange = ylim;
% 
% line([.45 .45], [gcayrange(1) gcayrange(2)], 'color', 'w', 'linewidth', 1, 'linestyle', ':')
% line([gcaxrange(1) gcaxrange(2)], [.55 .55], 'color', 'w', 'linewidth', 1, 'linestyle', ':')

% xlabel('covered length')
% ylabel('radon integral')

title(currTitle,'interpreter', 'none')

set(gca,  'YDir', 'normal', ...
          'XTick', [0 0.3 0.6 0.9], ...
          'YTick', [0 0.3 0.6 0.9], ...
          'XTickLabels', {'>0', '>0.3', '>0.6', '>0.9'}, ...
          'YTickLabels', {'>0', '>0.3', '>0.6', '>0.9'}, ...
          'fontsize', 10, 'linewidth', 1.5, 'TickDir', 'in', 'TickLength',[0.02, 0.01])

axis square

      
colormap(flipud(cm))


if okColorbar
    
    h= gca;
    originalSize = get(h, 'position');
    cb = colorbar;

    cb_size = get(cb, 'position');
    cb_size(4) = cb_size(4)*0.7;
    cb_size(1) = originalSize(1) + 0.9* originalSize(3);
    cb_size(2) = cb_size(2) + 0.5* cb_size(4);
    cb_size(3) = 1.5*cb_size(3);
    
    set(cb, 'position', cb_size, ...
            'YDir', 'reverse', ...
            'YTick', [log10(0.01) log10(0.05) log10(0.25)], ...
            'YTickLabel', {'P<0.01', '0.05', 'P>0.25'}, ...
            'fontsize', 10, 'linewidth', 1.5, 'TickDir', 'in', 'TickLength',[0.02, 0.01])
        
%     set(h, 'position', originalSize)
    
end



end