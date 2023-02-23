function plotsigMatrix(sigMatrix, currTitle, okColorbar)



wcthresholds = 0:0.1:0.9;
jdthresholds = 0.1:0.1:1;

cl = [-2 -0.6]; % color limit
cm = redblue; % red blue colormap 

imagesc(jdthresholds, wcthresholds, log10(sigMatrix), cl);  % 

hold on 


gcaxrange = xlim;
gcayrange = ylim;

line([.45 .45], [gcayrange(1) gcayrange(2)], 'color', 'w', 'linewidth', 1, 'linestyle', ':')
line([gcaxrange(1) gcaxrange(2)], [.55 .55], 'color', 'w', 'linewidth', 1, 'linestyle', ':')


% xlabel('jump distance')
% ylabel('weighted correlation')

title(currTitle,'interpreter', 'none')

set(gca,  'YDir', 'reverse', ...
          'XTick', [0.2 0.5 0.8], ...
          'YTick', [0 0.3 0.6 0.9], ...
          'XTickLabels', {'<0.2', '<0.5', '<0.8'}, ...
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
            'YTick', [log10(0.001) log10(0.01) log10(0.05) log10(0.25)], ...
            'YTickLabel', {'P < 0.001', '0.01', '0.05', 'P>0.25'}, ...
            'fontsize', 10, 'linewidth', 1.5, 'TickDir', 'in', 'TickLength',[0.02, 0.01])
        
        %             'YTick', [log10(0.001) log10(0.005) log10(0.05) log10(0.5) log10(1)], ...
%             'YTickLabel', {'P < 0.001', '0.005', '0.05', '0.5', 'P = 1'}, ...
        
%     set(h, 'position', originalSize)
    
end



end