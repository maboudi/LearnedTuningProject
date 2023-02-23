
plotwidth = 20;
plotheight = 29;

f=figure;
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
fontsize = 10;



% Data

leftmargin = 1.5;
rightmargin = 11;
bottommargin = 3;
topmargin = 2;
noc = 4;
nor = 12;
spacer = 0.5;
spacec = 0.5;


positions = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, spacer, spacec);


for ii = 1:(noc*nor)
    
    row = ceil(ii/noc);
    column = mod(ii-(row-1)*noc, noc*row);
    if column == 0; column = noc; end
    
    axes('position', positions{row, column},'XGrid', 'off', 'XMinorGrid', 'off', 'FontSize', 8, 'Box', 'off', 'Layer', 'top', 'linewidth', 1.5);
    
    imagesc(POSTbinnedPBEs.data{ii, 2})
    colormap(flipud(gray))
    
    set(gca, 'xtick', [], 'ytick', [])
    set(gca, 'xticklabel', [], 'yticklabel', [])
    
    caxis([0 3])
end




% Poisson

leftmargin = 11;
rightmargin = 1.5;
bottommargin = 3;
topmargin = 2;
noc = 4;
nor = 12;
spacer = 0.5;
spacec = 0.5;


positions = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, spacer, spacec);


for ii = 1:(noc*nor)
    
    row = ceil(ii/noc);
    column = mod(ii-(row-1)*noc, noc*row);
    if column == 0; column = noc; end
    
    axes('position', positions{row, column},'XGrid', 'off', 'XMinorGrid', 'off', 'FontSize', 8, 'Box', 'off', 'Layer', 'top', 'linewidth', 1.5);
    
    imagesc(POSTbinnedPBEs.p{ii, 2})
    colormap(flipud(gray))
    
    set(gca, 'xtick', [], 'ytick', [])
    set(gca, 'xticklabel', [], 'yticklabel', [])
    
    caxis([0 3])
end

