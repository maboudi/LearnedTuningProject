function plotBargraph2(PRE, PRE_null, POST, POST_null) 


PRE       = PRE(:);
PRE_null  = PRE_null(:);

POST      = POST(:);
POST_null = POST_null(:);



% STATS 

meanPRE       = mean(PRE);
stdPRE        = std(PRE);
meanPRE_null  = mean(PRE_null);
stdPRE_null   = std(PRE_null);

meanPOST      = mean(POST);
stdPOST       = std(POST);
meanPOST_null = mean(POST_null);
stdPOST_null  = std(POST_null);


% ranksum test: each period compared to null

pvalPRE = ranksum(PRE, PRE_null);
pvalPOST = ranksum(POST, POST_null);


% ranksum test: between periods

pvalPRE_POST = ranksum(PRE, POST);


% % plot the errorbars


figure;
set(gcf, 'Units', 'centimeters', 'position', [0 0 17.6 9])

barhandle = bar([meanPRE meanPRE_null; meanPOST meanPOST_null]);

set(barhandle(1), 'FaceColor', 'k', 'linewidth', 2)
set(barhandle(2), 'FaceColor', [.4 .4 .4], 'linewidth', 2)


hold on

h1 = errorbar((1:2)-0.15, [meanPRE meanPOST], [stdPRE stdPOST], '.');
set(h1, 'Color', 'k', 'LineWidth', 2, 'CapSize', 10, 'Marker', 'none')

h2 = errorbar((1:2)+0.15, [meanPRE_null meanPOST_null], [stdPRE_null stdPOST_null], '.');
set(h2, 'Color', 'k', 'LineWidth', 2, 'CapSize', 10, 'Marker', 'none')



% plot the significance markers

%%%%%%% individual periods

% PRE
sigMarker  = calSigMarker(pvalPRE);
lineY      = max([meanPRE+stdPRE meanPRE_null+stdPRE_null])+0.05;
lineX1     = 1-0.15;
lineX2     = 1+0.15;
textX      = 1;

addSigMarkers(lineX1, lineX2, lineY, textX, sigMarker)


% POST
sigMarker  = calSigMarker(pvalPOST);
lineY      = max([meanPOST+stdPOST meanPOST_null+stdPOST_null])+0.05;
lineX1     = 2-0.15;
lineX2     = 2+0.15;
textX      = 2;

addSigMarkers(lineX1, lineX2, lineY, textX, sigMarker)



%%%%%%% between periods

maxValue = max([meanPRE; meanPOST]);

maxValue = maxValue + max([stdPRE; stdPOST]);


% PRE & POST

sigMarker  = calSigMarker(pvalPRE_POST);    
lineY = maxValue + 0.2;

lineX1     = 1-0.15;
lineX2     = 2-0.15;
textX      = mean([lineX1 lineX2]);

addSigMarkers(lineX1, lineX2, lineY, textX, sigMarker)

set(gca, 'XTickLabel', {'PRE', 'RUN', 'POST'})

legend(barhandle, {'Data', 'shuffle'}, 'box', 'off')



end


function addSigMarkers(lineX1, lineX2, lineY, textX, sigMarker)

hold on
line([lineX1 lineX2], [lineY lineY], 'color', 'k', 'linewidth', 2)
hold on
text(textX, lineY+0.1, sigMarker, 'fontsize', 10, 'HorizontalAlignment', 'center')

end


function sigMarker = calSigMarker(pval)

if pval < 0.001
    sigMarker = '***';
elseif pval < 0.01
    sigMarker = '**';
elseif pval < 0.05
    sigMarker = '*';
else
    sigMarker = 'ns'
end

end
