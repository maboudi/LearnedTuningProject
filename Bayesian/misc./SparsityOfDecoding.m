function giniOfDecoding(BDseqscore, posteriorProbMatrix, onlyLineElements, begPosition, endPosition, method, textOnTop)


plotheight = 25; %% in cm
plotwidth  = 25;

nor = 1; %% number of events per row
noc = 3; %% number of events per column

mainPanelDimension = 3;
gapr = 3;     gapc = 4;

leftSpaceH = plotwidth - (noc*mainPanelDimension + (noc-1)*gapc);
leftSpaceV = plotheight - (nor*mainPanelDimension + (nor-1)*gapr);

leftedge = leftSpaceH/2;    rightedge = leftSpaceH/2;    topedge = leftSpaceV/2;    bottomedge = leftSpaceV/2;


fontsize = 5;
mainPanelsPos  = subplot_pos(plotwidth, plotheight, leftedge, rightedge, bottomedge, topedge, noc, nor, gapr, gapc); 


f = figure;
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);


f = fieldnames(onlyLineElements);

for ii = 1: length(f)
    
    period = f{ii};
    

% calculating sparisty based on the best fit line (in case of replay score)
% or the whole posterior probaility matrix in case of weighted correlation

periodLineElements = onlyLineElements.(f{ii}).data(:, 1);


nPBEs    = size(periodLineElements, 1);
nPosBins = size(periodLineElements{1}, 1);
bestFitBandWidth = 15; 
chanceLevel  = bestFitBandWidth/nPosBins;
chanceLevel2 = 1/nPosBins;

gini.(f{ii}).replayScore  = zeros(nPBEs, 1);
gini.(f{ii}).weightedCorr = zeros(nPBEs, 1);

%%% second method for replay score

sm = ones(15,1);

%%%

for pbe = 1: nPBEs
    
    % replay score
    pbeLineElements = periodLineElements{pbe};
    
    
    
    %%%%% second method
    temp = zeros(size(pbeLineElements));
    for jj = 1:size(pbeLineElements, 2)
        temp(:, jj) = conv(pbeLineElements(:, jj), sm, 'same');
    end
    
    temp(temp > 0) = 1;
    
    temp2 = temp .* posteriorProbMatrix.(f{ii}).data{pbe};
    temp2(isnan(temp)) = [];
    
    temp3 = sum(temp2, 2);
    
    firstBin = find(temp3 > chanceLevel2, 1, 'first');
    lastBin  = find(temp3 > chanceLevel2, 1, 'last');
    
    temp3 = temp3(firstBin: lastBin);
    
%     gini.(f{ii}).replayScore(pbe)  = ginicoeff(temp3);
    
    %%%%
    pbeLineElements = sum(pbeLineElements, 1);
    pbeLineElements(isnan(pbeLineElements)) = [];
    
%     pbeLineElements(~pbeLineElements) = [];
    
    
    
    % Removing the surrounding lower than chance level time bins
    
    firstBin = find(pbeLineElements > chanceLevel, 1, 'first');
    lastBin  = find(pbeLineElements > chanceLevel, 1, 'last');
    
    pbeLineElements = pbeLineElements(firstBin : lastBin);
    
    
    % weighted correlation
    temp = sum(posteriorProbMatrix.(f{ii}).data{pbe}, 2);
    
    firstBin = find(temp > chanceLevel2, 1, 'first');
    lastBin  = find(temp > chanceLevel2, 1, 'last');
    
    temp = temp(firstBin : lastBin);
    
    gini.(f{ii}).replayScore(pbe)  = ginicoeff(pbeLineElements);
    
%     if gini.(f{ii}).replayScore(pbe) > 0.6
%         gini.(f{ii}).replayScore(pbe)
%     end
    
    gini.(f{ii}).weightedCorr(pbe) = ginicoeff(temp);

end

mediangini.(f{ii}) = nanmean(gini.(f{ii}).(method));



% scatter plot: gini vs sequence score

% correlation

aa = BDseqscore.(f{ii}).data.(method).prctilescore(:, 1);
bb = gini.(f{ii}).(method);

nanIdx = find(isnan(bb));

aa(nanIdx) = [];
bb(nanIdx) = [];

nanIdx = find(isnan(bb));

[rho, pval] = corr(aa, bb);


ax1 = axes('position',mainPanelsPos{1,ii},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');

plot(BDseqscore.(f{ii}).data.(method).prctilescore(:, 1), gini.(f{ii}).(method), '.', 'markersize', 3)

xlim([0 100])
ylim([0 1])

set(gca, 'xtick', [0 100], 'ytick', [0 10], 'xticklabel', [], 'yticklabel', [], 'box' , 'off')

if pval < 0.001
text(10, 0.1, sprintf('rho=%.2f\np<1e%d', rho, ceil(log10(pval))), 'fontsize', 7)
else
   text(10, 0.1, sprintf('rho=%.2f\np=%.3f', rho, pval) , 'fontsize', 7)
end

title({period, sprintf('(n=%d)', length(posteriorProbMatrix.(f{ii}).data)), ''}, 'fontsize', 10)



% marginal distribution of gini 

leftHistwidth  = 1/plotwidth;
leftHistheight = mainPanelsPos{1,ii}(4);

leftHistPos(1) = mainPanelsPos{1,ii}(1) - leftHistwidth - 0.005;
leftHistPos(2) = mainPanelsPos{1,ii}(2);
leftHistPos(3) = leftHistwidth;
leftHistPos(4) = leftHistheight;

ax2 = axes('position', leftHistPos, 'XGrid', 'off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');


edges = 0:0.05:1;

binCenters = edges(1:end-1) + 0.025;

counts = hist(gini.(f{ii}).(method), edges);
counts(end - 1) = counts(end - 1) + counts(end);
counts(end) = [];

counts = fliplr(counts);

bar(binCenters, counts, 'facecolor', 'b', 'edgecolor', 'none', 'facealpha', 0.5)
line([1-mediangini.(f{ii}) 1-mediangini.(f{ii})], [0 max(counts)], 'color', 'r', 'linewidth', 2)

camroll(-90)

xlim([0 1])

set(gca, 'box', 'off', 'ytick', [], 'ycolor', 'none', 'xtick', [0 1], 'xticklabel', {sprintf('%.1f', 1), '0'})
xlabel('gini')




% median gini vs sequence score

bottomHistwidth  = mainPanelsPos{1,ii}(3);
bottomHistheight = 1/plotwidth;

bottomHistPos(1) = mainPanelsPos{1,ii}(1);
bottomHistPos(2) = mainPanelsPos{1,ii}(2) - bottomHistheight - 0.005;
bottomHistPos(3) = bottomHistwidth;
bottomHistPos(4) = bottomHistheight;


ax3 = axes('position', bottomHistPos, 'XGrid', 'off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');

[giniMed , binCenters] = randomFunc(aa, bb);

plot(binCenters, giniMed, '.-r', 'linewidth', 2)
line([0 100], [mediangini.(f{ii}) mediangini.(f{ii})], 'linestyle','--', 'color', 'r', 'linewidth', 1)


xlim([0 100])
ylim([min(giniMed)-0.03 max(giniMed)+0.03])
set(gca, 'box', 'off', 'ytick', [min(giniMed) max(giniMed)], 'yticklabel', {sprintf('%.2f', min(giniMed)), sprintf('%.2f', max(giniMed))}, 'xtick', [0 100], 'xticklabel', {'0', '100'})
xlabel('sequence score')
ylabel('median gini')

end
%%%


nor = 1; %% number of events per row
noc = 4; %% number of events per column

mainPanelDimension = 3;
gapr = 3;     gapc = 2;

leftSpaceH = plotwidth - (noc*mainPanelDimension + (noc-1)*gapc);
leftSpaceV = plotheight - (nor*mainPanelDimension + (nor-1)*gapr);

leftedge = leftSpaceH/2;    rightedge = leftSpaceH/2;    topedge = 18;    bottomedge = 5;

mainPanelsPos2  = subplot_pos(plotwidth, plotheight, leftedge, rightedge, bottomedge, topedge, noc, nor, gapr, gapc); 

begPosition1 = begPosition.RUN.data(41,:);
endPosition1 = endPosition.RUN.data(41,:);

ax1 = axes('position',mainPanelsPos2{1,1},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');

imagesc(posteriorProbMatrix.RUN.data{41}, [0 0.167]); set(gca, 'yDir', 'normal')
set(gca, 'xtick', [], 'ytick', [])
xlabel('time bin')
ylabel('position bin')

title('RUN, pbe 41, repscorep=99.5', 'fontsize', 8)
colormap('jet')

line([begPosition1(1) endPosition1(1)], [begPosition1(2) endPosition1(2)], 'linewidth', 1, 'color', 'w')
line([begPosition1(1) endPosition1(1)], [begPosition1(2)-7 endPosition1(2)-7], 'linestyle', ':', 'linewidth', 1, 'color', 'w')
line([begPosition1(1) endPosition1(1)], [begPosition1(2)+7 endPosition1(2)+7], 'linestyle', ':', 'linewidth', 1, 'color', 'w')


ax2 = axes('position',mainPanelsPos2{1,2},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');

temp = zeros(size(posteriorProbMatrix.RUN.data{41}));
for ii = 1: size(posteriorProbMatrix.RUN.data{41}, 2)
    temp(:, ii) = conv(onlyLineElements.RUN.data{41, 1}(:, ii), ones(15,1), 'same');
end
    
temp = posteriorProbMatrix.RUN.data{41}.*temp;

imagesc(temp, [0 0.167]); set(gca, 'yDir', 'normal')
set(gca, 'xtick', [], 'ytick', [])
xlabel('time bin')
ylabel('position bin')
colormap('jet')


line([begPosition1(1) endPosition1(1)], [begPosition1(2) endPosition1(2)], 'linewidth', 1, 'color', 'w')
line([begPosition1(1) endPosition1(1)], [begPosition1(2)-7 endPosition1(2)-7], 'linestyle', ':', 'linewidth', 1, 'color', 'w')
line([begPosition1(1) endPosition1(1)], [begPosition1(2)+7 endPosition1(2)+7], 'linestyle', ':', 'linewidth', 1, 'color', 'w')


ax3 = axes('position',mainPanelsPos2{1,3},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');


pbeLineElements = onlyLineElements.RUN.data{41, 1};

pbeLineElements = sum(pbeLineElements, 1);
% pbeLineElements(isnan(pbeLineElements)) = [];

firstBin = find(pbeLineElements > chanceLevel, 1, 'first');
lastBin  = find(pbeLineElements > chanceLevel, 1, 'last');

imagesc(onlyLineElements.RUN.data{41, 1}, [0 0.167*5])
set(gca, 'yDir', 'normal')
set(gca, 'xtick', [], 'ytick', [])
xlabel('time bin')
ylabel('position bin')
colormap('jet')

text(firstBin, find(onlyLineElements.RUN.data{41, 1}(:, firstBin)) + 10, '*', 'fontsize', 10, 'color', 'w')
text(lastBin, find(onlyLineElements.RUN.data{41, 1}(:, lastBin)) + 10, '*', 'fontsize', 10, 'color', 'w')


line([begPosition1(1) endPosition1(1)], [begPosition1(2) endPosition1(2)], 'linestyle', ':', 'linewidth', 1, 'color', 'w')


ax4 = axes('position',mainPanelsPos2{1,4},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');

pbeLineElements = pbeLineElements(firstBin : lastBin);

plot(pbeLineElements, 'b', 'linewidth', 2)
text(2, min(pbeLineElements)+range(pbeLineElements)/2, sprintf('gini coeff=%.2f', ginicoeff(pbeLineElements)))
xlim([1 length(pbeLineElements)])
xticks([1 length(pbeLineElements)])
xlabel('time bin')
ylabel('summed posterior')


annotation('arrow', [sum(mainPanelsPos2{1,1}([1 3]))+0.01 mainPanelsPos2{1,2}(1)- 0.03], [mainPanelsPos2{1,1}(2)+mainPanelsPos2{1,1}(4)/2 mainPanelsPos2{1,1}(2)+mainPanelsPos2{1,1}(4)/2], 'linewidth', 2)
annotation('arrow', [sum(mainPanelsPos2{1,2}([1 3]))+0.005 mainPanelsPos2{1,3}(1)- 0.03], [mainPanelsPos2{1,1}(2)+mainPanelsPos2{1,1}(4)/2 mainPanelsPos2{1,1}(2)+mainPanelsPos2{1,1}(4)/2], 'linewidth', 2)
annotation('arrow', [sum(mainPanelsPos2{1,3}([1 3]))+0.005 mainPanelsPos2{1,4}(1)- 0.03], [mainPanelsPos2{1,1}(2)+mainPanelsPos2{1,1}(4)/2 mainPanelsPos2{1,1}(2)+mainPanelsPos2{1,1}(4)/2], 'linewidth', 2)



nor = 1; %% number of events per row
noc = 4; %% number of events per column

mainPanelDimension = 3;
gapr = 3;     gapc = 2;

leftSpaceH = plotwidth - (noc*mainPanelDimension + (noc-1)*gapc);
leftSpaceV = plotheight - (nor*mainPanelDimension + (nor-1)*gapr);

leftedge = leftSpaceH/2;    rightedge = leftSpaceH/2;    topedge = 21;    bottomedge = 2;

mainPanelsPos2  = subplot_pos(plotwidth, plotheight, leftedge, rightedge, bottomedge, topedge, noc, nor, gapr, gapc); 

begPosition2 = begPosition.PRE.data(1528,:);
endPosition2 = endPosition.PRE.data(1528,:);

ax1 = axes('position',mainPanelsPos2{1,1},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');

imagesc(posteriorProbMatrix.PRE.data{1528}, [0 0.167]); set(gca, 'yDir', 'normal')
set(gca, 'xtick', [], 'ytick', [])
xlabel('time bin')
ylabel('position bin')
title('PRE, pbe 1528, repscorep=100', 'fontsize', 8)
colormap('jet')

line([begPosition2(1) endPosition2(1)], [begPosition2(2) endPosition2(2)], 'linewidth', 1, 'color', 'w')
line([begPosition2(1) endPosition2(1)], [begPosition2(2)-7 endPosition2(2)-7], 'linestyle', ':', 'linewidth', 1, 'color', 'w')
line([begPosition2(1) endPosition2(1)], [begPosition2(2)+7 endPosition2(2)+7], 'linestyle', ':', 'linewidth', 1, 'color', 'w')


ax2 = axes('position',mainPanelsPos2{1,2},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');

temp = zeros(size(posteriorProbMatrix.PRE.data{1528}));
for ii = 1: size(posteriorProbMatrix.PRE.data{1528}, 2)
    temp(:, ii) = conv(onlyLineElements.PRE.data{1528, 1}(:, ii), ones(15,1), 'same');
end
    
temp = posteriorProbMatrix.PRE.data{1528}.*temp;

imagesc(temp, [0 0.167]); set(gca, 'yDir', 'normal')
set(gca, 'xtick', [], 'ytick', [])
xlabel('time bin')
ylabel('position bin')
colormap('jet')


line([begPosition2(1) endPosition2(1)], [begPosition2(2) endPosition2(2)], 'linewidth', 1, 'color', 'w')
line([begPosition2(1) endPosition2(1)], [begPosition2(2)-7 endPosition2(2)-7], 'linestyle', ':', 'linewidth', 1, 'color', 'w')
line([begPosition2(1) endPosition2(1)], [begPosition2(2)+7 endPosition2(2)+7], 'linestyle', ':', 'linewidth', 1, 'color', 'w')


ax3 = axes('position',mainPanelsPos2{1,3},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');


pbeLineElements = onlyLineElements.PRE.data{1528, 1};

pbeLineElements = sum(pbeLineElements, 1);
% pbeLineElements(isnan(pbeLineElements)) = [];

firstBin = find(pbeLineElements > chanceLevel, 1, 'first');
lastBin  = find(pbeLineElements > chanceLevel, 1, 'last');

imagesc(onlyLineElements.PRE.data{1528, 1}, [0 0.167*5])
set(gca, 'yDir', 'normal')
set(gca, 'xtick', [], 'ytick', [])
xlabel('time bin')
ylabel('position bin')
colormap('jet')

text(firstBin, find(onlyLineElements.PRE.data{1528, 1}(:, firstBin)) + 10, '*', 'fontsize', 10, 'color', 'w')
text(lastBin, find(onlyLineElements.PRE.data{1528, 1}(:, lastBin)) + 10, '*', 'fontsize', 10, 'color', 'w')


line([begPosition2(1) endPosition2(1)], [begPosition2(2) endPosition2(2)], 'linestyle', ':', 'linewidth', 1, 'color', 'w')


ax4 = axes('position',mainPanelsPos2{1,4},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');

pbeLineElements = pbeLineElements(firstBin : lastBin);

plot(pbeLineElements, 'b', 'linewidth', 2)
text(2, min(pbeLineElements)+range(pbeLineElements)/2, sprintf('gini coeff=%.2f', ginicoeff(pbeLineElements)))
xlim([1 length(pbeLineElements)])
xticks([1 length(pbeLineElements)])
xlabel('time bin')
ylabel('summed posterior')


annotation('arrow', [sum(mainPanelsPos2{1,1}([1 3]))+0.01 mainPanelsPos2{1,2}(1)- 0.03], [mainPanelsPos2{1,1}(2)+mainPanelsPos2{1,1}(4)/2 mainPanelsPos2{1,1}(2)+mainPanelsPos2{1,1}(4)/2], 'linewidth', 2)
annotation('arrow', [sum(mainPanelsPos2{1,2}([1 3]))+0.005 mainPanelsPos2{1,3}(1)- 0.03], [mainPanelsPos2{1,1}(2)+mainPanelsPos2{1,1}(4)/2 mainPanelsPos2{1,1}(2)+mainPanelsPos2{1,1}(4)/2], 'linewidth', 2)
annotation('arrow', [sum(mainPanelsPos2{1,3}([1 3]))+0.005 mainPanelsPos2{1,4}(1)- 0.03], [mainPanelsPos2{1,1}(2)+mainPanelsPos2{1,1}(4)/2 mainPanelsPos2{1,1}(2)+mainPanelsPos2{1,1}(4)/2], 'linewidth', 2)






%%%
dim = [0.1 0.65 0.5 .1];

annotation('textbox',dim,'String',textOnTop, 'FitBoxToText','on', 'Interpreter', 'none');

print(gcf, 'filename3', '-dpdf', '-r0')



end


function  [var2Med , binCenters] = randomFunc(var1, var2)

edges = 0:10:100;
binCenters = edges(1:end-1) + 5;

[~, binIdx] = histc(var1, edges);
binIdx(binIdx == 11) = 10;

var2Med = zeros(10, 1);

for ii = 1: 10
    currAmounts = var2(binIdx == ii);
    
    var2Med(ii) = median(currAmounts);
    
end

end











