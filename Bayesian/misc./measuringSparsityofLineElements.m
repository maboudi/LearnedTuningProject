function SparsityOfDecoding(BDseqscore, posteriorProbMatrix, onlyLineElements)


plotheight = 25; %% in cm
plotwidth  = 25;

nor = 1; %% number of events per row
noc = 3; %% number of events per column

mainPanelDimension = 2;
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

sparsity.(f{ii}).replayScore  = zeros(nPBEs, 1);
sparsity.(f{ii}).weightedCorr = zeros(nPBEs, 1);

for pbe = 1: nPBEs
    
    % replay score
    pbeLineElements = periodLineElements{pbe};
    
    pbeLineElements(isnan(pbeLineElements)) = [];
    pbeLineElements(~pbeLineElements) = [];
    
    % Removing the surrounding lower than chance level time bins
    
    firstBin = find(pbeLineElements > chanceLevel, 1, 'first');
    lastBin  = find(pbeLineElements > chanceLevel, 1, 'last');
    
    pbeLineElements = pbeLineElements(firstBin : lastBin);
    
    
    % weighted correlation
    temp = sum(posteriorProbMatrix.(f{ii}).data{pbe}, 2);
    
    firstBin = find(temp > chanceLevel2, 1, 'first');
    lastBin  = find(temp > chanceLevel2, 1, 'last');
    
    temp = temp(firstBin : lastBin);
    
    sparsity.(f{ii}).replayScore(pbe)  = ginicoeff(pbeLineElements);
    sparsity.(f{ii}).weightedCorr(pbe) = ginicoeff(temp);

end

medianSparsity.(f{ii}) = nanmean(sparsity.(f{ii}).replayScore);



% scatter plot: sparsity vs sequence score

% correlation

aa = BDseqscore.(f{ii}).data.replayScore.prctilescore(:, 1);
bb = sparsity.(f{ii}).replayScore;

aa(nanIdx) = [];
bb(nanIdx) = [];

nanIdx = find(isnan(bb));

[rho, pval] = corr(aa, bb);


ax1 = axes('position',mainPanelsPos{1,ii},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');

plot(BDseqscore.(f{ii}).data.replayScore.prctilescore(:, 1), sparsity.(f{ii}).replayScore, '.', 'markersize', 5)

xlim([0 100])
ylim([0 1])

set(gca, 'xtick', [0 100], 'ytick', [0 10], 'xticklabel', [], 'yticklabel', [], 'box' , 'off')

text(20, 0.5, sprintf('rho=%.2f (p=%.3f)', rho, pval), 'fontsize', 14)

title({period, sprintf('(n=%d)', length(xData)), '', '', ''}, 'fontsize', 10)



% marginal distribution of sparsity 

leftHistwidth  = 1/plotwidth;
leftHistheight = mainPanelsPos{1,ii}(4);

leftHistPos(1) = mainPanelsPos{1,ii}(1) - leftHistwidth - 0.005;
leftHistPos(2) = mainPanelsPos{1,ii}(2);
leftHistPos(3) = leftHistwidth;
leftHistPos(4) = leftHistheight;

ax2 = axes('position', leftHistPos, 'XGrid', 'off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');


edges = 0:0.1:1;

binCenters = (0:0.1:1) + diff(edges([1 2]))/2;
binCenters(end) = [];
counts = hist(sparsity.(f{ii}).replayScore, edges);
counts(end) = [];

bar(binCenters, counts, 'facecolor', 'b', 'edgecolor', 'none', 'facealpha', 0.5)

camroll(-90)

xlim([0.05 1.05])

set(gca, 'box', 'off', 'ytick', [], 'ycolor', 'none', 'xtick', [0.05 1.05], 'xticklabel', {'0', '1'})
xlabel('sparsity')




% median sparsity vs sequence score

bottomHistwidth  = mainPanelsPos{1,ii}(3);
bottomHistheight = 1/plotwidth;

bottomHistPos(1) = mainPanelsPos{1,ii}(1);
bottomHistPos(2) = mainPanelsPos{1,ii}(2) - bottomHistheight - 0.005;
bottomHistPos(3) = bottomHistwidth;
bottomHistPos(4) = bottomHistheight;


ax3 = axes('position', bottomHistPos, 'XGrid', 'off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');

[sparsityMed , binCenters] = randomFunc(aa, bb);

plot(binCenters, sparsityMed, '.-r', 'linewidth', 3)
line([0 100], [medianSparsity.(f{ii}) medianSparsity.(f{ii})], 'linestyle','--', 'color', 'r', 'linewidth', 3)


xlim([0 100])

set(gca, 'box', 'off', 'ytick', [], 'yColor', 'none', 'xtick', [0 100], 'xticklabel', {'0', '100'})
xlabel(xAxis)


end

dim = [0.15 0.6 0.5 .1];

annotation('textbox',dim,'String',textOnTop, 'FitBoxToText','on', 'Interpreter', 'none');

print(gcf, 'filename', '-dpdf', '-r0')


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











