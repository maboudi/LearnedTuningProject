function [pathLens, pathLens_c, pathLens_ts, pathLens_p, pathLens_i] = pathLength(transmat, transmat_c, transmat_ts, transmat_p, transmat_i, ii)

% path analysis

transProbThresh = 0.2;
overlapThresh = 0.5;

[~, pathLens] = pathAnalysis(transmat(:,:,ii), transProbThresh, overlapThresh, [], [] , []); %% positionsRL positionsLR positionBins

[~, pathLens_c] = pathAnalysis(transmat_c(:,:,ii), transProbThresh, overlapThresh, [], [] , []);

[~, pathLens_ts] = pathAnalysis(transmat_ts(:,:,ii), transProbThresh, overlapThresh, [], [] , []);

[~, pathLens_p] = pathAnalysis(transmat_p(:,:,ii), transProbThresh, overlapThresh, [], [] , []);

[~, pathLens_i] = pathAnalysis(transmat_i(:,:,ii), transProbThresh, overlapThresh, [], [] , []);



% histograms

temp = [pathLens; pathLens_c; pathLens_ts; pathLens_p; pathLens_i];
bins = linspace(min(temp(:)), max(temp(:)), 15);

realH = hist(pathLens, bins);
UnitCycleH = hist(pathLens_c, bins);
timeSwapH = hist(pathLens_ts, bins);
poissH = hist(pathLens_p, bins);
identityH = hist(pathLens_i, bins);

maxCount = max([realH UnitCycleH timeSwapH poissH identityH]);
% statistical tests for comparing the distributions

% realVScycleP = ranksum(pathLens, pathLens_c);
% realVStimeSwapP = ranksum(pathLens, pathLens_ts);
% realVSpoissP = ranksum(pathLens, pathLens_p);
% realVSidentityP = ranksum(pathLens, pathLens_i);


% % smoothing and plot the distributions
% sigma = 2;
% halfwidth = 2 * sigma;
% smoothwin = gausswindow(sigma, halfwidth);
% 
% realH = conv(realH, smoothwin, 'same');
% UnitCycleH = conv(UnitCycleH, smoothwin, 'same');
% timeSwapH = conv(timeSwapH, smoothwin, 'same');
% poissH = conv(poissH, smoothwin, 'same');
% identityH = conv(identityH, smoothwin, 'same');

% realH = realH ./ sum(realH);
% UnitCycleH = UnitCycleH ./ sum(UnitCycleH);
% timeSwapH = timeSwapH ./ sum(timeSwapH);
% poissH = poissH ./ sum(poissH);
% identityH = identityH ./ sum(identityH);

% figure
%%%
hold on

bar_handle = bar(bins, [realH' UnitCycleH' timeSwapH' poissH' identityH']);

set(bar_handle(1),'faceColor','k', 'Edgecolor', 'none', 'BarWidth', 1)
set(bar_handle(2),'faceColor','c', 'Edgecolor', 'none', 'BarWidth', 1)
set(bar_handle(3),'faceColor','b', 'Edgecolor', 'none', 'BarWidth', 1)
set(bar_handle(4),'faceColor','g', 'Edgecolor', 'none', 'BarWidth', 1)
set(bar_handle(5),'faceColor','m', 'Edgecolor', 'none', 'BarWidth', 1)



% 
% plot(bins, realH, 'k', 'linewidth', 0.5)
% plot(bins, UnitCycleH, 'c','linewidth', 0.5)
% plot(bins, timeSwapH, 'b','linewidth', 0.5)
% plot(bins, poissH, 'g','linewidth', 0.5)
% plot(bins, identityH, 'm','linewidth', 0.5)

% 
% % significance bars
% [realPy, realPx]= max(realH);
% [UnitCyclePy, UnitCyclePx]= max(UnitCycleH);
% [timeSwapPy, timeSwapPx]= max(timeSwapH);
% [poissPy, poissPx]= max(poissH);
% [identityPy, identityPx]= max(identityH);
% 
% firstmax = max([realPy, UnitCyclePy, timeSwapPy, poissPy, identityPy])*1.02;
% textsize = 16;
% % Wilcoxon ranksum test
% if realVScycleP < 0.05
%     y = firstmax;
%     plot([bins(realPx) bins(UnitCyclePx)],[y y],'k','linewidth', 3)
%     
%     if realVScycleP < 0.01
%         psign = '**';
%     else
%         psign = '*';
%     end
%         
%     text(mean([bins(realPx), bins(UnitCyclePx)]), y+0.15, psign, 'fontsize', textsize)
% 
% end
% 
% if realVStimeSwapP < 0.05
%     y = firstmax+0.5;
%     plot([bins(realPx) bins(timeSwapPx)],[y y],'k','linewidth', 3)
%     
%     if realVStimeSwapP < 0.01
%         psign = '**';
%     else
%         psign = '*';
%     end
%         
%     text(mean([bins(realPx), bins(timeSwapPx)]), y+0.15, psign, 'fontsize', textsize)
% 
% end
% 
% if realVSpoissP < 0.05
%     y = firstmax+1;
%     plot([bins(realPx) bins(poissPx)],[y y],'k','linewidth', 3)
%     
%     if realVSpoissP < 0.01
%         psign = '**';
%     else
%         psign = '*';
%     end
%         
%     text(mean([bins(realPx), bins(poissPx)]), y+0.15, psign, 'fontsize', textsize)
% 
%     
% end
% 
% if realVSidentityP < 0.05
%     y = firstmax+1.5;
%     plot([bins(realPx) bins(identityPx)],[y y],'k','linewidth', 3)
%     
%     if realVSidentityP < 0.01
%         psign = '**';
%     else
%         psign = '*';
%     end
%         
%     text(mean([bins(realPx), bins(identityPx)]), y+0.15, psign, 'fontsize', textsize)
% 
% end

%%
facx = max(temp(:))/25;
facy = maxCount/8;

line([max(temp(:))-facx max(temp(:))+facx],[maxCount+3*facy maxCount+3*facy], 'color', 'k', 'linewidth', 1)
line([max(temp(:))-facx max(temp(:))+facx],[maxCount+2*facy maxCount+2*facy], 'color', 'c', 'linewidth', 1)
line([max(temp(:))-facx max(temp(:))+facx],[maxCount+1*facy maxCount+1*facy], 'color', 'b', 'linewidth', 1)
line([max(temp(:))-facx max(temp(:))+facx],[maxCount maxCount], 'color', 'g', 'linewidth', 1)
line([max(temp(:))-facx max(temp(:))+facx],[maxCount-1*facy maxCount-1*facy], 'color', 'm', 'linewidth', 1)

text(max(temp(:))+1.2*facx, maxCount+3*facy, 'Real', 'fontsize', 8)
text(max(temp(:))+1.2*facx, maxCount+2*facy, 'Incoherent', 'fontsize', 8)
text(max(temp(:))+1.2*facx, maxCount+1*facy, 'Coherent', 'fontsize', 8)
text(max(temp(:))+1.2*facx, maxCount, 'Poisson', 'fontsize', 8)
text(max(temp(:))+1.2*facx, maxCount-1*facy, 'Unit identity', 'fontsize', 8)

ylim([0 maxCount+3*facy])

hold off

set(gca, 'fontsize', 10)
xlabel('Path length', 'fontsize', 10)
ylabel('Counts', 'fontsize', 10)
title('Connectivity', 'fontsize', 10)

