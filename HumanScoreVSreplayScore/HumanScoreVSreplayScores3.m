
% clear
% close all

load('rndPBEs_new.mat')
load('rndPBEs.mat')


for ii = 1:length(rndPBEs_new)
    rndPBEs_new(ii).PBEno = rndPBEs(ii).PBEno;    
end

rndPBEs = rndPBEs_new;


subjects = {'KM'; 'BG2'; 'NK'; 'JA'; 'UK'; 'PH';'EP'; 'KD'};

% subjects = {'JA'};

for ss = 1: numel(subjects)
    
   load(sprintf('HumanScores_%s', subjects{ss}), 'scores')
   

   [scores2, uniqPBEidx] = scoreConsistency(scores, [rndPBEs.PBEno]);
   scores2(:,3) = abs(scores2(:,1) - scores2(:, 2));
   
   
   
   
   % null distributions of consistency 
   
   all_surr_diff = [];
   
   for ii = 1:1000
       
      surr_scores  = scores(randperm(numel(scores)));
      surr_scores2 = scoreConsistency(surr_scores, [rndPBEs.PBEno]);
      
      surr_diff = abs(surr_scores2(:,1) - surr_scores2(:, 2));
      surr_diff = surr_diff(~isnan(surr_diff));
      
      all_surr_diff = [all_surr_diff; surr_diff];
      
   end
   
   
   prctile = nan(size(scores2, 1), 1);
   
   for pbe = 1:size(scores2, 1)
       
       if isnan(scores2(pbe, 3))
          continue
       end
       
       prctile(pbe) = numel(find(all_surr_diff >= scores2(pbe, 3)))/numel(all_surr_diff)*100;
       
   end
    
   scores2(:, 4) = prctile;
   scores2(:, 5) = ceil(nanmedian(scores2(:, 1:2), 2));
   
   scores_all.(subjects{ss}) = scores2;
    
end


%%

rndPBEs = rndPBEs(uniqPBEidx);

rt_ts = [rndPBEs.rt_ts]';
rt_ui = [rndPBEs.rt_ui]';
rt_pf = [rndPBEs.rt_pf]';
rt_ds = [rndPBEs.rt_ds]';

wc_ts = [rndPBEs.wc_ts]';
wc_ui = [rndPBEs.wc_ui]';
wc_pf = [rndPBEs.wc_pf]';
wc_ds = [rndPBEs.wc_ds]';

PRE  = [rndPBEs.PRE]';
RUN  = [rndPBEs.RUN]';
POST = [rndPBEs.POST]';





%% plotting the figure for each subject separately


for ii = 8%:numel(subjects)
    
currSubject = subjects{ii};
scores = scores_all.(currSubject)(:, 5);

plotwidth = 20;
plotheight = 29;

f=figure;
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
fontsize = 10;



% overall distribution of scores

leftmargin = 1.5;
rightmargin = 13;
bottommargin = 22;
topmargin = 1.5;
noc = 1;
nor = 1;
spacer = 0;
spacec = 0;

positions = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, spacer, spacec);
axes('position',positions{1},'XGrid','off','XMinorGrid','off','FontSize',8,'Box','off','Layer','top','linewidth', 1.5);

binEdges = 0:9;
h_PRE = histc(scores(PRE==1), binEdges);
h_RUN = histc(scores(RUN==1), binEdges);
h_POST = histc(scores(POST==1), binEdges);

bar(binEdges, [h_PRE h_RUN h_POST], 'stacked')

legend('PRE', 'RUN', 'POST', 'location', 'best')
xlabel('human replay score')
ylabel('number of PBEs')
set(gca, 'box', 'off', 'linewidth', 1.5)


% % scatterplot of algorithmic replay scores v human replays score

rt_ts_med = zeros(10,1);
rt_ui_med = zeros(10,1);
wc_ts_med = zeros(10,1);
wc_ui_med = zeros(10,1);

rt_pf_med = zeros(10,1);
rt_ds_med = zeros(10,1);
wc_pf_med = zeros(10,1);
wc_ds_med = zeros(10,1);

for ihscore = 1:10
    
    currhscore = ihscore-1;
    rt_ts_med(ihscore) = nanmedian(rt_ts(scores == currhscore));
    rt_ui_med(ihscore) = nanmedian(rt_ui(scores == currhscore));
    wc_ts_med(ihscore) = nanmedian(wc_ts(scores == currhscore));
    wc_ui_med(ihscore) = nanmedian(wc_ui(scores == currhscore));
    
    rt_pf_med(ihscore) = nanmedian(rt_pf(scores == currhscore));
    rt_ds_med(ihscore) = nanmedian(rt_ds(scores == currhscore));
    wc_pf_med(ihscore) = nanmedian(wc_pf(scores == currhscore));
    wc_ds_med(ihscore) = nanmedian(wc_ds(scores == currhscore));
end
    


leftmargin = 9;
rightmargin = 1.5;
bottommargin = 21;
topmargin = 1.5;
noc = 4;
nor = 2;
spacer = 1;
spacec = 0.5;

positions = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, spacer, spacec);


axes('position',positions{1,1},'XGrid','off','XMinorGrid','off','FontSize',8,'Box','off','Layer','top', 'linewidth', 1.5);

hold on
scatter(scores+randn(size(scores))*0.1, rt_ts, 5, 'filled')
scatter(0:9, rt_ts_med, 15, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'Marker', 'd')

xlim([-1 10])
ylim([-5 105])

title('RT-TS', 'fontsize', fontsize)


axes('position',positions{1,2},'XGrid','off','XMinorGrid','off','FontSize',8,'Box','off','Layer','top', 'linewidth', 1.5);

hold on
scatter(scores+randn(size(scores))*0.1, rt_ui, 5, 'filled')
scatter(0:9, rt_ui_med, 15, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'Marker', 'd')

xlim([-1 10])
ylim([-5 105])

title('RT-UI', 'fontsize', fontsize)



axes('position',positions{1,3},'XGrid','off','XMinorGrid','off','FontSize',8,'Box','off','Layer','top', 'linewidth', 1.5);

hold on
scatter(scores+randn(size(scores))*0.1, wc_ts, 5, 'filled')
scatter(0:9, wc_ts_med, 15, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'Marker', 'd')
xlim([-1 10])
ylim([-5 105])

title('WC-TS', 'fontsize', fontsize)



axes('position',positions{1,4},'XGrid','off','XMinorGrid','off','FontSize',8,'Box','off','Layer','top', 'linewidth', 1.5);

hold on
scatter(scores+randn(size(scores))*0.1, wc_ui, 5, 'filled')
scatter(0:9, wc_ui_med, 15, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'Marker', 'd')
xlim([-1 10])
ylim([-5 105])

title('WC-UI', 'fontsize', fontsize)



axes('position',positions{2,1},'XGrid','off','XMinorGrid','off','FontSize',8,'Box','off','Layer','top', 'linewidth', 1.5);

hold on
scatter(scores+randn(size(scores))*0.1, rt_pf, 5, 'filled')
scatter(0:9, rt_pf_med, 15, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'Marker', 'd')
xlim([-1 10])
ylim([-5 105])

ylabel('replay score %', 'fontsize', fontsize)
xlabel('human replay score', 'fontsize', fontsize)
title('RT-PF', 'fontsize', fontsize)



axes('position',positions{2,2},'XGrid','off','XMinorGrid','off','FontSize',8,'Box','off','Layer','top', 'linewidth', 1.5);

hold on
scatter(scores+randn(size(scores))*0.1, rt_ds, 5, 'filled')
scatter(0:9, rt_ds_med, 15, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'Marker', 'd')
xlim([-1 10])
ylim([-5 105])

title('RT-DS', 'fontsize', fontsize)


axes('position',positions{2,3},'XGrid','off','XMinorGrid','off','FontSize',8,'Box','off','Layer','top', 'linewidth', 1.5);

hold on
scatter(scores+randn(size(scores))*0.1, wc_pf, 5, 'filled')
scatter(0:9, wc_pf_med, 15, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'Marker', 'd')
xlim([-1 10])
ylim([-5 105])

title('WC-PF', 'fontsize', fontsize)



axes('position',positions{2,4},'XGrid','off','XMinorGrid','off','FontSize',8,'Box','off','Layer','top', 'linewidth', 1.5);

hold on
scatter(scores+randn(size(scores))*0.1, wc_ds, 5, 'filled')
scatter(0:9, wc_ds_med, 15, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'Marker', 'd')
xlim([-1 10])
ylim([-5 105])

title('WC-DS', 'fontsize', fontsize)



% % uncertainty within the scores from each participant 

leftmargin = 1.5;
rightmargin = 13;
bottommargin = 19;
topmargin = 9.5;
noc = 1;
nor = 1;
spacer = 0;
spacec = 0;
positions = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, spacer, spacec);

axes('position',positions{1},'XGrid','off','XMinorGrid','off','FontSize',8,'Box','off','Layer','top', 'linewidth', 1.5);


currScores = scores_all.(subjects{ii});
currScores = currScores(~isnan(currScores(:,3)), 1:2);

allDiff = scores_all.(subjects{ii})(:, 3:4); 

[~, sortIdx] = sort(currScores(:, 1), 'ascend');

imagesc(currScores(sortIdx, :)')
colormap('jet')

set(gca, 'fontsize', fontsize, 'linewidth', 1.5, 'YTick', [])


    title('consistency of scoring', 'fontsize', fontsize)

    set(gca, 'ytick', [1 2],'yticklabel', {'1st'; '2nd'})



xlabel('PBE no.', 'fontsize', fontsize) 

caxis([1 8])


% the histogram of the differences between the two repetitions 
leftmargin = 1.5;
rightmargin = 13;
bottommargin = 14;
topmargin = 11;
noc = 1;
nor = 1;
spacer = 0;
spacec = 0;
positions = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, spacer, spacec);

axes('position', positions{1},'XGrid','off','XMinorGrid','off','FontSize',8,'Box','off','Layer','top', 'linewidth', 1.5);

allDiff = allDiff(~isnan(allDiff(:, 1)), :);


count = zeros(max(allDiff(:,1)), 1);
for jj = 1:max(allDiff(:,1))+1
  
    count(jj) = numel(find(allDiff(:, 1) == jj-1));
        
end

count = count/numel(allDiff(:, 1));

bar(0:max(allDiff(:, 1)), count, 'FaceColor', 'k', 'FaceAlpha', 0.7, 'EdgeColor', 'none')

xlabel('difference btw rep.s', 'fontsize', fontsize)
ylabel('ratio of PBEs', 'fontsize', fontsize)

set(gca, 'box', 'off', 'linewidth', 1.5)





% % distribution of algorithmic replay scores for two categories with high
% and low human replay scores



leftmargin = 1.5;
rightmargin = 1;
bottommargin = 4;
topmargin = 18;
noc = 8;
nor = 3;
spacer = 1.2;
spacec = 0.7;

positions = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, spacer, spacec);


thresholds = 5:7;

for ithresh = 1:numel(thresholds)
   
    currThresh = thresholds(ithresh);
    
    axes('position', positions{ithresh, 1},'XGrid','off','XMinorGrid','off','FontSize',8,'Box','off','Layer','top', 'linewidth', 1);

        
    rt_ts_low = rt_ts(scores < currThresh);
    rt_ts_hi  = rt_ts(scores >= currThresh);
    hold on
    try
    customPlot2(rt_ts_low, rt_ts_hi)
    if ithresh == 1
        title({'RT-TS'}, 'fontsize', fontsize)
    end
    ylabel(sprintf('thresh = %d', currThresh), 'fontsize', fontsize, 'fontweight', 'normal')
    catch
        currThresh
        size(rt_ts_low)
    end
    
    if ithresh == 3 
        xlabel('replay Score %', 'fontsize', fontsize)
    end
    
    axes('position', positions{ithresh, 2},'XGrid','off','XMinorGrid','off','FontSize',8,'Box','off','Layer','top', 'linewidth', 1);

    rt_ui_low = rt_ui(scores < currThresh);
    rt_ui_hi  = rt_ui(scores >= currThresh);
    hold on
    
    
    customPlot2(rt_ui_low, rt_ui_hi)
    if ithresh == 1
        title({'RT-UI'}, 'fontsize', fontsize)
    end
    
    
    
    
    axes('position', positions{ithresh, 3},'XGrid','off','XMinorGrid','off','FontSize',8,'Box','off','Layer','top', 'linewidth', 1);
   
    rt_pf_low = rt_pf(scores < currThresh);
    rt_pf_hi  = rt_pf(scores >= currThresh);
    hold on
    try
    customPlot2(rt_pf_low, rt_pf_hi)
    if ithresh == 1
        title({'RT-PF'}, 'fontsize', fontsize)
    end
    catch
        currThresh
        size(rt_pf_low)
    end
    
    
    
    
    axes('position', positions{ithresh, 4},'XGrid','off','XMinorGrid','off','FontSize',8,'Box','off','Layer','top', 'linewidth', 1);

    rt_ds_low = rt_ds(scores < currThresh);
    rt_ds_hi  = rt_ds(scores >= currThresh);
    hold on
    
    try
    customPlot2(rt_ds_low, rt_ds_hi)   
    if ithresh == 1
        title({'RT-DS'}, 'fontsize', fontsize)
    end
    catch
        currThresh
        size(rt_ds_low)
    end
    
    
    
    axes('position', positions{ithresh, 5},'XGrid','off','XMinorGrid','off','FontSize',8,'Box','off','Layer','top', 'linewidth', 1);

    wc_ts_low = wc_ts(scores < currThresh);
    wc_ts_hi  = wc_ts(scores >= currThresh);
    hold on
    try
    customPlot2(wc_ts_low, wc_ts_hi)
    if ithresh == 1
        title({'WC-TS'}, 'fontsize', fontsize)
    end
    catch
        currThresh
        size(wc_ts_low)
    end
    
    
    
    
    axes('position', positions{ithresh, 6},'XGrid','off','XMinorGrid','off','FontSize',8,'Box','off','Layer','top', 'linewidth', 1);

    wc_ui_low = wc_ui(scores < currThresh);
    wc_ui_hi  = wc_ui(scores >= currThresh);
    hold on
    try
    customPlot2(wc_ui_low, wc_ui_hi)
    if ithresh == 1
        title({'WC-UI'}, 'fontsize', fontsize)
    end

    catch
        currThresh
        size(wc_ui_low)
    end
    

    
    axes('position', positions{ithresh, 7},'XGrid','off','XMinorGrid','off','FontSize',8,'Box','off','Layer','top', 'linewidth', 1);

    wc_pf_low = wc_pf(scores < currThresh);
    wc_pf_hi  = wc_pf(scores >= currThresh);
    hold on
    try
    customPlot2(wc_pf_low, wc_pf_hi)
    if ithresh == 1
        title({'WC-PF'}, 'fontsize', fontsize)
    end
    catch
        size(wc_pf_low)
    end
    
    axes('position', positions{ithresh, 8},'XGrid','off','XMinorGrid','off','FontSize',8,'Box','off','Layer','top', 'linewidth', 1);

    wc_ds_low = wc_ds(scores < currThresh);
    wc_ds_hi  = wc_ds(scores >= currThresh);
    hold on
    try
    customPlot2(wc_ds_low, wc_ds_hi)
    if ithresh == 1
        title({'WC-DS'}, 'fontsize', fontsize)
    end

    catch
        currThresh
        size(wc_ds_low)
    end
    
    
end


% % diff. btw high and low logarithmic replay scores vs threshold


leftmargin = 9;
rightmargin = 1.5;
bottommargin = 14;
topmargin = 10;
noc = 1;
nor = 1;
spacer = 0;
spacec = 0;

positions = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, spacer, spacec);

axes('position', positions{1},'XGrid','off','XMinorGrid','off','FontSize',8,'Box','off','Layer','top', 'linewidth', 1.5);
hold on

thresholds = 2:8;
nThreshs = numel(thresholds);

rt_ts_diff = zeros(nThreshs, 1);
rt_ui_diff = zeros(nThreshs, 1);
wc_ts_diff = zeros(nThreshs, 1);
wc_ui_diff = zeros(nThreshs, 1);

rt_pf_diff = zeros(nThreshs, 1);
rt_ds_diff = zeros(nThreshs, 1);
wc_pf_diff = zeros(nThreshs, 1);
wc_ds_diff = zeros(nThreshs, 1);


incluPeriodIdx = PRE==1| RUN==1 | POST==1;

for ithresh = 1:nThreshs
   
    currThresh = thresholds(ithresh);
    
    rt_ts_diff(ithresh) = nanmedian(rt_ts(scores >= currThresh & incluPeriodIdx)) - nanmedian(rt_ts(scores < currThresh & incluPeriodIdx));
    rt_ui_diff(ithresh) = nanmedian(rt_ui(scores >= currThresh & incluPeriodIdx)) - nanmedian(rt_ui(scores < currThresh & incluPeriodIdx));
    wc_ts_diff(ithresh) = nanmedian(wc_ts(scores >= currThresh & incluPeriodIdx)) - nanmedian(wc_ts(scores < currThresh & incluPeriodIdx));
    wc_ui_diff(ithresh) = nanmedian(wc_ui(scores >= currThresh & incluPeriodIdx)) - nanmedian(wc_ui(scores < currThresh & incluPeriodIdx));
    
    rt_pf_diff(ithresh) = nanmedian(rt_pf(scores >= currThresh & incluPeriodIdx)) - nanmedian(rt_pf(scores < currThresh & incluPeriodIdx));
    rt_ds_diff(ithresh) = nanmedian(rt_ds(scores >= currThresh & incluPeriodIdx)) - nanmedian(rt_ds(scores < currThresh & incluPeriodIdx));
    wc_pf_diff(ithresh) = nanmedian(wc_pf(scores >= currThresh & incluPeriodIdx)) - nanmedian(wc_pf(scores < currThresh & incluPeriodIdx));
    wc_ds_diff(ithresh) = nanmedian(wc_ds(scores >= currThresh & incluPeriodIdx)) - nanmedian(wc_ds(scores < currThresh & incluPeriodIdx));

end



plot(thresholds, rt_ts_diff, 'color', [204 0 0]/255, 'linewidth', 2)

p = plot(thresholds, rt_ui_diff, 'color', [204 0 0]/255, 'linewidth', 2);
p.Color(4) = 0.5;

plot(thresholds, rt_pf_diff, 'color', [204 0 0]/255, 'linewidth', 2, 'linestyle', ':')

p = plot(thresholds, rt_ds_diff, 'color', [204 0 0]/255, 'linewidth', 2, 'linestyle', ':');
p.Color(4) = 0.5;


%

plot(thresholds, wc_ts_diff, 'color', [0 153 0]/255, 'linewidth', 2)

p = plot(thresholds, wc_ui_diff, 'color', [0 153 0]/255, 'linewidth', 2);
p.Color(4) = 0.5;


plot(thresholds, wc_pf_diff, 'color', [0 153 0]/255, 'linewidth', 2, 'linestyle', ':')

p = plot(thresholds, wc_ds_diff, 'color', [0 153 0]/255, 'linewidth', 2, 'linestyle', ':');
p.Color(4) = 0.5;



xlabel('human score thresh.', 'fontsize', fontsize)
ylabel('median diff. in replay score (high-low)', 'fontsize', fontsize)
legend('RT-TS', 'RT-UI', 'RT-PF', 'RT-DS', 'WC-TS', 'WC-UI', 'WC-PF', 'WC-DS', 'fontsize', 7, 'location', 'eastout')

grid on

xlim([2 8])


print(gcf, sprintf('plots_%s', currSubject), '-dpdf', '-r0')

end

%%
% 
% figure;
% hold on
% 
% subjects = {'KM'; 'BG'; 'NK'; 'JA'; 'UK'};
% colors = {'b'; 'g';'r'; 'c'; 'm'};
% 
% for ss = 1:numel(subjects)
%     
%     subjectNames = subjects{ss};
%     scores = scores_all.(subjectNames)(:, 5);
%     
%     includIdx = find(scores >= 0);
%     
%     r_rt_ts.(subjectNames) = corr(scores(includIdx), rt_ts(includIdx), 'type', 'pearson'); 
%     r_rt_ui.(subjectNames) = corr(scores(includIdx), rt_ui(includIdx), 'type', 'pearson');
%     r_rt_pf.(subjectNames) = corr(scores(includIdx), rt_pf(includIdx), 'type', 'pearson'); 
%     r_rt_ds.(subjectNames) = corr(scores(includIdx), rt_ds(includIdx), 'type', 'pearson');
%     
%     r_wc_ts.(subjectNames) = corr(scores(includIdx), wc_ts(includIdx), 'type', 'pearson');
%     r_wc_ui.(subjectNames) = corr(scores(includIdx), wc_ui(includIdx), 'type', 'pearson');
%     r_wc_pf.(subjectNames) = corr(scores(includIdx), wc_pf(includIdx), 'type', 'pearson');
%     r_wc_ds.(subjectNames) = corr(scores(includIdx), wc_ds(includIdx), 'type', 'pearson');
%     
%     
%     scatter(1, r_rt_ts.(subjectNames), 10, colors{ss}, 'filled')
%     scatter(2, r_rt_ui.(subjectNames), 10, colors{ss}, 'filled')
%     scatter(3, r_rt_pf.(subjectNames), 10, colors{ss}, 'filled')
%     scatter(4, r_rt_ds.(subjectNames), 10, colors{ss}, 'filled')
%     scatter(5, r_wc_ts.(subjectNames), 10, colors{ss}, 'filled')
%     scatter(6, r_wc_ui.(subjectNames), 10, colors{ss}, 'filled')
%     scatter(7, r_wc_pf.(subjectNames), 10, colors{ss}, 'filled')
%     scatter(8, r_wc_ds.(subjectNames), 10, colors{ss}, 'filled')
%     
% end


%% functions


function customPlot(lowDist, hiDist)

binEdges = linspace(-0.01, 100.01, 10);
binCenters = binEdges(1:end-1)+diff(binEdges(1:2))/2;

h_low = histc(lowDist, binEdges); h_low(end) = [];
h_low = h_low/sum(h_low);
low_median = nanmedian(lowDist);

h_hi  = histc(hiDist,  binEdges); h_hi(end) = [];
h_hi  = h_hi/sum(h_hi);
hi_median  = nanmedian(hiDist);

% pval = ranksum(hiDist, lowDist, 'tail', 'right');
difference = hi_median - low_median;


hold on
bar(binCenters, h_low, 'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.2)
bar(binCenters, h_hi, 'FaceColor', [32, 131, 188]/255, 'EdgeColor', 'none', 'FaceAlpha', 0.5)
yl = ylim(gca);

line([low_median low_median], [yl(1) yl(2)], 'linewidth', 1.5, 'color', 'k')
line([hi_median hi_median], [yl(1) yl(2)], 'linewidth', 1.5, 'color', [32, 131, 188]/255)

text(0, yl(2)*0.6, sprintf('diff=%.1f', difference), 'fontsize', 7, 'color', [0.2 0.2 0.2])

% xlim([0 100])
ylim([0 1.1*yl(2)])

% try
%     xticks([low_median hi_median])
% catch
%     xticks([hi_median low_median])
% end

xticks([0 100])
% xtickangle(45)



end


function customPlot2(lowDist, hiDist)


binEdges = linspace(-0.01, 100.01, 100);
binCenters = binEdges(1:end-1)+diff(binEdges(1:2))/2;

% within-PBE time swap
h_low  = calpdf(lowDist, binEdges);
low_median = nanmedian(lowDist);


h_hi   = calpdf(hiDist, binEdges);
hi_median  = nanmedian(hiDist);

difference = hi_median - low_median;



hold on

area(binCenters, h_low, 'facecolor', 'k', 'edgecolor', 'none', 'facealpha', 0.2)
area(binCenters, h_hi, 'facecolor', 'b', 'edgecolor', 'none', 'facealpha', 0.3)


% plot(binCenters, h_low, 'color', 'k', 'linestyle', '-', 'linewidth', 1.5);
% plot(binCenters, h_hi, 'color', [32, 131, 188]/255, 'linestyle', '-', 'linewidth', 1.5);


yl = ylim(gca);

line([low_median low_median], [yl(1) yl(2)], 'linewidth', 1, 'color', 'k')
line([hi_median hi_median], [yl(1) yl(2)], 'linewidth', 1, 'color', 'b')

text(5, yl(2)*0.6, sprintf('diff=%.1f', difference), 'fontsize', 7, 'color', [0.2 0.2 0.2])

% xlim([0 100])
ylim([0 1.1*yl(2)])
xticks([0 100])


end

function [scores2, uniqPBEidx] = scoreConsistency(scores, PBEno)

% 
uniqs = unique(PBEno);
uniqs = sort(uniqs, 'ascend');

uniqPBEs  = nan(numel(uniqs), 2);
scores2   = nan(numel(uniqs), 2);

uniqPBEidx = zeros(numel(uniqs), 1);

for ii = 1:numel(uniqs)

   temp = find(PBEno == uniqs(ii));
   uniqPBEs(ii, 1:numel(temp)) = temp;
   scores2(ii, 1:numel(temp))  = scores(temp);
   
   uniqPBEidx(ii) = find(PBEno == uniqs(ii), 1, 'first');

end


end


function h = calpdf(seqScore, bins)

h = histc(seqScore, bins);
h(end-1) = h(end-1)+ h(end);
h(end) = [];

h = h/sum(h);
h = conv(h, ones(10,1)/10, 'same');
% h = cumsum(h);

end










