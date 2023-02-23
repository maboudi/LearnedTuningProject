
clear

load('rndPBEs_new.mat')
load('rndPBEs.mat')


% idx = [];
for ii = 1:length(rndPBEs_new)
    
    rndPBEs_new(ii).PBEno = rndPBEs(ii).PBEno;
%     if ~isnan(rndPBEs2(ii).startT)
%         idx = [idx; ii];
%     end
    
end

rndPBEs = rndPBEs_new;%(idx);


subjects = {'KM'; 'BG'; 'NK'; 'JA'; 'UK'; 'PH'; 'EP'; 'KD'};


for ss = 1: numel(subjects)
    
   load(sprintf('HumanScores_%s', subjects{ss}), 'scores')
   
%    if strcmp(subjects{ss}, 'NK')
%        scores(scores == 4) = 0;
%    end
   
%    scores = scores(idx);
   

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

rndPBEs = rndPBEs(uniqPBEidx); % get rid of the repeated PBEs
scores = [];

rndPBEs_all = [];
for ss = 1:numel(subjects)
    
    currentScores = scores_all.(subjects{ss})(:, 5);
%     includ = find(currentScores >= 1);
    
    scores = [scores; currentScores];
    rndPBEs_all = [rndPBEs_all; rndPBEs(:)];
end



rndPBEs = rndPBEs_all;


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


%%

figure;
set(gca, 'linewidth', 2, 'fontsize', 12)

binEdges = 0:9;
h_PRE = histc(scores(PRE==1), binEdges);
h_RUN = histc(scores(RUN==1), binEdges);
h_POST = histc(scores(POST==1), binEdges);

bar(binEdges, [h_PRE h_RUN h_POST], 'stacked')

legend('PRE', 'RUN', 'POST')
xlabel('human replay score', 'fontsize', 12)
ylabel('number of PBEs', 'fontsize', 12)




[r_rt_ts, p_rt_ts] = corr(scores, rt_ts, 'type', 'spearman'); 
[r_rt_ui, p_rt_ui] = corr(scores, rt_ui, 'type', 'spearman'); 
[r_wc_ts, p_wc_ts] = corr(scores, wc_ts, 'type', 'spearman');
[r_wc_ui, p_wc_ui] = corr(scores, wc_ui, 'type', 'spearman');


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
    


figure;

set(gca, 'position', [400 800 1500 300])

subplot(2,4,1)
set(gca, 'linewidth', 1.5)
hold on
scatter(scores+randn(size(scores))*0.1, rt_ts, 5, 'filled')
scatter(0:9, rt_ts_med, 15, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'Marker', 'd')

% xlim([-1 10])
% ylim([-5 105])

ylabel('replay score %', 'fontsize', 12)
xlabel('human replay score', 'fontsize', 12)
title('RT-TS', 'fontsize', 12)

subplot(2,4,2)
set(gca, 'linewidth', 1.5)

hold on
scatter(scores+randn(size(scores))*0.1, rt_ui, 5, 'filled')
scatter(0:9, rt_ui_med, 15, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'Marker', 'd')

xlim([-1 10])
ylim([-5 105])

title('RT-UI', 'fontsize', 12)

subplot(2,4,3)
set(gca, 'linewidth', 1.5)

hold on
scatter(scores+randn(size(scores))*0.1, wc_ts, 5, 'filled')
scatter(0:9, wc_ts_med, 15, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'Marker', 'd')
xlim([-1 10])
ylim([-5 105])

title('WC-TS', 'fontsize', 12)

subplot(2,4,4)
set(gca, 'linewidth', 1.5)

hold on
scatter(scores+randn(size(scores))*0.1, wc_ui, 5, 'filled')
scatter(0:9, wc_ui_med, 15, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'Marker', 'd')
xlim([-1 10])
ylim([-5 105])

title('WC-UI', 'fontsize', 12)


%

subplot(2,4,5)
set(gca, 'linewidth', 1.5)

hold on
scatter(scores+randn(size(scores))*0.1, rt_pf, 5, 'filled')
scatter(0:9, rt_pf_med, 15, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'Marker', 'd')
xlim([-1 10])
ylim([-5 105])

ylabel('replay score %', 'fontsize', 12)
xlabel('human replay score', 'fontsize', 12)
title('RT-PF', 'fontsize', 12)

subplot(2,4,6)
set(gca, 'linewidth', 1.5)

hold on
scatter(scores+randn(size(scores))*0.1, rt_ds, 5, 'filled')
scatter(0:9, rt_ds_med, 15, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'Marker', 'd')
xlim([-1 10])
ylim([-5 105])

title('RT-DS', 'fontsize', 12)

subplot(2,4,7)
set(gca, 'linewidth', 1.5)

hold on
scatter(scores+randn(size(scores))*0.1, wc_pf, 5, 'filled')
scatter(0:9, wc_pf_med, 15, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'Marker', 'd')
xlim([-1 10])
ylim([-5 105])

title('WC-PF', 'fontsize', 12)

subplot(2,4,8)
set(gca, 'linewidth', 1.5)

hold on
scatter(scores+randn(size(scores))*0.1, wc_ds, 5, 'filled')
scatter(0:9, wc_ds_med, 15, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'Marker', 'd')
xlim([-1 10])
ylim([-5 105])

title('WC-DS', 'fontsize', 12)

%% divergence between the scores across participants

referenceScores = scores_all.KM(:, 5);

[~, sortIdx] = sort(referenceScores, 'ascend');


sortedScores = zeros(numel(subjects), numel(referenceScores));

for ss = 1:numel(subjects)
    
    currSubject = subjects{ss};
    
    currScores = scores_all.(subjects{ss})(sortIdx, 5);
    
    sortedScores(ss, :) = currScores;
end




cm = colormap('jet');
cm(1, :) = [1 1 1];

scoreCorr = zeros(numel(subjects));
for ss = 1:numel(subjects)
    for uu = 1:numel(subjects) 
        
        if uu >= ss; continue; end
        scores1 = scores_all.(subjects{ss})(:, 5);
        scores2 = scores_all.(subjects{uu})(:, 5);
                
        scoreCorr(ss, uu) = corr(scores1, scores2, 'type', 'spearman');
        
    end
end


figure; 

ax1 = subplot(2,8,1:8);
imagesc(sortedScores)
set(gca, 'fontsize', 12, 'linewidth', 2, 'YTick', 1:numel(subjects), 'YTickLabel', subjects)
ylabel('participants', 'fontsize', 14);
xlabel('PBE no.', 'fontsize', 14)

colormap(ax1, 'jet');

cb = colorbar;
cb.LineWidth = 2;
cb.FontSize = 12;
cb.Label.String = 'human replay score';
cb.Ticks = 1:9;


ax2 = subplot(2,8, 9);
imagesc(scoreCorr(2:numel(subjects), 1:numel(subjects)-1))
set(gca, 'fontsize', 12, 'linewidth', 2, 'YTick', 1:numel(subjects)-1, 'YTickLabel', subjects(2:numel(subjects)), 'XTick', 1:numel(subjects)-1, 'XTickLabel', subjects(1:numel(subjects)-1), 'box', 'off')

colormap(ax2, cm);
cb2 = colorbar;

cb2 = colorbar;
cb2.LineWidth = 2;
cb2.FontSize = 12;
cb2.Label.String = 'correlation coeff.';
cb2.Ticks = 1:9;



%% uncertainty within each subject

figure;


allDiff = [];

for ss = 1:numel(subjects)
   
    currScores = scores_all.(subjects{ss});
    currScores = currScores(~isnan(currScores(:,3)), 1:2);
    
    allDiff = [allDiff; scores_all.(subjects{ss})(:, 3:4)]; 
    
    [~, sortIdx] = sort(currScores(:, 1), 'ascend');
    
    subplot(numel(subjects)+4,1,ss)
    
    imagesc(currScores(sortIdx, :)')
    colormap('jet')
    
    set(gca, 'fontsize', 12, 'linewidth', 2, 'YTick', [])
    
    if ss == 1
       
        title('consistency of scoring', 'fontsize', 14)
        
        set(gca, 'ytick', [1 2],'yticklabel', {'first time'; 'second time'})
    end
    
    
    if ss == numel(subjects)
        xlabel('PBE no.', 'fontsize', 12) 
    end
    
    
    if ss == ceil(numel(subjects)/2)
        ylabel('participants', 'fontsize', 12)
    end
    caxis([1 8])
end

allDiff = allDiff(~isnan(allDiff(:, 1)), :);

subplot(numel(subjects)+4, 1, [ss+2:ss+4])


count = zeros(max(allDiff(:,1)), 1);
for ii = 1:max(allDiff(:,1))+1
  
    count(ii) = numel(find(allDiff(:, 1) == ii-1));
        
end

count = count/numel(allDiff(:, 1));

bar(0:max(allDiff(:, 1)), count, 'FaceColor', 'k', 'FaceAlpha', 0.8, 'EdgeColor', 'none')


xlabel('difference bw repetitions', 'fontsize', 14)
ylabel('ratio of PBEs x participants', 'fontsize', 14)

set(gca, 'linewidth', 2)
% 
% subplot(numel(subjects)+7, 1, [ss+5:ss+7])
% hist(allDiff(:, 2));
% 
% 
% xlabel('difference bw repitions %', 'fontsize', 14)
% ylabel('ratio of PBEs x participants', 'fontsize', 14)
% 
% set(gca, 'linewidth', 2)

%% plot distributions of replay scores by applying thresholds on human replay scores

thresholds = 5:7;

figure;
set(gcf, 'position', [300 300 1105 550])
for ithresh = 1:numel(thresholds)
   
    currThresh = thresholds(ithresh);
    
    subplot(8, numel(thresholds), ithresh)
    set(gca, 'linewidth', 1.5)
        
    rt_ts_low = rt_ts(scores < currThresh);
    rt_ts_hi  = rt_ts(scores >= currThresh);
    hold on
    try
    customPlot(rt_ts_low, rt_ts_hi)
    if ithresh == 1
        ylabel({'RT-TS'; ''; 'ratio of PBEs'}, 'fontsize', 12)
    end
    title(sprintf('thresh = %d', currThresh), 'fontsize', 14, 'fontweight', 'normal')
    catch
        currThresh
        size(rt_ts_low)
    end
    
    subplot(8, numel(thresholds), ithresh+numel(thresholds))
    set(gca, 'linewidth', 1.5)

    rt_ui_low = rt_ui(scores < currThresh);
    rt_ui_hi  = rt_ui(scores >= currThresh);
    hold on
    
    try
    customPlot(rt_ui_low, rt_ui_hi)
    if ithresh == 1
        ylabel({'RT-UI'; ''; 'ratio of PBEs'}, 'fontsize', 12)
    end
    catch
        currThresh
        size(rt_ui_low)
    end
    
    
    
    
    subplot(8, numel(thresholds), ithresh+numel(thresholds)*2)
    set(gca, 'linewidth', 1.5)
        
    rt_pf_low = rt_pf(scores < currThresh);
    rt_pf_hi  = rt_pf(scores >= currThresh);
    hold on
    try
    customPlot(rt_pf_low, rt_pf_hi)
    if ithresh == 1
        ylabel({'RT-PF'; ''; 'ratio of PBEs'}, 'fontsize', 12)
    end
    catch
        currThresh
        size(rt_pf_low)
    end
    
    subplot(8, numel(thresholds), ithresh+numel(thresholds)*3)
    set(gca, 'linewidth', 1.5)

    rt_ds_low = rt_ds(scores < currThresh);
    rt_ds_hi  = rt_ds(scores >= currThresh);
    hold on
    
    try
    customPlot(rt_ds_low, rt_ds_hi)   
    if ithresh == 1
        ylabel({'RT-DS'; ''; 'ratio of PBEs'}, 'fontsize', 12)
    end
    catch
        currThresh
        size(rt_ds_low)
    end
    
    
    
    subplot(8, numel(thresholds), ithresh+numel(thresholds)*4)
    set(gca, 'linewidth', 1.5)

    wc_ts_low = wc_ts(scores < currThresh);
    wc_ts_hi  = wc_ts(scores >= currThresh);
    hold on
    try
    customPlot(wc_ts_low, wc_ts_hi)
    if ithresh == 1
        ylabel({'WC-TS'; ''; 'ratio of PBEs'}, 'fontsize', 12)
    end
    catch
        currThresh
        size(wc_ts_low)
    end
    
    subplot(8, numel(thresholds), ithresh+numel(thresholds)*5)
    set(gca, 'linewidth', 1.5)

    wc_ui_low = wc_ui(scores < currThresh);
    wc_ui_hi  = wc_ui(scores >= currThresh);
    hold on
    try
    customPlot(wc_ui_low, wc_ui_hi)
    if ithresh == 1
        ylabel({'WC-UI'; ''; 'ratio of PBEs'}, 'fontsize', 12)
    end

%     xlabel('replay Score %', 'fontsize', 12)
    catch
        currThresh
        size(wc_ui_low)
    end
    
    %
    
    subplot(8, numel(thresholds), ithresh+numel(thresholds)*6)
    set(gca, 'linewidth', 1.5)

    wc_pf_low = wc_pf(scores < currThresh);
    wc_pf_hi  = wc_pf(scores >= currThresh);
    hold on
    try
    customPlot(wc_pf_low, wc_pf_hi)
    if ithresh == 1
        ylabel({'WC-PF'; ''; 'ratio of PBEs'}, 'fontsize', 12)
    end
    catch
        size(wc_pf_low)
    end
    
    subplot(8, numel(thresholds), ithresh+numel(thresholds)*7)
    set(gca, 'linewidth', 1.5)

    wc_ds_low = wc_ds(scores < currThresh);
    wc_ds_hi  = wc_ds(scores >= currThresh);
    hold on
    try
    customPlot(wc_ds_low, wc_ds_hi)
    if ithresh == 1
        ylabel({'WC-DS'; ''; 'ratio of PBEs'}, 'fontsize', 12)
    end

    xlabel('replay Score %', 'fontsize', 12)
    catch
        currThresh
        size(wc_ds_low)
    end
    
    
end


%% correlation between the median difference and threshold on human scores


thresholds = 2:9;
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


figure;
hold on
set(gca, 'linewidth', 1.5)


plot(thresholds, rt_ts_diff, 'color', [204 0 0]/255, 'linewidth', 3)

p = plot(thresholds, rt_ui_diff, 'color', [204 0 0]/255, 'linewidth', 3);
p.Color(4) = 0.5;

plot(thresholds, rt_pf_diff, 'color', [204 0 0]/255, 'linewidth', 3, 'linestyle', ':')

p = plot(thresholds, rt_ds_diff, 'color', [204 0 0]/255, 'linewidth', 3, 'linestyle', ':');
p.Color(4) = 0.5;


%

plot(thresholds, wc_ts_diff, 'color', [0 153 0]/255, 'linewidth', 3)

p = plot(thresholds, wc_ui_diff, 'color', [0 153 0]/255, 'linewidth', 3);
p.Color(4) = 0.5;


plot(thresholds, wc_pf_diff, 'color', [0 153 0]/255, 'linewidth', 3, 'linestyle', ':')

p = plot(thresholds, wc_ds_diff, 'color', [0 153 0]/255, 'linewidth', 3, 'linestyle', ':');
p.Color(4) = 0.5;



xlabel('human score thresh.', 'fontsize', 12)
ylabel('median diff. in replay score (high-low)', 'fontsize', 12)
legend('RT-TS', 'RT-UI', 'RT-PF', 'RT-DS', 'WC-TS', 'WC-UI', 'WC-PF', 'WC-DS', 'fontsize', 10, 'location', 'best')




%% functions
function customPlot(lowDist, hiDist)

binEdges = linspace(-0.01, 100.01, 20);
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
bar(binCenters, h_low, 'FaceColor', 'none', 'EdgeColor', 'k')
bar(binCenters, h_hi, 'FaceColor', [32, 131, 188]/255, 'EdgeColor', 'none', 'FaceAlpha', 0.7)
yl = ylim(gca);

line([low_median low_median], [yl(1) yl(2)], 'linewidth', 3, 'color', 'k')
line([hi_median hi_median], [yl(1) yl(2)], 'linewidth', 3, 'color', [32, 131, 188]/255)

text(10, yl(2)*0.4, sprintf('diff = %.1f', difference), 'fontsize', 12, 'color', [0.2 0.2 0.2])

% xlim([0 100])
ylim([0 1.1*yl(2)])

try
    xticks([low_median hi_median])
catch
    xticks([hi_median low_median])
end
xtickangle(90)



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
