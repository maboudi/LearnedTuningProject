
rndPBEs_all = [];
load('rndPBEs.mat')
rndPBEs_all = [rndPBEs_all; rndPBEs(:)];

load('rndPBEs.mat')
rndPBEs_all = [rndPBEs_all; rndPBEs(:)];

rndPBEs = rndPBEs_all;

hscores = [];
load('HumanScores_KouroshNK.mat')
scores(scores == 4) = 0;
hscores = [hscores; scores];
% cmpr(:,1) = scores;

load('HumanScores_Kourosh.mat')

hscores = [hscores; scores];
% cmpr(:,2) = scores;

% cmpr = cmpr + randn(size(cmpr))*0.3;


% 
% indx = [rndPBEs.PBEno]';
% 
% uniqs = unique(indx);
% uniqs = sort(uniqs, 'ascend');
% 
% uniqPBEs = nan(numel(uniqs), 2);
% scores = nan(numel(uniqs), 2);
% 
% for ii = 1:numel(uniqs)
%     
%     temp = find(indx == uniqs(ii));
%    uniqPBEs(ii, 1:numel(temp)) = temp;
%    scores(ii, 1:numel(temp))   = hscores(temp);
%    
% end
% 
% hscores  = scores(:,1);
% rndPBEs  = rndPBEs(uniqPBEs(:,1));



rt_ts = [rndPBEs.rt_ts]';
rt_ui = [rndPBEs.rt_ui]';
wc_ts = [rndPBEs.wc_ts]';
wc_ui = [rndPBEs.wc_ui]';


PRE  = [rndPBEs.PRE]';
RUN  = [rndPBEs.RUN]';
POST = [rndPBEs.POST]';


% hscores = scores;

figure;
set(gca, 'linewidth', 2, 'fontsize', 12)

binEdges = 0:9;

h = histc(hscores, binEdges);

bar(binEdges, h, 'FaceColor', [.2 .2 .2], 'EdgeColor', 'none')
xlabel('human replay score', 'fontsize', 12)
ylabel('number of PBEs', 'fontsize', 12)


figure;

binEdges = 0:9;
h_PRE = histc(hscores(PRE==1), binEdges);
h_RUN = histc(hscores(RUN==1), binEdges);
h_POST = histc(hscores(POST==1), binEdges);

bar(binEdges, [h_PRE h_RUN h_POST], 'stacked')

legend('PRE', 'RUN', 'POST')
xlabel('human replay score', 'fontsize', 12)
ylabel('number of PBEs', 'fontsize', 12)




[r_rt_ts, p_rt_ts] = corr(hscores, rt_ts, 'type', 'spearman'); 
[r_rt_ui, p_rt_ui] = corr(hscores, rt_ui, 'type', 'spearman'); 
[r_wc_ts, p_wc_ts] = corr(hscores, wc_ts, 'type', 'spearman');
[r_wc_ui, p_wc_ui] = corr(hscores, wc_ui, 'type', 'spearman');


rt_ts_med = zeros(10,1);
rt_ui_med = zeros(10,1);
wc_ts_med = zeros(10,1);
wc_ui_med = zeros(10,1);

for ihscore = 1:10
    
    currhscore = ihscore-1;
    rt_ts_med(ihscore) = nanmedian(rt_ts(hscores == currhscore));
    rt_ui_med(ihscore) = nanmedian(rt_ui(hscores == currhscore));
    wc_ts_med(ihscore) = nanmedian(wc_ts(hscores == currhscore));
    wc_ui_med(ihscore) = nanmedian(wc_ui(hscores == currhscore));
end
    


figure;

set(gca, 'position', [400 800 1500 300])

subplot(1,4,1)
set(gca, 'linewidth', 1.5)

hold on
scatter(hscores+randn(size(hscores))*0.1, rt_ts, 5, 'filled')
scatter(0:9, rt_ts_med, 15, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'Marker', 'd')
xlim([-1 10])
ylim([-5 105])

ylabel('replay score %', 'fontsize', 12)
xlabel('human replay score', 'fontsize', 12)
title('RT-TS', 'fontsize', 12)

subplot(1,4,2)
set(gca, 'linewidth', 1.5)

hold on
scatter(hscores+randn(size(hscores))*0.1, rt_ui, 5, 'filled')
scatter(0:9, rt_ui_med, 15, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'Marker', 'd')
xlim([-1 10])
ylim([-5 105])

title('RT-UI', 'fontsize', 12)

subplot(1,4,3)
set(gca, 'linewidth', 1.5)

hold on
scatter(hscores+randn(size(hscores))*0.1, wc_ts, 5, 'filled')
scatter(0:9, wc_ts_med, 15, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'Marker', 'd')
xlim([-1 10])
ylim([-5 105])

title('WC-TS', 'fontsize', 12)

subplot(1,4,4)
set(gca, 'linewidth', 1.5)

hold on
scatter(hscores+randn(size(hscores))*0.1, wc_ui, 5, 'filled')
scatter(0:9, wc_ui_med, 15, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'Marker', 'd')
xlim([-1 10])
ylim([-5 105])

title('WC-UI', 'fontsize', 12)




%% plot distributions of replay scores by applying thresholds on human replay scores

thresholds = 5:7;

figure;
set(gcf, 'position', [300 300 1105 550])
for ithresh = 1:numel(thresholds)
   
    currThresh = thresholds(ithresh);
    
    subplot(4, numel(thresholds), ithresh)
    set(gca, 'linewidth', 1.5)

    rt_ts_low = rt_ts(hscores < currThresh);
    rt_ts_hi  = rt_ts(hscores >= currThresh);
    hold on
    try
    customPlot(rt_ts_low, rt_ts_hi)
    if ithresh == 1
        ylabel({'RT-TS'; ''; 'ratio of PBEs'}, 'fontsize', 12)
    end
    title(sprintf('thresh = %d', currThresh), 'fontsize', 14, 'fontweight', 'normal')
    catch
        size(rt_ts_low)
    end
    
    subplot(4, numel(thresholds), ithresh+numel(thresholds))
    set(gca, 'linewidth', 1.5)

    rt_ui_low = rt_ui(hscores < currThresh);
    rt_ui_hi  = rt_ui(hscores >= currThresh);
    hold on
    
    try
    customPlot(rt_ui_low, rt_ui_hi)
    if ithresh == 1
        ylabel({'RT-UI'; ''; 'ratio of PBEs'}, 'fontsize', 12)
    end
    catch
        size(rt_ui_low)
    end
    
    subplot(4, numel(thresholds), ithresh+numel(thresholds)*2)
    set(gca, 'linewidth', 1.5)

    wc_ts_low = wc_ts(hscores < currThresh);
    wc_ts_hi  = wc_ts(hscores >= currThresh);
    hold on
    try
    customPlot(wc_ts_low, wc_ts_hi)
    if ithresh == 1
        ylabel({'WC-TS'; ''; 'ratio of PBEs'}, 'fontsize', 12)
    end
    catch
        size(wc_ts_low)
    end
    
    subplot(4, numel(thresholds), ithresh+numel(thresholds)*3)
    set(gca, 'linewidth', 1.5)

    wc_ui_low = wc_ui(hscores < currThresh);
    wc_ui_hi  = wc_ui(hscores >= currThresh);
    hold on
    try
    customPlot(wc_ui_low, wc_ui_hi)
    if ithresh == 1
        ylabel({'WC-UI'; ''; 'ratio of PBEs'}, 'fontsize', 12)
    end

    xlabel('replay Score %', 'fontsize', 12)
    catch
        size(wc_ui_low)
    end
    
end


%% correlation between the median difference and threshold on human scores


thresholds = 2:8;
nThreshs = numel(thresholds);

rt_ts_diff = zeros(nThreshs, 1);
rt_ui_diff = zeros(nThreshs, 1);
wc_ts_diff = zeros(nThreshs, 1);
wc_ui_diff = zeros(nThreshs, 1);

incluPeriodIdx = PRE==1| RUN==1 | POST==1;

for ithresh = 1:nThreshs
   
    currThresh = thresholds(ithresh);
    
    rt_ts_diff(ithresh) = nanmedian(rt_ts(hscores >= currThresh & incluPeriodIdx)) - nanmedian(rt_ts(hscores < currThresh & incluPeriodIdx));
    rt_ui_diff(ithresh) = nanmedian(rt_ui(hscores >= currThresh & incluPeriodIdx)) - nanmedian(rt_ui(hscores < currThresh & incluPeriodIdx));
    wc_ts_diff(ithresh) = nanmedian(wc_ts(hscores >= currThresh & incluPeriodIdx)) - nanmedian(wc_ts(hscores < currThresh & incluPeriodIdx));
    wc_ui_diff(ithresh) = nanmedian(wc_ui(hscores >= currThresh & incluPeriodIdx)) - nanmedian(wc_ui(hscores < currThresh & incluPeriodIdx));

end


figure;
hold on
set(gca, 'linewidth', 1.5)


plot(thresholds, rt_ts_diff, 'color', [204 0 0]/255, 'linewidth', 3)

p = plot(thresholds, rt_ui_diff, 'color', [204 0 0]/255, 'linewidth', 3);
p.Color(4) = 0.5;


plot(thresholds, wc_ts_diff, 'color', [0 153 0]/255, 'linewidth', 3)

p = plot(thresholds, wc_ui_diff, 'color', [0 153 0]/255, 'linewidth', 3);
p.Color(4) = 0.5;


xlabel('human score thresh.', 'fontsize', 12)
ylabel('median diff. in replay score (high-low)', 'fontsize', 12)
legend('RT-TS', 'RT-UI', 'WC-TS', 'WC-UI', 'fontsize', 10, 'location', 'best')




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

xticks([low_median hi_median])
xtickangle(90)



end

