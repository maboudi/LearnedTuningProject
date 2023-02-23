
sessionNames = fieldnames(PBEinfo_pooled);

nSessions  = numel(sessionNames);

rndRUNPBEs   = [];
rndSLEEPPBEs = [];

for ii = 1:nSessions
    
   
    PBEinfo = PBEinfo_pooled.(sessionNames{ii});
        
    runIdx = [PBEinfo.RUN]';

    rippleOverlap = [PBEinfo.rippleOverlap]';
    
    nREM   = [PBEinfo.NREM]';
    
    
    nFiringUnits = [PBEinfo.nFiringUnits]';
    duration     = [PBEinfo.duration]';
    
    rt_ts = [PBEinfo.rt_ts]';
    rt_ui = [PBEinfo.rt_ui]';
    wc_ts = [PBEinfo.wc_ts]';
    wc_ui = [PBEinfo.wc_ui]';
    
    % % run select PBEs
    
%     qualifiedPBEs = find(runIdx == 1 & rippleOverlap == 1 & nFiringUnits >= 5 & duration >= 100); 
    qualifiedPBEs = find(runIdx == 1 & rippleOverlap == 1 & nFiringUnits >= 5 & duration >= 100 & duration <= 600 & (rt_ts > 90 | rt_ui > 90 | wc_ts > 90 | wc_ui> 90));
    
    subset = qualifiedPBEs(randperm(numel(qualifiedPBEs)));
    subset = subset(1:5);
    
    rndRUNPBEs = [rndRUNPBEs PBEinfo(subset)];
    
    
    % % sleep select PBEs
    
    qualifiedPBEs   = find(runIdx ~= 1 & nREM == 1 &  nFiringUnits >= 5 & duration >= 100 & duration <= 600 & (rt_ts > 90 | rt_ui > 90 | wc_ts > 90 | wc_ui > 90));
    
    subset = qualifiedPBEs(randperm(numel(qualifiedPBEs)));
    subset = subset(1:50);
    
    rndSLEEPPBEs = [rndSLEEPPBEs PBEinfo(subset)];
    
end

rndPBEs = [rndRUNPBEs rndSLEEPPBEs];
rndPBEs = rndPBEs(randperm(numel(rndPBEs)));


%% fifty first PBEs are going to be thrown away (expected instability in the beginnig of the scoring)
%% Additionally we repeat 50 PBEs to check stability for the rest of the events.

for ii = 1:550
    rndPBEs(ii).PBEno = ii;
end

firstSet = rndPBEs(1:50);

secondSet = rndPBEs(51:550);
repeated_subset = secondSet(randperm(numel(secondSet), 50));

secondSet_wRep = [secondSet repeated_subset];
secondSet_wRep = secondSet_wRep(randperm(numel(secondSet_wRep)));


rndPBEs = [firstSet secondSet_wRep];


basePath = '~/Documents/HMM_project/PBEs_wRep';
mkdir(basePath)

for ii = 1: numel(rndPBEs)
    
    
    currPPM = rndPBEs(ii).posteriorProbMat;
    [nPosBins nTimeBins] = size(currPPM);
    
    
    maxProb   = max(currPPM(:));
    minProb   = min(currPPM(:));
    Clim_low  = minProb;
    Clim_high = minProb + 0.5 *(maxProb - minProb);
    %%% plot the posterior probability matrix
    figure;
    set(gcf, 'visible', 'off')
    set(gcf, 'paperUnits', 'inches', 'PaperPosition', [0 0 3 3])
    imagesc(1:nTimeBins, (1:nPosBins)*2, currPPM, [Clim_low Clim_high])
    set(gca,'YDir','normal')
    colormap('hot')
    
    xticks([3:3:size(currPPM, 2)])
    yticks(30:30:size(currPPM, 1)*30)
    
    print(gcf, fullfile(basePath, sprintf('pbe_%d', ii)), '-dpng')
    close

end

%%

rt_ts = [rndPBEs.rt_ts]';
rt_ui = [rndPBEs.rt_ui]';
wc_ts = [rndPBEs.wc_ts]';
wc_ui = [rndPBEs.wc_ui]';


figure; 

subplot(2,6,1)
[xx , yy]=ecdf(rt_ts);
plot(yy, xx, 'linewidth', 3, 'color', 'k')
xlabel('radonT-timeSwap', 'fontsize', 12)
ylabel('cdf', 'fontsize', 12)

subplot(2,6,2)
[xx , yy]=ecdf(rt_ui);
plot(yy, xx, 'linewidth', 3, 'color', 'k')
xlabel('radonT-unitID', 'fontsize', 12)
ylabel('cdf', 'fontsize', 12)

subplot(2,6,3)
[xx , yy]=ecdf(wc_ts);
plot(yy, xx, 'linewidth', 3, 'color', 'k')
xlabel('wCorr-timeSwap', 'fontsize', 12)
ylabel('cdf', 'fontsize', 12)

subplot(2,6,4)
[xx , yy]=ecdf(wc_ui);
plot(yy, xx, 'linewidth', 3, 'color', 'k')
xlabel('wCorr-unitID', 'fontsize', 12)
ylabel('cdf', 'fontsize', 12)



subplot(2,6,7)
plot(rt_ts, wc_ts, '.k', 'markersize', 5)
xlabel('radonT-timeSwap', 'fontsize', 12)
ylabel('wCorr-timeSwap', 'fontsize', 12)


subplot(2,6,8)
plot(rt_ui, wc_ui, '.k', 'markersize', 5)
xlabel('radonT-unitID', 'fontsize', 12)
ylabel('wCorr-unitID', 'fontsize', 12)


subplot(2,6,9)
plot(rt_ts, rt_ui, '.k', 'markersize', 5)
xlabel('radonT-timeSwap', 'fontsize', 12)
ylabel('radonT-unitID', 'fontsize', 12)


subplot(2,6,10)
plot(wc_ts, wc_ui, '.k', 'markersize', 5)
xlabel('wCorr-timeSwap', 'fontsize', 12)
ylabel('wCorr-unitID', 'fontsize', 12)


subplot(2,6,11)
plot(wc_ts, rt_ui, '.k', 'markersize', 5)
xlabel('wCorr-timeSwap', 'fontsize', 12)
ylabel('radonT-unitID', 'fontsize', 12)


subplot(2,6,12)
plot(wc_ts, rt_ui, '.k', 'markersize', 5)
xlabel('wCorr-unitID', 'fontsize', 12)
ylabel('radonT-timeSwap', 'fontsize', 12)


%%

runIdx = [rndPBEs.RUN]';
preIdx = [rndPBEs.PRE]';
postIdx = [rndPBEs.POST]';

rt_ts = [rndPBEs.rt_ts]';
rt_ui = [rndPBEs.rt_ui]';
wc_ts = [rndPBEs.wc_ts]';
wc_ui = [rndPBEs.wc_ui]';


figure; 

subplot(2,6,1)
hold on
[xx , yy]=ecdf(rt_ts(runIdx == 1));
p = plot(yy, xx, 'linewidth', 3, 'color', 'k');
p.Color(4) = 0.5;

[xx , yy]=ecdf(rt_ts(preIdx == 1));
p = plot(yy, xx, 'linewidth', 3, 'color', 'b');
p.Color(4) = 0.5;

[xx , yy]=ecdf(rt_ts(postIdx == 1));
p = plot(yy, xx, 'linewidth', 3, 'color', 'r');
p.Color(4) = 0.5;

xlabel('radonT-timeSwap', 'fontsize', 12)
ylabel('cdf', 'fontsize', 12)




subplot(2,6,2)
hold on
[xx , yy]=ecdf(rt_ui(runIdx == 1));
p=plot(yy, xx, 'linewidth', 3, 'color', 'k');
p.Color(4) = 0.5;

[xx , yy]=ecdf(rt_ui(preIdx == 1));
p=plot(yy, xx, 'linewidth', 3, 'color', 'b');
p.Color(4) = 0.5;

[xx , yy]=ecdf(rt_ui(postIdx == 1));
p=plot(yy, xx, 'linewidth', 3, 'color', 'r');
p.Color(4) = 0.5;

xlabel('radonT-unitID', 'fontsize', 12)
ylabel('cdf', 'fontsize', 12)




subplot(2,6,3)
hold on
[xx , yy]=ecdf(wc_ts(runIdx == 1));
p=plot(yy, xx, 'linewidth', 3, 'color', 'k');
p.Color(4) = 0.5;

[xx , yy]=ecdf(wc_ts(preIdx == 1));
p=plot(yy, xx, 'linewidth', 3, 'color', 'b');
p.Color(4) = 0.5;

[xx , yy]=ecdf(wc_ts(postIdx == 1));
p=plot(yy, xx, 'linewidth', 3, 'color', 'r');
p.Color(4) = 0.5;

xlabel('wCorr-timeSwap', 'fontsize', 12)
ylabel('cdf', 'fontsize', 12)




subplot(2,6,4)
hold on
[xx , yy]=ecdf(wc_ui(runIdx == 1));
p=plot(yy, xx, 'linewidth', 3, 'color', 'k');
p.Color(4) = 0.5;

[xx , yy]=ecdf(wc_ui(preIdx == 1));
p=plot(yy, xx, 'linewidth', 3, 'color', 'b');
p.Color(4) = 0.5;

[xx , yy]=ecdf(wc_ui(postIdx == 1));
p=plot(yy, xx, 'linewidth', 3, 'color', 'r');
p.Color(4) = 0.5;

xlabel('wCorr-unitID', 'fontsize', 12)
ylabel('cdf', 'fontsize', 12)






subplot(2,6,7)
hold on

plot(rt_ts(runIdx == 1), wc_ts(runIdx == 1), '.k', 'markersize', 5)
plot(rt_ts(preIdx == 1), wc_ts(preIdx == 1), '.b', 'markersize', 5)
plot(rt_ts(postIdx == 1), wc_ts(postIdx == 1), '.r', 'markersize', 5)

xlabel('radonT-timeSwap', 'fontsize', 12)
ylabel('wCorr-timeSwap', 'fontsize', 12)


subplot(2,6,8)
hold on

plot(rt_ui(runIdx == 1), wc_ui(runIdx == 1), '.k', 'markersize', 5)
plot(rt_ui(preIdx == 1), wc_ui(preIdx == 1), '.b', 'markersize', 5)
plot(rt_ui(postIdx == 1), wc_ui(postIdx == 1), '.r', 'markersize', 5)

xlabel('radonT-unitID', 'fontsize', 12)
ylabel('wCorr-unitID', 'fontsize', 12)


subplot(2,6,9)
hold on

plot(rt_ts(runIdx == 1), rt_ui(runIdx == 1), '.k', 'markersize', 5)
plot(rt_ts(preIdx == 1), rt_ui(preIdx == 1), '.b', 'markersize', 5)
plot(rt_ts(postIdx == 1), rt_ui(postIdx == 1), '.r', 'markersize', 5)

xlabel('radonT-timeSwap', 'fontsize', 12)
ylabel('radonT-unitID', 'fontsize', 12)


subplot(2,6,10)
hold on

plot(wc_ts(runIdx == 1), wc_ui(runIdx == 1), '.k', 'markersize', 5)
plot(wc_ts(preIdx == 1), wc_ui(preIdx == 1), '.b', 'markersize', 5)
plot(wc_ts(postIdx == 1), wc_ui(postIdx == 1), '.r', 'markersize', 5)

xlabel('wCorr-timeSwap', 'fontsize', 12)
ylabel('wCorr-unitID', 'fontsize', 12)


subplot(2,6,11)
hold on

plot(wc_ts(runIdx == 1), rt_ui(runIdx == 1), '.k', 'markersize', 5)
plot(wc_ts(preIdx == 1), rt_ui(preIdx == 1), '.b', 'markersize', 5)
plot(wc_ts(postIdx == 1), rt_ui(postIdx == 1), '.r', 'markersize', 5)

xlabel('wCorr-timeSwap', 'fontsize', 12)
ylabel('radonT-unitID', 'fontsize', 12)


subplot(2,6,12)
hold on

plot(wc_ts(runIdx == 1), rt_ui(runIdx == 1), '.k', 'markersize', 5)
plot(wc_ts(preIdx == 1), rt_ui(preIdx == 1), '.b', 'markersize', 5)
plot(wc_ts(postIdx == 1), rt_ui(postIdx == 1), '.r', 'markersize', 5)

xlabel('wCorr-unitID', 'fontsize', 12)
ylabel('radonT-timeSwap', 'fontsize', 12)


%% for calculation of correlation of replay scores calcualted using different methods


sessionNames = fieldnames(PBEinfo_pooled);

nSessions  = numel(sessionNames);

rt_ts_pooled = struct('PRE', [], 'RUN', [], 'POST', []);
rt_ui_pooled = struct('PRE', [], 'RUN', [], 'POST', []);

wc_ts_pooled = struct('PRE', [], 'RUN', [], 'POST', []);
wc_ui_pooled = struct('PRE', [], 'RUN', [], 'POST', []);



for ii = 10%1:nSessions
    
   
    PBEinfo = PBEinfo_pooled.(sessionNames{ii});
        
    runIdx = [PBEinfo.RUN]';
    preIdx = [PBEinfo.PRE]';
    postIdx = [PBEinfo.POST]';
    
    rippleOverlap = [PBEinfo.rippleOverlap]';
    
    NREM   = [PBEinfo.NREM]';
    
    
    nFiringUnits = [PBEinfo.nFiringUnits]';
    duration     = [PBEinfo.duration]';
    
    rt_ts = [PBEinfo.rt_ts]';
    rt_ui = [PBEinfo.rt_ui]';
    wc_ts = [PBEinfo.wc_ts]';
    wc_ui = [PBEinfo.wc_ui]';
    
    % % RUN
    
    qualifiedPBEs = find(runIdx == 1 & rippleOverlap == 1 & nFiringUnits >= 5 & duration >= 80 & duration <= 1000);
    rt_ts_pooled.RUN = [rt_ts_pooled.RUN; rt_ts(qualifiedPBEs)];
    rt_ui_pooled.RUN = [rt_ui_pooled.RUN; rt_ui(qualifiedPBEs)];
    wc_ts_pooled.RUN = [wc_ts_pooled.RUN; wc_ts(qualifiedPBEs)];
    wc_ui_pooled.RUN = [wc_ui_pooled.RUN; wc_ui(qualifiedPBEs)];
    
    % % PRE
    qualifiedPBEs   = find(preIdx == 1 & NREM == 1 & nFiringUnits >= 5 & duration >= 80 & duration <= 1000 );
    rt_ts_pooled.PRE = [rt_ts_pooled.PRE; rt_ts(qualifiedPBEs)];
    rt_ui_pooled.PRE = [rt_ui_pooled.PRE; rt_ui(qualifiedPBEs)];
    wc_ts_pooled.PRE = [wc_ts_pooled.PRE; wc_ts(qualifiedPBEs)];
    wc_ui_pooled.PRE = [wc_ui_pooled.PRE; wc_ui(qualifiedPBEs)];
    
    
    % % POST
    qualifiedPBEs   = find(postIdx == 1 &  NREM == 1 & nFiringUnits >= 5 & duration >= 80 & duration <= 1000 );
    rt_ts_pooled.POST = [rt_ts_pooled.POST; rt_ts(qualifiedPBEs)];
    rt_ui_pooled.POST = [rt_ui_pooled.POST; rt_ui(qualifiedPBEs)];
    wc_ts_pooled.POST = [wc_ts_pooled.POST; wc_ts(qualifiedPBEs)];
    wc_ui_pooled.POST = [wc_ui_pooled.POST; wc_ui(qualifiedPBEs)];

end

ts_correlation.PRE = corr(rt_ts_pooled.PRE, wc_ts_pooled.PRE);
ts_correlation.RUN = corr(rt_ts_pooled.RUN, wc_ts_pooled.RUN);
ts_correlation.POST = corr(rt_ts_pooled.POST, wc_ts_pooled.POST);


ui_correlation.PRE = corr(rt_ui_pooled.PRE, wc_ui_pooled.PRE);
ui_correlation.RUN = corr(rt_ui_pooled.RUN, wc_ui_pooled.RUN);
ui_correlation.POST = corr(rt_ui_pooled.POST, wc_ui_pooled.POST);



% % time swap

periodNames = {'PRE'; 'RUN'; 'POST'};
for iperiod = 1:3

    currPeriod = periodNames{iperiod};
    
    x = rt_ts_pooled.(currPeriod);
    y = wc_ts_pooled.(currPeriod);
    
    X = [ones(length(x), 1) x];
    
    b = X\y;
    
    ycal = X*b;
    
    Rsq_ts.(currPeriod) = 1 - sum((y - ycal).^2)/sum((y - mean(y)).^2);
    
end


periodNames = {'PRE'; 'RUN'; 'POST'};
for iperiod = 1:3

    currPeriod = periodNames{iperiod};
    
    x = rt_ui_pooled.(currPeriod);
    y = wc_ui_pooled.(currPeriod);
    
    X = [ones(length(x), 1) x];
    
    b = X\y;
    
    ycal = X*b;
    
    Rsq_ui.(currPeriod) = 1 - sum((y - ycal).^2)/sum((y - mean(y)).^2);
    
end


%% Distribution of replay Scores



figure;
set(gcf, 'Units', 'centimeters', 'position', [0 0 17.6 9])

cl = [204 0 0]/255;
% cl = [0 153 0]/255;


% % PRE
subplot(3,3,[4 7])
plotseqscoredist(rt_ts_pooled.PRE, rt_ui_pooled.PRE, [], [], 0, sprintf('%s\n  (n=%d)', 'PRE', length(rt_ts_pooled.PRE)), cl)
ylabel('ratio of PBEs', 'fontsize', 14)

% % RUN
subplot(3,3,[5 8])
plotseqscoredist(rt_ts_pooled.RUN, rt_ui_pooled.RUN, [], [], 0, sprintf('%s\n  (n=%d)', 'RUN', length(rt_ts_pooled.RUN)), cl)

% % POST
subplot(3,3,[6 9])
plotseqscoredist(rt_ts_pooled.POST, rt_ui_pooled.POST, [], [], 0, sprintf('%s\n  (n=%d)', 'POST', length(rt_ts_pooled.POST)), cl) 


% legend
legendSub = subplot(3,3,3);

hold on

%         p_pts = plot(1, nan, 'color', 'k', 'linestyle', ':', 'linewidth', 3);
p_pts = plot(1, nan, 'color', [.7 .7 .7], 'linestyle', ':', 'linewidth', 2);
p     = plot(1, nan, 'color', cl, 'linestyle', '-', 'linewidth', 3);
p_h   = plot(1, nan, 'color', cl, 'linestyle', '-', 'linewidth', 3);
p_h.Color(4) = 0.5;

h_legend = legend([p, p_h, p_pts],{'Time swap', 'Unit-ID shuffle', 'diagonal'}, 'Location', 'South');


set(h_legend, 'fontsize', 14)
legend boxoff 

set(legendSub, 'Visible', 'off');

% filename = fullfile(subfolder, [BDscoreMethod '-' shuffleMethod]);

% print(gcf, filename, '-dpdf')
        

% close all



