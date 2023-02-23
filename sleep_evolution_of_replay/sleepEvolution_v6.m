clear
clc

parentDir = '/data/concat_GreatLakes_datasets_temp_beforeAssemblyTunings';

files = dir(parentDir);
isdir = [files.isdir];
subfolders = files(isdir);

sessionNames = {subfolders(3:end).name};

periodNames  = {'pre'; 'run'; 'post'};

binDur = [900 900 900]; % in sec
stepSize = [450 450 450];
    
longest_preDur = -1;
longest_runDur = -1;
longest_postDur = -1;

for ii = 9%1:numel(sessionNames)
    
    
    sessionDir = fullfile(parentDir, sessionNames{ii});
    
    load(fullfile(sessionDir, 'BayesianDecoding', [sessionNames{ii} '.PBEInfo_replayScores.mat']), 'PBEInfo_replayScores')
    PBEInfo = PBEInfo_replayScores;
    
    PBEInfo = PBEInfo(ismember({PBEInfo.brainState}, {'NREM'; 'QW'}));
    
    clear PBEInfo_replayScores
    
    load(fullfile(sessionDir, 'spikes', [sessionNames{ii} '.spikes.mat']), 'fileInfo')
    behavior = fileInfo.behavior;
    
    peakT   = [PBEInfo.peakT]';
    peakMUA = [PBEInfo.peakMUA]';
    radonI  = [PBEInfo.radonIntegral]';
    wCorr   = [PBEInfo.weightedCorr]';

    wc_ui   = [PBEInfo.wc_ui]';
    wc_ts   = [PBEInfo.wc_ts]';
    rt_ts   = [PBEInfo.rt_ts]';
    rt_ui   = [PBEInfo.rt_ui]';

    
    
    %% time course of different metrics

    startT.pre = behavior.time(1,1);
    endT.pre   = behavior.time(2,1);
    dur.pre    = endT.pre - startT.pre;
    if dur.pre > longest_preDur
        longest_preDur = dur.pre;
    end
    

    startT.run = behavior.time(2,1);
    endT.run   = behavior.time(2,2);
    dur.run    = endT.run - startT.run;
    if dur.run > longest_runDur
        longest_runDur = dur.run;
    end
    
    
    startT.post = behavior.time(2,2);
    endT.post   = behavior.time(3,2);
    dur.post    = endT.post - startT.post;
    
    if dur.post > longest_postDur
        longest_postDur = dur.post;
    end
    

    
    var_of_interest = wc_ts;
    

    for iperiod = 1:3

        currPeriod = periodNames{iperiod};
        
        currSessionName = sessionNames{ii};
        currSessionName(currSessionName == '-') = '_';

        [varCollect.(currSessionName).(currPeriod), var_of_interest_mean.(currSessionName).(currPeriod) , var_of_interest_sem.(currSessionName).(currPeriod), binCenters.(currSessionName).(currPeriod)] = timebinning(var_of_interest, peakT, binDur(iperiod), stepSize(iperiod), [startT.(currPeriod) endT.(currPeriod)]);

    end

end


%%
sessionNames2 = fieldnames(varCollect);
nSessions = numel(sessionNames2);


colors= colormap('jet');
% colors = colors(randperm(length(colors)), :);

cl = zeros(nSessions, 3);
for ii = 1:nSessions
    cl(ii, :) = colors(floor(64/nSessions * ii), :); 
end


%%

% figure configuration

plotheight = 15; %% in cm
plotwidth = 20;
fontsize = 12;

nor = 2;
noc = 2; %% number of events per column

leftedge = 2;    rightedge = 2;    topedge = 2;    bottomedge = 2;

gapc = 1;
gapr = 2;


proportions = [1 longest_runDur/longest_preDur longest_postDur/longest_preDur];

subheight = (plotheight - topedge - bottomedge - gapr * (nor - 1))/2;


% PRE

xfirst = leftedge;
yfirst = bottomedge;

subwidth.pre  = (plotwidth - leftedge - rightedge - gapc * (noc - 1))*longest_preDur/(longest_preDur + longest_postDur);
plotPos.pre = [xfirst/plotwidth yfirst/plotheight subwidth.pre/plotwidth subheight/plotheight];

% RUN

xfirst = leftedge;
yfirst = bottomedge + subheight + gapr;

subwidth.run = (plotwidth - leftedge - rightedge - gapc * (noc - 1))*longest_runDur/(longest_preDur+longest_postDur);
plotPos.run  = [xfirst/plotwidth yfirst/plotheight subwidth.run/plotwidth subheight/plotheight];


% POST

xfirst = leftedge + subwidth.pre + gapc;
yfirst = bottomedge;

subwidth.post = (plotwidth - leftedge - rightedge - gapc * (noc - 1))*longest_postDur/(longest_preDur + longest_postDur);
plotPos.post  = [xfirst/plotwidth yfirst/plotheight subwidth.post/plotwidth subheight/plotheight];



f=figure;
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);




%% PRE

PRE_nBins = zeros(nSessions, 1);
for ii = 1:nSessions
    PRE_nBins(ii) = numel(varCollect.(sessionNames2{ii}).pre);
end

[nBins_max_PRE, idx] = max(PRE_nBins);
nBins_min_PRE = max(PRE_nBins);



allSessions_PRE = cell(nSessions, nBins_min_PRE);
for ii = 1:nSessions    
    currSessionData = varCollect.(sessionNames2{ii}).pre;
    allSessions_PRE(ii, end-numel(currSessionData)+1 : end) = currSessionData;
end


var_of_interest_mean_pooled = zeros(nBins_min_PRE, 1);
var_of_interest_sem_pooled  = zeros(nBins_min_PRE, 1);

for ibin = 1:nBins_min_PRE
    
    var_of_interest_mean_pooled(ibin) = nanmean(cell2mat(allSessions_PRE(:, ibin)));
    var_of_interest_sem_pooled(ibin)  = nanstd(cell2mat(allSessions_PRE(:, ibin)))/numel(cell2mat(allSessions_PRE(:, ibin)));
end

binCenters_pooled = (1:nBins_min_PRE) + (nBins_max_PRE - nBins_min_PRE);
binCenters_pooled = binCenters.(sessionNames2{idx}).pre(binCenters_pooled);


ax(1) = axes('position',plotPos.pre,'YGrid','on','FontSize',fontsize,'Box','on','Layer','top');

hold on


for ii = 1:nSessions
    
    
    currbinCenters = binCenters.(sessionNames2{ii}).pre - binCenters.(sessionNames2{ii}).pre(end) ;
    
    displayName = sessionNames2{ii};
    displayName(displayName == '_') = '-';
    
    customPlot(var_of_interest_mean.(sessionNames2{ii}).pre, var_of_interest_sem.(sessionNames2{ii}).pre, currbinCenters/3600, cl(ii, :), 0.4, displayName)

end

currbinCenters = binCenters_pooled - binCenters.(sessionNames2{idx}).pre(end);

displayName = 'pooled';
customPlot(var_of_interest_mean_pooled, var_of_interest_sem_pooled, currbinCenters/3600, 'k', 0.7, displayName)

yticks(0:10:100)

xlim([currbinCenters(1)/3600 currbinCenters(end)/3600])

xlabel('time (hr)')
ylabel('replay score %')
title('PRE')



%% MAZE

RUN_nBins = zeros(nSessions, 1);
for ii = 1:nSessions
    RUN_nBins(ii) = numel(varCollect.(sessionNames2{ii}).run);
end

[nBins_max_RUN, idx] = max(RUN_nBins);
nBins_min_RUN = max(RUN_nBins);



allSessions_RUN = cell(nSessions, nBins_min_RUN);

for ii = 1:nSessions
    
    currSessionData = varCollect.(sessionNames2{ii}).run;
    allSessions_RUN(ii, 1:numel(currSessionData)) = currSessionData;
    
end



var_of_interest_mean_pooled = zeros(nBins_min_RUN, 1);
var_of_interest_sem_pooled  = zeros(nBins_min_RUN, 1);

for ibin = 1:nBins_min_RUN
    
    var_of_interest_mean_pooled(ibin) = nanmean(cell2mat(allSessions_RUN(:, ibin)));
    var_of_interest_sem_pooled(ibin)  = nanstd(cell2mat(allSessions_RUN(:, ibin)))/numel(cell2mat(allSessions_RUN(:, ibin)));
end

binCenters_pooled = (1:nBins_min_RUN) + (nBins_max_RUN - nBins_min_RUN);
binCenters_pooled = binCenters.(sessionNames2{idx}).run(binCenters_pooled);




ax2 = axes('position',plotPos.run,'YGrid','on','FontSize',fontsize,'Box','on','Layer','top');


hold on


for ii = 1:nSessions
    
    
    currbinCenters = binCenters.(sessionNames2{ii}).run - binCenters.(sessionNames2{ii}).run(1);
    
    displayName = sessionNames2{ii};
    displayName(displayName == '_') = '-';
    
    customPlot(var_of_interest_mean.(sessionNames2{ii}).run, var_of_interest_sem.(sessionNames2{ii}).run, currbinCenters/3600, cl(ii, :), 0.4, displayName)

end

currbinCenters = binCenters_pooled - binCenters_pooled(1);

displayName = 'pooled';
customPlot(var_of_interest_mean_pooled, var_of_interest_sem_pooled, currbinCenters/3600, 'k', 0.7, displayName)

yl = ylim;
yticks(floor(linspace(yl(1), yl(2), 4)))
% yticks(0:10:100)

xlim([currbinCenters(1)/3600 currbinCenters(end)/3600])

xlabel('time(hr)')
ylabel('replay score %')
title('MAZE')

legend('location', 'eastout')


%% POST

POST_nBins = zeros(nSessions, 1);
for ii = 1:nSessions
    POST_nBins(ii) = numel(varCollect.(sessionNames2{ii}).post);
end

[nBins_max_POST, idx] = max(POST_nBins);
nBins_min_POST = max(POST_nBins);



allSessions_POST = cell(nSessions, nBins_min_POST);

for ii = 1:nSessions
    
    currSessionData = varCollect.(sessionNames2{ii}).post;
    allSessions_POST(ii, 1:numel(currSessionData)) = currSessionData;
    
end


var_of_interest_mean_pooled = zeros(nBins_min_POST, 1);
var_of_interest_sem_pooled  = zeros(nBins_min_POST, 1);

for ibin = 1:nBins_min_POST
    
    var_of_interest_mean_pooled(ibin) = nanmean(cell2mat(allSessions_POST(:, ibin)));
    var_of_interest_sem_pooled(ibin)  = nanstd(cell2mat(allSessions_POST(:, ibin)))/numel(cell2mat(allSessions_POST(:, ibin)));
end

binCenters_pooled = 1:nBins_min_POST;
binCenters_pooled = binCenters.(sessionNames2{idx}).post(binCenters_pooled);


ax(2)= axes('position',plotPos.post,'YGrid','on','FontSize',fontsize,'Box','on','Layer','top');


hold on

for ii = 1:nSessions
    
    
    currbinCenters = binCenters.(sessionNames2{ii}).post - binCenters.(sessionNames2{ii}).post(1);
    
    displayName = sessionNames2{ii};
    displayName(displayName == '_') = '-';
    
    customPlot(var_of_interest_mean.(sessionNames2{ii}).post, var_of_interest_sem.(sessionNames2{ii}).post, currbinCenters/3600, cl(ii, :), 0.4, displayName)

end

currbinCenters = binCenters_pooled - binCenters_pooled(1);

displayName = 'pooled';
customPlot(var_of_interest_mean_pooled, var_of_interest_sem_pooled, currbinCenters/3600, 'k', 0.7, displayName)

xlim([currbinCenters(1)/3600 currbinCenters(end)/3600])

yticks(0:10:100)
yticklabels([])

xlabel('time(hr)')
title('POST')

linkaxes(ax, 'y')


%% functions


function [varCollect, var_of_interest_mean, var_of_interest_sem, binCenters] = timebinning(var_of_interest, peakT, binDur, stepDur, period)

totalT = period(2) - period(1);
nBins  = floor((totalT - binDur)/stepDur) + 1;

startT     = period(1) + (0:nBins-1)*stepDur;
endT       = startT + binDur;
binCenters = startT + binDur/2;

varCollect = cell(nBins, 1);
var_of_interest_mean = zeros(nBins, 1);
var_of_interest_sem = zeros(nBins, 1);

for ii = 1:nBins
    
   idx = peakT > startT(ii) & peakT <= endT(ii);
   
   varCollect{ii} = var_of_interest(idx);
   var_of_interest_mean(ii) = nanmean(varCollect{ii});
   var_of_interest_sem(ii)  = nanstd(varCollect{ii})/numel(find(idx));
   
    
end

end


function customPlot(var_of_interest_mean, var_of_interest_sem, binCenters, cl, alph, displayName)


 fill([binCenters fliplr(binCenters)],[var_of_interest_mean'+var_of_interest_sem' fliplr(var_of_interest_mean'-var_of_interest_sem')],cl, 'linestyle','-', 'EdgeColor', cl, 'FaceAlpha', alph, 'EdgeAlpha', alph, 'displayName', displayName);


end
