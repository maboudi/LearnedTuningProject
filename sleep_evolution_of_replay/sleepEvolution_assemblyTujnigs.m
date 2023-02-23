clear
clc

parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/assemblyTuning_finalResults';

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

for ii = 7:9%:numel(sessionNames)
    
    
    sessionDir = fullfile(parentDir, sessionNames{ii});


    
    fileName = [sessionNames{ii} '.assemblyTuning_vs_time2.mat'];
    try
        load(fullfile(sessionDir, 'assemblyTunings_entirePost', fileName), 'assemblyTuningPFcorr_time_zscore', 'binCenters')
    catch
        load(fullfile(sessionDir, 'assemblyTunings', fileName), 'assemblyTuningPFcorr_time_zscore', 'binCenters')
    end

    
    load(fullfile(sessionDir, 'spikes', [sessionNames{ii} '.spikes.mat']), 'fileInfo')
    behavior = fileInfo.behavior;
    
    
    
    currSessionName = sessionNames{ii};
    currSessionName(currSessionName == '-') = '_';
   
    
    for iperiod = 1:3
        
        currPeriod = periodNames{iperiod};
        
       [nUnits, nWins] = size(assemblyTuningPFcorr_time_zscore.(currPeriod));
    
%         median_assemblyTuning.(currSessionName).(currPeriod) = nanmedian(assemblyTuningPFcorr_time_zscore.(currPeriod));
        median_assemblyTuning.(currSessionName).(currPeriod) = nanmean(assemblyTuningPFcorr_time_zscore.(currPeriod));

        
        binCenters2.(currSessionName).(currPeriod) = binCenters.(currPeriod);

        lq_assemblyTuning.(currSessionName).(currPeriod) = zeros(nWins, 1);
        hq_assemblyTuning.(currSessionName).(currPeriod) = zeros(nWins, 1);
        assemblytuning.(currSessionName).(currPeriod)    = cell(nWins, 1);
%         sem_assemblyTuning.(currSessionName).(currPeriod) = zeros(nWins, 1);

        for iwin = 1:nWins
            assemblytuning.(currSessionName).(currPeriod){iwin} = assemblyTuningPFcorr_time_zscore.(currPeriod)(:, iwin);
%             lq_assemblyTuning.(currSessionName).(currPeriod)(iwin) = prctile(assemblyTuningPFcorr_time_zscore.(currPeriod)(:, iwin), 25);
%             hq_assemblyTuning.(currSessionName).(currPeriod)(iwin) = prctile(assemblyTuningPFcorr_time_zscore.(currPeriod)(:, iwin), 75);
            
            sem_assemblyTuning = nanstd(assemblyTuningPFcorr_time_zscore.(currPeriod)(:, iwin))/nUnits;

            lq_assemblyTuning.(currSessionName).(currPeriod)(iwin) = median_assemblyTuning.(currSessionName).(currPeriod)(iwin) - sem_assemblyTuning;
            hq_assemblyTuning.(currSessionName).(currPeriod)(iwin) = median_assemblyTuning.(currSessionName).(currPeriod)(iwin) + sem_assemblyTuning;
            
        end
    end
    
    binCenters = binCenters2;
    
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
    

end


%%
sessionNames2 = fieldnames(median_assemblyTuning);
nSessions     = numel(sessionNames2);


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

subwidth.run = (plotwidth - leftedge - rightedge - gapc * (noc - 1))*longest_runDur/(longest_preDur + longest_postDur);
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
    PRE_nBins(ii) = numel(assemblytuning.(sessionNames2{ii}).pre);
end

[nBins_max_PRE, idx] = max(PRE_nBins);
nBins_min_PRE = max(PRE_nBins);



allSessions_PRE = cell(nSessions, nBins_min_PRE);
for ii = 1:nSessions    
    currSessionData = assemblytuning.(sessionNames2{ii}).pre;
    allSessions_PRE(ii, end-numel(currSessionData)+1 : end) = currSessionData;
end


assemblyTuning_median_pooled = zeros(nBins_min_PRE, 1);
assemblyTuning_lq_pooled     = zeros(nBins_min_PRE, 1);
assemblyTuning_hq_pooled     = zeros(nBins_min_PRE, 1);

for ibin = 1:nBins_min_PRE
    
    assemblyTuning_median_pooled(ibin) = nanmedian(cell2mat(allSessions_PRE(:, ibin)));
    assemblyTuning_lq_pooled(ibin)     = prctile(cell2mat(allSessions_PRE(:, ibin)), 25);
    assemblyTuning_hq_pooled(ibin)     = prctile(cell2mat(allSessions_PRE(:, ibin)), 75);
end

binCenters_pooled = (1:nBins_min_PRE) + (nBins_max_PRE - nBins_min_PRE);
binCenters_pooled = binCenters.(sessionNames2{idx}).pre(binCenters_pooled);


ax(1) = axes('position',plotPos.pre,'YGrid','on','FontSize',fontsize,'Box','on','Layer','top');

hold on


for ii = 1:nSessions
    
    
    currbinCenters = binCenters.(sessionNames2{ii}).pre - binCenters.(sessionNames2{ii}).pre(end) ;
    
    displayName = sessionNames2{ii};
    displayName(displayName == '_') = '-';
    
    customPlot(median_assemblyTuning.(sessionNames2{ii}).pre, lq_assemblyTuning.(sessionNames2{ii}).pre, hq_assemblyTuning.(sessionNames2{ii}).pre, currbinCenters/3600, cl(ii, :), 0.4, displayName)

end

currbinCenters = binCenters_pooled - binCenters.(sessionNames2{idx}).pre(end);

% displayName = 'pooled';
% customPlot(assemblyTuning_median_pooled, assemblyTuning_lq_pooled, assemblyTuning_hq_pooled, currbinCenters/3600, 'k', 0.7, displayName)

% yticks(0:10:100)

xlim([currbinCenters(1)/3600 currbinCenters(end)/3600])

xlabel('time (hr)')
ylabel('assembly tuning')
title('PRE')



%% MAZE

RUN_nBins = zeros(nSessions, 1);
for ii = 1:nSessions
    RUN_nBins(ii) = numel(assemblytuning.(sessionNames2{ii}).run);
end

[nBins_max_RUN, idx] = max(RUN_nBins);
nBins_min_RUN = max(RUN_nBins);



allSessions_RUN = cell(nSessions, nBins_min_RUN);

for ii = 1:nSessions
    
    currSessionData = assemblytuning.(sessionNames2{ii}).run;
    allSessions_RUN(ii, 1:numel(currSessionData)) = currSessionData;
    
end


assemblyTuning_median_pooled = zeros(nBins_min_RUN, 1);
assemblyTuning_lq_pooled     = zeros(nBins_min_RUN, 1);
assemblyTuning_hq_pooled     = zeros(nBins_min_RUN, 1);

for ibin = 1:nBins_min_RUN
    
    assemblyTuning_median_pooled(ibin) = nanmedian(cell2mat(allSessions_RUN(:, ibin)));
    assemblyTuning_lq_pooled(ibin)     = prctile(cell2mat(allSessions_RUN(:, ibin)), 25);
    assemblyTuning_hq_pooled(ibin)     = prctile(cell2mat(allSessions_RUN(:, ibin)), 75);
end


binCenters_pooled = (1:nBins_min_RUN) + (nBins_max_RUN - nBins_min_RUN);
binCenters_pooled = binCenters.(sessionNames2{idx}).run(binCenters_pooled);




ax2 = axes('position',plotPos.run,'YGrid','on','FontSize',fontsize,'Box','on','Layer','top');


hold on


for ii = 1:nSessions
    
    
    currbinCenters = binCenters.(sessionNames2{ii}).run - binCenters.(sessionNames2{ii}).run(1);
    
    displayName = sessionNames2{ii};
    displayName(displayName == '_') = '-';
    
    customPlot(median_assemblyTuning.(sessionNames2{ii}).run, lq_assemblyTuning.(sessionNames2{ii}).run, hq_assemblyTuning.(sessionNames2{ii}).run, currbinCenters/3600, cl(ii, :), 0.4, displayName)

end

currbinCenters = binCenters_pooled - binCenters_pooled(1);

% displayName = 'pooled';
% customPlot(assemblyTuning_median_pooled, assemblyTuning_lq_pooled, assemblyTuning_hq_pooled, currbinCenters/3600, 'k', 0.7, displayName)

yl = ylim;
% yticks(linspace(yl(1), yl(2), 4))
% yticks(0:10:100)

xlim([currbinCenters(1)/3600 currbinCenters(end)/3600])

xlabel('time(hr)')
ylabel('assembly tuning')
title('MAZE')

legend('location', 'eastout')



%% POST

POST_nBins = zeros(nSessions, 1);
for ii = 1:nSessions
    POST_nBins(ii) = numel(assemblytuning.(sessionNames2{ii}).post);
end

[nBins_max_POST, idx] = max(POST_nBins);
nBins_min_POST = max(POST_nBins);



allSessions_POST = cell(nSessions, nBins_min_POST);

for ii = 1:nSessions
    
    currSessionData = assemblytuning.(sessionNames2{ii}).post;
    allSessions_POST(ii, 1:numel(currSessionData)) = currSessionData;
    
end


assemblyTuning_median_pooled = zeros(nBins_min_POST, 1);
assemblyTuning_lq_pooled     = zeros(nBins_min_POST, 1);
assemblyTuning_hq_pooled     = zeros(nBins_min_POST, 1);

for ibin = 1:nBins_min_RUN
    
    assemblyTuning_median_pooled(ibin) = nanmedian(cell2mat(allSessions_POST(:, ibin)));
    assemblyTuning_lq_pooled(ibin)     = prctile(cell2mat(allSessions_POST(:, ibin)), 25);
    assemblyTuning_hq_pooled(ibin)     = prctile(cell2mat(allSessions_POST(:, ibin)), 75);
end

binCenters_pooled = 1:nBins_min_POST;
binCenters_pooled = binCenters.(sessionNames2{idx}).post(binCenters_pooled);


ax(2)= axes('position',plotPos.post,'YGrid','on','FontSize',fontsize,'Box','on','Layer','top');


hold on

for ii = 1:nSessions
    
    
    currbinCenters = binCenters.(sessionNames2{ii}).post - binCenters.(sessionNames2{ii}).post(1);
    
    displayName = sessionNames2{ii};
    displayName(displayName == '_') = '-';
    
    customPlot(median_assemblyTuning.(sessionNames2{ii}).post, lq_assemblyTuning.(sessionNames2{ii}).post, hq_assemblyTuning.(sessionNames2{ii}).post, currbinCenters/3600, cl(ii, :), 0.4, displayName)

end

currbinCenters = binCenters_pooled - binCenters_pooled(1);

% displayName = 'pooled';
% customPlot(assemblyTuning_median_pooled, assemblyTuning_lq_pooled, assemblyTuning_hq_pooled, currbinCenters/3600, 'k', 0.7, displayName)

xlim([currbinCenters(1)/3600 currbinCenters(end)/3600])

% yticks(0:)
yticklabels([])

xlabel('time(hr)')
title('POST')

linkaxes(ax, 'y')


%%


function customPlot(var_of_interest_mean, var_of_interest_lq, var_of_interest_hq, binCenters, cl, alph, displayName)


fill([binCenters fliplr(binCenters)],[var_of_interest_lq' fliplr(var_of_interest_hq')],cl, 'linestyle','-', 'EdgeColor', cl, 'FaceAlpha', alph, 'EdgeAlpha', alph, 'displayName', displayName);
% hold on
% plot(binCenters, var_of_interest_mean, 'linestyle', '-', 'color', cl)    

end

