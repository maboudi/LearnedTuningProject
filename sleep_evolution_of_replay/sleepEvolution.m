% load('/home/kouroshmaboudi/Documents/HMM_project/HumanScoreData/PBEinfo_pooled.mat')

clearvars -except PBEinfo_pooled
bathPath = '/home/kouroshmaboudi/Documents/HMM_project/assemblyTuning_Feb10th_2';


sessionNames = fieldnames(PBEinfo_pooled);
periodNames  = {'PRE'; 'RUN'; 'POST'};

binDur = [900 300 900]; % in sec
stepSize = [450 200 450];


for ii = [1 2 3 4 5 6 9 10]
    
PBEinfo = PBEinfo_pooled.(sessionNames{ii});

try
currSessionFolder = fullfile(bathPath, sessionNames{ii});
cd(currSessionFolder)
catch
 temp = sessionNames{ii};
 temp(temp == '_') = '-';
 currSessionFolder = fullfile(bathPath, temp);
 cd(currSessionFolder)
end

load('allVariables.mat', 'behavior')

runIdx  = [PBEinfo.RUN]';
preIdx  = [PBEinfo.PRE]';
postIdx  = [PBEinfo.POST]';


peakT   = [PBEinfo.peakT]';
peakMUA = [PBEinfo.peakMUA]';
radonI  = [PBEinfo.radonIntegral]';
wCorr   = [PBEinfo.weightedCorr]';

wc_ui   = [PBEinfo.wc_ui]';
wc_ts   = [PBEinfo.wc_ts]';
rt_ts   = [PBEinfo.rt_ts]';
rt_ui   = [PBEinfo.rt_ui]';


if ii >= 6
    
    peakT = (peakT- behavior.time(1,1))/1e6;
    behavior.time = (behavior.time - behavior.time(1,1))/1e6;
    
end




%% time course of different metrics

startT.PRE = behavior.time(1,1);
endT.PRE   = behavior.time(2,1);

startT.RUN = behavior.time(2,1);
endT.RUN   = behavior.time(2,2);

startT.POST = behavior.time(2,2);
endT.POST   = behavior.time(3,2);


% % 

var_of_interest = wc_ui;


for iperiod = 1:3
    
    currPeriod = periodNames{iperiod};
    
    includedVaribles = [PBEinfo.(currPeriod)]';

    [varCollect.(sessionNames{ii}).(currPeriod), var_of_interest_mean.(sessionNames{ii}).(currPeriod) , var_of_interest_sem.(sessionNames{ii}).(currPeriod), binCenters.(sessionNames{ii}).(currPeriod)] = timebinning(var_of_interest, peakT, binDur(iperiod), stepSize(iperiod), [startT.(currPeriod) endT.(currPeriod)]);

end


end


sessionNames2 = fieldnames(varCollect);
nSessions = numel(sessionNames2);

%% PRE

PRE_nBins = zeros(nSessions, 1);
for ii = 1:nSessions
    PRE_nBins(ii) = numel(varCollect.(sessionNames2{ii}).PRE);
end

[nBins_max_PRE, idx] = max(PRE_nBins);
nBins_min_PRE = min(PRE_nBins);



allSessions_PRE = cell(nSessions, nBins_min_PRE);

for ii = 1:nSessions
    
    currSessionData = varCollect.(sessionNames2{ii}).PRE;
    allSessions_PRE(ii, :) = currSessionData(end-nBins_min_PRE+1 : end);
    
end


var_of_interest_mean_pooled = zeros(nBins_min_PRE, 1);
var_of_interest_sem_pooled  = zeros(nBins_min_PRE, 1);

for ibin = 1:nBins_min_PRE
    
    var_of_interest_mean_pooled(ibin) = nanmean(cell2mat(allSessions_PRE(:, ibin)));
    var_of_interest_sem_pooled(ibin)  = nanstd(cell2mat(allSessions_PRE(:, ibin)))/numel(cell2mat(allSessions_PRE(:, ibin)));
end

binCenters_pooled = (1:nBins_min_PRE) + (nBins_max_PRE - nBins_min_PRE);
binCenters_pooled = binCenters.(sessionNames{idx}).PRE(binCenters_pooled);

figure;
ax(1) = subplot(1,2,1);
hold on


for ii = 1:nSessions
    
    
    currbinCenters = binCenters_pooled(end)-binCenters.(sessionNames2{ii}).PRE+stepSize(1) - binCenters.(sessionNames{idx}).PRE(end);
    customPlot(var_of_interest_mean.(sessionNames2{ii}).PRE, var_of_interest_sem.(sessionNames2{ii}).PRE, currbinCenters/3600, [204 0 0]/255, 0.1)

end

currbinCenters = binCenters_pooled - binCenters.(sessionNames{idx}).PRE(end);
plot(currbinCenters/3600, var_of_interest_mean_pooled, 'color', [204 0 0]/255, 'linewidth', 3)

% customPlot(var_of_interest_mean_pooled, var_of_interest_sem_pooled, binCenters_pooled/3600, 'r', 1)


%% POST

POST_nBins = zeros(nSessions, 1);
for ii = 1:nSessions
    POST_nBins(ii) = numel(varCollect.(sessionNames2{ii}).POST);
end

[nBins_max_POST, idx] = max(POST_nBins);
nBins_min_POST = min(POST_nBins);



allSessions_POST = cell(nSessions, nBins_min_POST);

for ii = 1:nSessions
    
    currSessionData = varCollect.(sessionNames2{ii}).POST;
    allSessions_POST(ii, :) = currSessionData(1:nBins_min_POST);
    
end


var_of_interest_mean_pooled = zeros(nBins_min_POST, 1);
var_of_interest_sem_pooled  = zeros(nBins_min_POST, 1);

for ibin = 1:nBins_min_POST
    
    var_of_interest_mean_pooled(ibin) = nanmean(cell2mat(allSessions_POST(:, ibin)));
    var_of_interest_sem_pooled(ibin)  = nanstd(cell2mat(allSessions_POST(:, ibin)))/numel(cell2mat(allSessions_POST(:, ibin)));
end

binCenters_pooled = 1:nBins_min_POST;
binCenters_pooled = binCenters.(sessionNames{idx}).POST(binCenters_pooled);


ax(2)= subplot(1,2,2);
hold on

for ii = 1:nSessions
    
    
    currbinCenters = binCenters.(sessionNames2{ii}).POST - binCenters.(sessionNames2{ii}).POST(1) ;
    customPlot(var_of_interest_mean.(sessionNames2{ii}).POST, var_of_interest_sem.(sessionNames2{ii}).POST, currbinCenters/3600, [204 0 0]/255, 0.1)

end

currbinCenters = binCenters_pooled - binCenters_pooled(1);
plot(currbinCenters/3600, var_of_interest_mean_pooled, 'color', [204 0 0]/255, 'linewidth', 3)


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


function customPlot(var_of_interest_mean, var_of_interest_sem, binCenters, cl, alph)


 fill([binCenters fliplr(binCenters)],[var_of_interest_mean'+var_of_interest_sem' fliplr(var_of_interest_mean'-var_of_interest_sem')],cl, 'linestyle','-', 'EdgeColor', cl, 'FaceAlpha', alph, 'EdgeAlpha', alph);


end
