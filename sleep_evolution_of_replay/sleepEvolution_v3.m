% load('/home/kouroshmaboudi/Documents/HMM_project/HumanScoreData/PBEinfo_pooled.mat')


% % in this version indiviual sessions and mean across session are plotted
% versus time during the sleep period


clearvars -except PBEinfo_pooled
bathPath = '/home/kouroshmaboudi/Documents/HMM_project/assemblyTuning_Feb12th';


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

var_of_interest = rt_ts;
cl = 'k';

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

% % added
for ii= 1:nSessions
    
    nNancells_added = nBins_max_PRE - PRE_nBins(ii);
    nancell = cell(nNancells_added, 1);
    nancell(:) = {nan};
    
    varCollect.(sessionNames2{ii}).PRE = [nancell; varCollect.(sessionNames2{ii}).PRE];
    var_of_interest_mean.(sessionNames2{ii}).PRE = [nan(nNancells_added, 1); var_of_interest_mean.(sessionNames2{ii}).PRE];
    var_of_interest_sem.(sessionNames2{ii}).PRE  = [nan(nNancells_added, 1); var_of_interest_sem.(sessionNames2{ii}).PRE];
end
% %


allSessions_PRE = cell(nSessions, nBins_max_PRE);

for ii = 1:nSessions
    
    currSessionData = varCollect.(sessionNames2{ii}).PRE;
    allSessions_PRE(ii, :) = currSessionData;
    
end


var_of_interest_mean_pooled = zeros(nBins_max_PRE, 1);
var_of_interest_sem_pooled  = zeros(nBins_max_PRE, 1);

for ibin = 1:nBins_max_PRE
    
    var_of_interest_mean_pooled(ibin) = nanmean(cell2mat(allSessions_PRE(:, ibin)));
    var_of_interest_sem_pooled(ibin)  = nanstd(cell2mat(allSessions_PRE(:, ibin)))/numel(cell2mat(allSessions_PRE(:, ibin)));
end


figure;
ax(1) = subplot(1,2,1);
hold on


currbinCenters = binCenters.(sessionNames2{idx}).PRE - binCenters.(sessionNames2{idx}).PRE(end);

for ii = 1:nSessions
    
    customPlot(var_of_interest_mean.(sessionNames2{ii}).PRE, var_of_interest_sem.(sessionNames2{ii}).PRE, currbinCenters/3600, cl, 0.1)

end

plot(currbinCenters/3600, var_of_interest_mean_pooled, 'color', cl, 'linewidth', 3)

% customPlot(var_of_interest_mean_pooled, var_of_interest_sem_pooled, binCenters_pooled/3600, 'r', 1)
xlim([currbinCenters(1)/3600 currbinCenters(end)/3600])
xlabel('time(hr) before RUN')

rr = ylim;
text(-3, rr(2) - 0.2*diff(rr), 'PRE', 'fontsize', 14)
ylabel('radon integral', 'fontsize', 12)


%% POST

POST_nBins = zeros(nSessions, 1);
for ii = 1:nSessions
    POST_nBins(ii) = numel(varCollect.(sessionNames2{ii}).POST);
end

[nBins_max_POST, idx] = max(POST_nBins);
nBins_min_POST = min(POST_nBins);

% % added
for ii= 1:nSessions
    
    nNancells_added = nBins_max_POST - POST_nBins(ii);
    nancell = cell(nNancells_added, 1);
    nancell(:) = {nan};
    
    varCollect.(sessionNames2{ii}).POST = [varCollect.(sessionNames2{ii}).POST; nancell];
    var_of_interest_mean.(sessionNames2{ii}).POST = [var_of_interest_mean.(sessionNames2{ii}).POST; nan(nNancells_added, 1)];
    var_of_interest_sem.(sessionNames2{ii}).POST  = [var_of_interest_sem.(sessionNames2{ii}).POST; nan(nNancells_added, 1)];
end
% %


allSessions_POST = cell(nSessions, nBins_max_POST);

for ii = 1:nSessions
    
    currSessionData = varCollect.(sessionNames2{ii}).POST;
    allSessions_POST(ii, :) = currSessionData(1:nBins_max_POST);
    
end


var_of_interest_mean_pooled = zeros(nBins_max_POST, 1);
var_of_interest_sem_pooled  = zeros(nBins_max_POST, 1);

for ibin = 1:nBins_max_POST
    
    var_of_interest_mean_pooled(ibin) = nanmean(cell2mat(allSessions_POST(:, ibin)));
    var_of_interest_sem_pooled(ibin)  = nanstd(cell2mat(allSessions_POST(:, ibin)))/numel(cell2mat(allSessions_POST(:, ibin)));
end


ax(2)= subplot(1,2,2);
hold on


currbinCenters = binCenters.(sessionNames2{idx}).POST - binCenters.(sessionNames2{idx}).POST(1);

for ii = 1:nSessions
    customPlot(var_of_interest_mean.(sessionNames2{ii}).POST, var_of_interest_sem.(sessionNames2{ii}).POST, currbinCenters/3600, cl, 0.1)
end

plot(currbinCenters/3600, var_of_interest_mean_pooled, 'color', cl, 'linewidth', 3)

xlim([currbinCenters(1)/3600 currbinCenters(end)/3600])

% 
% figure; 
% customPlot(var_of_interest_mean_pooled, var_of_interest_sem_pooled, currbinCenters/3600, [0 153 0]/255, 0.1)



linkaxes(ax, 'y')
pre_panel = get(ax(1), 'position');
post_panel = get(ax(2), 'position');

xlabel('time(hr) after RUN')

set(ax(2), 'position', [post_panel(1) post_panel(2) post_panel(3)*nBins_max_POST/nBins_max_PRE post_panel(4)])

rr = ylim;
text(1.5, rr(2) - 0.2*diff(rr), 'POST', 'fontsize', 14)

set(gcf, 'position', [100 100 613 350])

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
