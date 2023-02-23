% load('PBEinfo_pooled.mat')


% % two variables are plotted on two y-axes  


clearvars -except PBEinfo_pooled
bathPath = '/home/kouroshmaboudi/Documents/HMM_project/assemblyTuning_Feb12th';


sessionNames = fieldnames(PBEinfo_pooled);
periodNames  = {'PRE'; 'RUN'; 'POST'};

ii = 1;

ll = PBEinfo_pooled.(sessionNames{ii});

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

runIdx  = [ll.RUN]';
preIdx  = [ll.PRE]';
postIdx = [ll.POST]';

NREMidx = [ll.NREM]';
rippleOverlap = [ll.rippleOverlap]';


peakT   = [ll.peakT]';
peakMUA = [ll.peakMUA]';
radonI  = [ll.radonIntegral]';
wCorr   = [ll.weightedCorr]';

wc_ui   = [ll.wc_ui]';
wc_ts   = [ll.wc_ts]';
rt_ts   = [ll.rt_ts]';
rt_ui   = [ll.rt_ui]';


if ii >= 6
    
    peakT = (peakT- behavior.time(1,1))/1e6;

    behavior.time = (behavior.time - behavior.time(1,1))/1e6;
    
end

% peakT = peakT(NREMidx == 1);
% radonI = radonI(NREMidx == 1);
% wCorr = wCorr(NREMidx == 1);
% wc_ui = wc_ui(NREMidx == 1);
% wc_ts = wc_ts(NREMidx == 1);
% rt_ts = rt_ts(NREMidx == 1);
% rt_ui = rt_ui(NREMidx == 1);

binDur = [1800 400 1800]; % in sec
stepSize = [450 100 450];


%% time course of 

pre_startT = behavior.time(1,1);
pre_endT   = behavior.time(1,2);

run_startT = behavior.time(2,1);
run_endT   = behavior.time(2,2);

post_startT = behavior.time(3,1);
post_endT   = behavior.time(3,2);


var_of_interest = rt_ts;


figure; 
hold on

minY = 0;
maxY = 1;

% patch([pre_startT pre_endT pre_endT pre_startT]/3600, [minY  minY maxY maxY], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
patch([run_startT run_endT run_endT run_startT]/3600, [minY  minY maxY maxY], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
% patch([post_startT post_endT post_endT post_startT]/3600, [minY  minY maxY maxY], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none')

yyaxis left


% PRE
[varCollect, var_of_interest_mean, var_of_interest_sem, binCenters] = timebinning(var_of_interest, peakT, binDur(1), stepSize(1), [pre_startT pre_endT]);
customPlot(var_of_interest_mean, var_of_interest_sem, binCenters/3600, [204 0 0]/255, 0.7)
minY = min(var_of_interest_mean);
maxY = max(var_of_interest_mean);


% RUN
[varCollect, var_of_interest_mean, var_of_interest_sem, binCenters] = timebinning(var_of_interest, peakT, binDur(2), stepSize(2), [run_startT run_endT]);
customPlot(var_of_interest_mean, var_of_interest_sem, binCenters/3600, [204 0 0]/255, 0.7)
minY = min(minY, min(var_of_interest_mean));
maxY = max(maxY, max(var_of_interest_mean));

% POST
[varCollect, var_of_interest_mean, var_of_interest_sem, binCenters] = timebinning(var_of_interest, peakT, binDur(3), stepSize(3), [post_startT post_endT]);
customPlot(var_of_interest_mean, var_of_interest_sem, binCenters/3600, [204 0 0]/255, 0.7)
minY = min(minY, min(var_of_interest_mean));
maxY = max(maxY, max(var_of_interest_mean));

minY1 = minY- 0.1*(maxY-minY);
maxY1 = maxY+ 0.1*(maxY-minY);

ylim([minY maxY])
xlim([pre_startT post_endT]/3600)


yyaxis right
var_of_interest = wc_ts;

% PRE
[varCollect, var_of_interest_mean, var_of_interest_sem, binCenters] = timebinning(var_of_interest, peakT, binDur(1), stepSize(1), [pre_startT pre_endT]);
customPlot(var_of_interest_mean, var_of_interest_sem, binCenters/3600, [0 153 0]/255, 0.7)
minY = min(var_of_interest_mean);
maxY = max(var_of_interest_mean);

% RUN
[varCollect, var_of_interest_mean, var_of_interest_sem, binCenters] = timebinning(var_of_interest, peakT, binDur(2), stepSize(2), [run_startT run_endT]);
customPlot(var_of_interest_mean, var_of_interest_sem, binCenters/3600, [0 153 0]/255, 0.7)
minY = min(minY, min(var_of_interest_mean));
maxY = max(maxY, max(var_of_interest_mean));
% POST
[varCollect, var_of_interest_mean, var_of_interest_sem, binCenters] = timebinning(var_of_interest, peakT, binDur(3), stepSize(3), [post_startT post_endT]);
customPlot(var_of_interest_mean, var_of_interest_sem, binCenters/3600, [0 153 0]/255, 0.7)
minY = min(minY, min(var_of_interest_mean));
maxY = max(maxY, max(var_of_interest_mean));

minY2 = minY- 0.1*(maxY-minY);
maxY2 = maxY+ 0.1*(maxY-minY);

ylim([minY maxY])
xlim([pre_startT post_endT]/3600)

rr = gca;
rr.YAxis(1).Color = [204 0 0]/255;
rr.YAxis(2).Color = [0 153 0]/255;


xlabel('time(hr)', 'fontsize', 12)
rr.YAxis(1).Label.String   = 'radon integral';
rr.YAxis(1).Label.FontSize = 12;
rr.YAxis(1).Label.Color = [204 0 0]/255;
rr.YAxis(1).Limits = [minY1 maxY1];

rr.YAxis(2).Label.String   = 'weighted correlation';
rr.YAxis(2).Label.FontSize = 12;
rr.YAxis(2).Label.Color = [0 153 0]/255;
rr.YAxis(2).Limits = [minY2 maxY2];

% plot(binCenters, cellfun(@numel, varCollect), '-', 'linewidth', 3, 'color', 'k')


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


tt = fill([binCenters fliplr(binCenters)],[var_of_interest_mean'+var_of_interest_sem' fliplr(var_of_interest_mean'-var_of_interest_sem')],cl, 'linestyle','-', 'EdgeColor', cl, 'FaceAlpha', alph);


end
