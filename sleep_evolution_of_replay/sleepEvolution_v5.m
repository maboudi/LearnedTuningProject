% load('PBEinfo_pooled.mat')
 
% % plot quartiles of replay scores versus time

clearvars -except PBEinfo_pooled
bathPath = '~/Documents/HMM_project/GreatLakes_firstRun_Nov2020/Grosmark_originalClusters';


sessionNames = fieldnames(PBEinfo_pooled);
periodNames  = {'PRE'; 'RUN'; 'POST'};


figure;
for ii = 1:14

ll = PBEinfo_pooled.(sessionNames{ii});

try
    currSessionFolder = fullfile(bathPath, sessionNames{ii}, 'BayesianDecoding');
    cd(currSessionFolder)
catch
     temp = sessionNames{ii};
     temp(temp == '_') = '-';
     currSessionFolder = fullfile(bathPath, temp, 'BayesianDecoding');
     cd(currSessionFolder)
end

load('allVariables.mat', 'behavior')

runIdx  = [ll.RUN]';
preIdx  = [ll.PRE]';
postIdx = [ll.POST]';


NREMidx = [ll.NREM]'; % filtering used during PRE and POST behavioral epochs
QWidx   = [ll.QW]';
rippleOverlap = [ll.rippleOverlap]'; % filtering PBEs during RUN


peakT   = [ll.peakT]';
peakMUA = [ll.peakMUA]';
radonI  = [ll.radonIntegral]';
wCorr   = [ll.weightedCorr]';

wc_ui   = [ll.wc_ui]';
wc_ts   = [ll.wc_ts]';
rt_ts   = [ll.rt_ts]';
rt_ui   = [ll.rt_ui]';



% adjusting the timestamps in Hiro's dataset
if ii >= 9    
    peakT = (peakT- behavior.time(1,1))/1e6;
    behavior.time = (behavior.time - behavior.time(1,1))/1e6;
end

% 
% peakT = peakT(NREMidx == 1 | QWidx == 1);
% radonI = radonI(NREMidx == 1 | QWidx == 1);
% wCorr = wCorr(NREMidx == 1 | QWidx == 1);
% wc_ui = wc_ui(NREMidx == 1 | QWidx == 1);
% wc_ts = wc_ts(NREMidx == 1 | QWidx == 1);
% rt_ts = rt_ts(NREMidx == 1 | QWidx == 1);
% rt_ui = rt_ui(NREMidx == 1 | QWidx == 1);


binDur = [1800 900 1800]; % in sec
stepSize = [450 300 450];

if ii>= 9 % on eof hiro's sessions
    binDur(2) = 1800;
    stepSize(2) = 450;
end




%% time course of 

startT.PRE = behavior.time(1,1);
endT.PRE   = behavior.time(2,1);

startT.RUN = behavior.time(2,1);
endT.RUN   = behavior.time(2,2);

startT.POST = behavior.time(2,2);
endT.POST   = behavior.time(3,2);



subplot(3,5,ii)
hold on

minY = 0;
maxY = 1;

% patch([startT.RUN endT.RUN endT.RUN startT.RUN]/3600, [minY  minY maxY maxY], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')

alpha_level = [1 0.8 0.6 0.4];

cl = colormap('jet');

cl = cl(15*(1:4), :);

qtlNames = {'qtl 1'; 'qtl 2'; 'qtl 3'; 'qtl 4'};

for iperiod = 1:3
    
    currPeriod = periodNames{iperiod};
    
    
    nBins = floor((diff([startT.(currPeriod) endT.(currPeriod)]) - binDur(iperiod)) / stepSize(iperiod)) + 1;
    
    counts.(currPeriod) = zeros(4,nBins);
    binCenters.(currPeriod) = zeros(4,nBins);
    
    if iperiod == 2
    
        idx = rippleOverlap == 1;
    else
        idx = NREMidx == 1 | QWidx == 1;
    end
    
    var_of_interest = rt_ui(idx);
    peakTimes = peakT(idx);
    
    
    for iq = 1:4
        
        currQrtlidx = var_of_interest >= (iq-1)*25 & var_of_interest <= iq*25; %   need to revise this later
        
        currQrtlValues = var_of_interest(currQrtlidx);
        currPeakTimes = peakTimes(currQrtlidx);
        
        [~, var_of_interest_mean, var_of_interest_sem, counts.(currPeriod)(iq, :),  binCenters.(currPeriod)(iq, :)] = timebinning(currQrtlValues, currPeakTimes, binDur(iperiod), stepSize(iperiod), [startT.(currPeriod) endT.(currPeriod)]);
%         customPlot(var_of_interest_mean, var_of_interest_sem, binCenters/3600, [204 0 0]/255, alpha_level(iq))

    end

end


for iperiod = 1:3
    
    
%     if ii < 6 && iperiod == 2
%         continue
%     end
    currPeriod = periodNames{iperiod};
    
    currCounts = counts.(currPeriod);
    currRatios = currCounts./sum(currCounts, 1);
    
%     bar(binCenters.(currPeriod)(1,:)/3600, currRatios', 'stacked')
    
    
    for iq = 1:4

        p = plot(binCenters.(currPeriod)(iq, :)/3600, currCounts(iq, :), 'color', cl(iq, :), 'linewidth', 2, 'DisplayName', qtlNames{iq});
%         p.Color(4) = alpha_level(iq);
                
    end
    
end



% ylim([minY maxY])
xlim([startT.PRE endT.POST]/3600)

title(sessionNames{ii}, 'fontsize', 12, 'Interpreter', 'none')
xlabel('time (hour)', 'fontsize', 12)
ylabel('ratio of PBEs w/in each quartile of RT-UI')

% plot(binCenters, cellfun(@numel, varCollect), '-', 'linewidth', 3, 'color', 'k')


end

%% functions


function [varCollect, var_of_interest_mean, var_of_interest_sem, counts, binCenters] = timebinning(var_of_interest, peakT, binDur, stepDur, period)

totalT = period(2) - period(1);
nBins  = floor((totalT - binDur)/stepDur) + 1;

startT     = period(1) + (0:nBins-1)*stepDur;
endT       = startT + binDur;
binCenters = startT + binDur/2;

varCollect = cell(nBins, 1);
var_of_interest_mean = zeros(nBins, 1);
var_of_interest_sem = zeros(nBins, 1);
counts = zeros(nBins, 1);

for ii = 1:nBins
    
   idx = peakT > startT(ii) & peakT <= endT(ii);
   
   varCollect{ii} = var_of_interest(idx);
   var_of_interest_mean(ii) = nanmean(varCollect{ii});
   var_of_interest_sem(ii)  = nanstd(varCollect{ii})/numel(find(idx));
   counts(ii) = numel(varCollect{ii});
    
end

end


function customPlot(var_of_interest_mean, var_of_interest_sem, binCenters, cl, alph)


tt = fill([binCenters fliplr(binCenters)],[var_of_interest_mean'+var_of_interest_sem' fliplr(var_of_interest_mean'-var_of_interest_sem')],cl, 'linestyle','-', 'EdgeColor', cl, 'FaceAlpha', alph);


end
