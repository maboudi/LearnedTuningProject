function assemblyTuning_vs_time_woAVG_cleanerVersion(sessionNumber)


% we might need to consider only stable units when calculating the assembly
% tuning for a specific unit

% sz = getenv('SLURM_CPUS_PER_TASK');
% p  = parpool(str2double(sz));


% addpath(genpath('/nfs/turbo/umms-kdiba/NCMLproject/ReplayPreplayAnalyses'))


% parentDir = '/nfs/turbo/umms-kdiba/Kourosh/NCMLproject/concat_GreatLakes_datasets_temp';
parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/assemblyTuning_finalResults/';


rr = dir(parentDir);
sessionName = rr(sessionNumber+2).name;


basePath    = fullfile(parentDir, sessionName);


load(fullfile(basePath, 'BayesianDecoding', [sessionName '.PBEInfo_replayScores.mat']), 'PBEInfo_replayScores')
load(fullfile(basePath, 'spikes', [sessionName '.spikes.mat']), 'spikes_pyr', 'fileInfo')
load(fullfile(basePath, 'spikes', [sessionName '.clusterQuality.mat']), 'clusterQuality')

% load assemblyTunings

s = load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTuning_vs_time4.mat']), 'binCenters', 'assemblyTuningPFcorr_time', 'asTuning_prefPos_time', 'assemblyTunings_time');
AT_prefLocation_avg   = s.asTuning_prefPos_time;
% AT_placeFieldCorr_avg = s.assemblyTuningPFcorr_time;
binCenters_avg        = s.binCenters;
assemblyTunings_time  = s.assemblyTunings_time;


% behavioral epochs
behavior = fileInfo.behavior;

startT.pre  = behavior.time(1,1); endT.pre  = behavior.time(2,1); 
startT.run  = behavior.time(2,1); endT.run  = behavior.time(2,2); 
startT.post = behavior.time(2,2); endT.post = behavior.time(3,2); 


spikes = spikes_pyr;
nPosBins = numel(spikes(1).spatialTuning_smoothed.uni);


nUnits  = numel(spikes); 



% non-directional spatial tuning
spatialTunings_merge = zeros(nUnits, nPosBins);
peakPFfiring = zeros(nUnits, 1);
for iUnit = 1:nUnits
    spatialTunings_merge(iUnit, :) = spikes(iUnit).spatialTuning_smoothed.uni;
    peakPFfiring(iUnit) = spikes(iUnit).peakFR.uni;
    
end


% we are going to analyze assembly tuning separately for each epoch, so we
% need to reorganize the PBEs


% PBEInfo = PBEInfo(acceptedIdx); % only PBEs qulaified in terms of duration, number of participant units, the brain satte during which they occurred
PBEInfo = PBEInfo_replayScores;
PBEInfo = PBEInfo(ismember({PBEInfo.brainState}, {'QW'; 'NREM'}));



%% assembly Tunings considering all of the PBEs
% Spatial information and correlation with place fields of assesmbly tunings
% and corresponding z-scores by comparing the values againts distibutions calculated based on unit ID shuffle surrogates


% ifshuffle = 0;
% [AT_prefLocation.data, AT_placeFieldCorr.data, AT_placeFieldKL.data, unit_firings, concatPBEbinCenters] = calculateAssemblyTuning_woAVG(PBEInfo, spikes, clusterQuality, ifshuffle); % 



% ifshuffle = 1;
% [AT_prefLocation.ui, AT_placeFieldCorr.ui, AT_placeFieldKL.ui] = calculateAssemblyTuning_woAVG(PBEInfo, spikes, clusterQuality, ifshuffle); % 


%% 


ntPBEs = numel(PBEInfo);
binDur = 0.02;
unit_firings = cell(ntPBEs, 1);
PBEbinCenters = cell(ntPBEs, 1);

for ipbe = 1:ntPBEs    
    unit_firings{ipbe} = PBEInfo(ipbe).fr_20msbin; 
    PBEbinCenters{ipbe} = PBEInfo(ipbe).startT + (1:size(unit_firings{ipbe}, 2))*binDur - binDur/2;  
end

unit_firings = cell2mat(unit_firings');
concatPBEbinCenters = cell2mat(PBEbinCenters');




%% calculate the KL distance for the averaged tunings

epochNames = {'pre'; 'run'; 'post'};
dataTypes  = {'data'; 'ui'};

spatialTuning = spatialTunings_merge';
spatialTuning = spatialTuning ./ repmat(sum(spatialTuning, 1), [size(spatialTuning, 1) 1]);
spatialTuning = spatialTuning + eps;


for iEpoch = 1:3
    
    currEpoch = epochNames{iEpoch};
     
    nTimeBins = size(assemblyTunings_time.(currEpoch).data, 3);
    
    for idata = 1:2
        curDataType = dataTypes{idata};
        asTuning_PF_KLdistance.(currEpoch).(curDataType) = zeros(nUnits, nTimeBins);
        asTuning_PF_KSstat.(currEpoch).(curDataType)     = zeros(nUnits, nTimeBins);
        
        asTuning_PFcorrelation.(currEpoch).(curDataType) = zeros(nUnits, nTimeBins);

        for iUnit = 1:nUnits

            learnedTuning = squeeze(assemblyTunings_time.(currEpoch).(curDataType)(iUnit, :, :, 1));
            learnedTuning = learnedTuning ./ repmat(sum(learnedTuning, 1), [size(learnedTuning, 1) 1]);
            learnedTuning = learnedTuning + eps;

            asTuning_PF_KLdistance.(currEpoch).(curDataType)(iUnit, :)  = sum(learnedTuning .* (log(learnedTuning) - log(repmat(spatialTuning(:, iUnit), [1 size(learnedTuning, 2)]))), 1);
            
            try
                for iBin = 1:size(learnedTuning, 2)
                    
                    x1 = learnedTuning(:, iBin);
                    x2 = spatialTuning(:, iUnit);
                                       
                    asTuning_PF_KSstat.(currEpoch).(curDataType)(iUnit, iBin) = max(abs(cumsum(x1) - cumsum(x2)));
                    

%                     % Compute the asymptotic P-value approximation and accept or
%                     % reject the null hypothesis on the basis of the P-value.
%                     %
% 
%                     n1     =  length(x1);
%                     n2     =  length(x2);
%                     n      =  n1 * n2 /(n1 + n2);
%                     lambda =  max((sqrt(n) + 0.12 + 0.11/sqrt(n)) * KSstatistic , 0);
% 
%                     if tail ~= 0        % 1-sided test.
% 
%                        pValue  =  exp(-2 * lambda * lambda);
% 
%                     else                % 2-sided test (default).
%                     %
%                     %  Use the asymptotic Q-function to approximate the 2-sided P-value.
%                     %
%                        j       =  (1:101)';
%                        pValue  =  2 * sum((-1).^(j-1).*exp(-2*lambda*lambda*j.^2));
%                        pValue  =  min(max(pValue, 0), 1);
% 
%                     end
    

                end
            catch
            end
            
            asTuning_PFcorrelation.(currEpoch).(curDataType)(iUnit, :) = corr(learnedTuning, spatialTuning(:, iUnit));
        end
    end
end


%%

binDur  = 900; % 15 minutes
stepDur = 300; % 5 minutes


for iEpoch = 1:3
    
    currEpoch = epochNames{iEpoch};

    fprintf(['\nProcessing ' currEpoch ' ..'])

    
    totalT = endT.(currEpoch) - startT.(currEpoch);
    nBins  = floor((totalT - binDur)/stepDur) + 1;

    binStarts   = startT.(currEpoch) + (0:nBins-1)*stepDur;
    binEnds     = binStarts + binDur;
    binCenters.(currEpoch)  = binStarts + binDur/2;


    unit_fr.(currEpoch) = nan(nUnits, nBins);
    for iBin = 1:nBins
        
        idx = concatPBEbinCenters > binStarts(iBin) & concatPBEbinCenters < binEnds(iBin);
        
        for iUnit = 1:nUnits
          
            unit_fr.(currEpoch)(iUnit, iBin) = numel(find(unit_firings(iUnit, idx)))/(numel(idx) * 0.02);

        end
    end
end
 


%% plots

placeCells = find(peakPFfiring >= 2);
placeCells = placeCells(randperm(numel(placeCells)));

exampleUnits = [81; 24; placeCells(1:18)];



%% plot the matrix of correlation between learned tunnigs in different time windows


plotheight = 15;
plotwidth  = 15;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);



binCenters_avg = s.binCenters;



binCenters_avg.pre = binCenters_avg.pre/3600;
binCenters_avg.run = binCenters_avg.run/3600;
binCenters_avg.post = binCenters_avg.post/3600;


binCenters = [binCenters_avg.pre binCenters_avg.run binCenters_avg.post];
tickLabels = 0:2:binCenters(end);

for iUnit = 1:numel(exampleUnits)
    
    
    currUnit = exampleUnits(iUnit);

    learnedTuning = cell(3,1); 
    for iEpoch = 1:3

        currEpoch = epochNames{iEpoch};

        learnedTuning{iEpoch} = squeeze(assemblyTunings_time.(currEpoch).data(currUnit, :, :, 1));
        learnedTuning{iEpoch} = learnedTuning{iEpoch} ./ repmat(max(learnedTuning{iEpoch}, [], 1), [size(learnedTuning{iEpoch}, 1) 1]);

    end
    
    tmpConcat = cell2mat(learnedTuning');

    subplot(5,4,iUnit)
    
    corrMat = corr(tmpConcat, tmpConcat);
    corrMat(isnan(corrMat)) = 0;
    imagesc(binCenters, binCenters, corrMat)
    caxis([-1 1])
    colormap('jet')
    
    
    line([binCenters_avg.run(1) binCenters_avg.run(1)], binCenters([1 end]), 'linewidth', 1, 'color', 'k')
    line([binCenters_avg.run(end) binCenters_avg.run(end)], binCenters([1 end]), 'linewidth', 1, 'color', 'k')
    line(binCenters([1 end]), [binCenters_avg.run(1) binCenters_avg.run(1)], 'linewidth', 1, 'color', 'k')
    line(binCenters([1 end]), [binCenters_avg.run(end) binCenters_avg.run(end)], 'linewidth', 1, 'color', 'k')

    
    text(0, -0.5, sprintf('unit %d', currUnit), 'fontsize', 8)
    
    set(gca, 'box', 'on', 'fontsize', 8)
    
    xticks(tickLabels)
    yticks(tickLabels)
    axis square
    
end
    
    

    


%%

placeCells = find(peakPFfiring >= 2);

for iUnit = 33%nUnits

    
if ~ismember(iUnit, placeCells)
    continue
end

plotheight = 10;
plotwidth  = 15;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);



% % firing rate of the unit during the PBEs
ax(1) = subplot(14,1,2);
hold on
for iEpoch = 1:3
    currEpoch = epochNames{iEpoch};

    plot(binCenters.(currEpoch), unit_fr.(currEpoch)(iUnit, :), 'linewidth', 2, 'color', [0 0 0 0.7]);
end
xlim([0 binCenters_avg.post(end)])
xticklabels([])


set(gca, 'fontsize', fontsize)


yl = ylim(ax(1));

text(startT.run, yl(1)+2*diff(yl), 'run', 'fontsize', 10)
text(mean([startT.pre endT.pre]), yl(1)+2*diff(yl), 'pre', 'fontsize', 10)
text(mean([startT.post endT.post]), yl(1)+2*diff(yl), 'post', 'fontsize', 10)

text(startT.pre, yl(1)+3*diff(yl),  sprintf('unit-%d', iUnit), 'fontweight', 'bold')




% % learned tunings plotted for individual time winodws

ax(2) = subplot(14,1, 3:7);
hold on

currSpatialTuning = spikes(iUnit).spatialTuning_smoothed.uni;

learnedTuning = cell(3,1); 
for iEpoch = 1:3
    
    currEpoch = epochNames{iEpoch};
    
    learnedTuning{iEpoch} = squeeze(assemblyTunings_time.(currEpoch).data(iUnit, :, :, 1));
    learnedTuning{iEpoch} = learnedTuning{iEpoch} ./ repmat(max(learnedTuning{iEpoch}, [], 1), [size(learnedTuning{iEpoch}, 1) 1]);
    
    imagesc(binCenters_avg.(currEpoch), 1:nPosBins, learnedTuning{iEpoch}); colormap('jet')

end

plot(ax(2), 500*currSpatialTuning + startT.run - diff([startT.run endT.run])/2 , 1:numel(currSpatialTuning), 'linewidth', 2, 'color', [0 0 0 0.5], 'DisplayName', 'track spatial tuning')


yl = ylim(ax(2));

line(500*[0 1] + startT.run - diff([startT.run endT.run])/2 , [yl(2)+5 yl(2)+5], 'linewidth', 2, 'color', [0 0 0 0.5])
text(startT.run - diff([startT.run endT.run])/2, yl(2)+15, '1Hz', 'fontsize', 6, 'color', 'k')

ylim([0 nPosBins+10])


[maxFR, peakPos] = max(currSpatialTuning);

b1 = find(currSpatialTuning < maxFR * 0.25 & (1:nPosBins)' < peakPos, 1, 'last');
b2 = find(currSpatialTuning < maxFR * 0.25 & (1:nPosBins)' > peakPos, 1, 'first');



xlim([0 binCenters_avg.post(end)])
xl = xlim;
try
    line(xl, [b1, b1], 'color', 'w', 'linestyle', '--', 'linewidth', 0.5)
catch
end

try
    line(xl, [b2, b2], 'color', 'w', 'linestyle', '--', 'linewidth', 0.5)
catch
end

xticklabels([])

set(gca, 'fontsize', fontsize)




% % spatial bin correlation between the learned tunings and track spatial bin tuning

ax(3) = subplot(14, 1, 8:9);
hold on

for iEpoch = 1:3
    
    currEpoch = epochNames{iEpoch};
    plot(binCenters_avg.(currEpoch), asTuning_PFcorrelation.(currEpoch).data(iUnit, :), 'linewidth', 2, 'color', [0 0 0 0.7], 'DisplayName', 'avg of multiple time bins');
    
    xl = binCenters_avg.(currEpoch)([1 end]);
    line(xl, [nanmedian(asTuning_PFcorrelation.(currEpoch).ui(iUnit, :)) nanmedian(asTuning_PFcorrelation.(currEpoch).ui(iUnit, :))], 'color', [0 0 0 0.3], 'linewidth', 1, 'DisplayName', 'pooled ui median_indiv bins');
    line(xl, [nanprctile(asTuning_PFcorrelation.(currEpoch).ui(iUnit, :), 25) nanprctile(asTuning_PFcorrelation.(currEpoch).ui(iUnit, :), 25)], 'color', [0 0 0 0.3], 'linewidth', 0.5, 'linestyle', '--')
    line(xl, [nanprctile(asTuning_PFcorrelation.(currEpoch).ui(iUnit, :), 75) nanprctile(asTuning_PFcorrelation.(currEpoch).ui(iUnit, :), 75)], 'color', [0 0 0 0.3], 'linewidth', 0.5, 'linestyle', '--')

end
xticklabels([])
xlim([0 binCenters_avg.post(end)])

set(gca, 'fontsize', fontsize)




% % KL distance: bits of infomration lost when I use learned tunings to
% describ the information encoded in the track spatail tunings

ax(4) = subplot(14, 1, 10:11);
hold on

for iEpoch = 1:3
    
    currEpoch = epochNames{iEpoch};

    plot(binCenters_avg.(currEpoch), asTuning_PF_KLdistance.(currEpoch).data(iUnit, :), 'linewidth', 2, 'color', [0 0 0 0.7], 'DisplayName', 'avg of multiple time bins');
    
    xl = binCenters_avg.(currEpoch)([1 end]);
    line(xl, [nanmedian(asTuning_PF_KLdistance.(currEpoch).ui(iUnit, :)) nanmedian(asTuning_PF_KLdistance.(currEpoch).ui(iUnit, :))], 'color', [0 0 0 0.3], 'linewidth', 1, 'DisplayName', 'pooled ui median: 20msbin-resolved');
    line(xl, [nanprctile(asTuning_PF_KLdistance.(currEpoch).ui(iUnit, :), 25) nanprctile(asTuning_PF_KLdistance.(currEpoch).ui(iUnit, :), 25)], 'color', [0 0 0 0.3], 'linewidth', 0.5, 'linestyle', '--')
    line(xl, [nanprctile(asTuning_PF_KLdistance.(currEpoch).ui(iUnit, :), 75) nanprctile(asTuning_PF_KLdistance.(currEpoch).ui(iUnit, :), 75)], 'color', [0 0 0 0.3], 'linewidth', 0.5, 'linestyle', '--')

end

xticklabels([])
xlim([0 binCenters_avg.post(end)])

set(gca, 'fontsize', fontsize)




% % KS stats: the max difference between two CDFs

ax(5) = subplot(14, 1, 12:13);
hold on

for iEpoch = 1:3
    
    currEpoch = epochNames{iEpoch};

    plot(binCenters_avg.(currEpoch), asTuning_PF_KSstat.(currEpoch).data(iUnit, :), 'linewidth', 2, 'color', [0 0 0 0.7], 'DisplayName', 'avg of multiple time bins');
    
    xl = binCenters_avg.(currEpoch)([1 end]);
    line(xl, [nanmedian(asTuning_PF_KSstat.(currEpoch).ui(iUnit, :)) nanmedian(asTuning_PF_KSstat.(currEpoch).ui(iUnit, :))], 'color', [0 0 0 0.3], 'linewidth', 1, 'DisplayName', 'pooled ui median: 20msbin-resolved');
    line(xl, [nanprctile(asTuning_PF_KSstat.(currEpoch).ui(iUnit, :), 25) nanprctile(asTuning_PF_KSstat.(currEpoch).ui(iUnit, :), 25)], 'color', [0 0 0 0.3], 'linewidth', 0.5, 'linestyle', '--')
    line(xl, [nanprctile(asTuning_PF_KSstat.(currEpoch).ui(iUnit, :), 75) nanprctile(asTuning_PF_KSstat.(currEpoch).ui(iUnit, :), 75)], 'color', [0 0 0 0.3], 'linewidth', 0.5, 'linestyle', '--')

end


xlim([0 binCenters_avg.post(end)])

set(gca, 'fontsize', fontsize)


linkaxes(ax, 'x')

ylabel(ax(1), 'PBE fr(Hz)', 'fontsize', fontsize)
ylabel(ax(2), 'position bin', 'fontsize', fontsize)
ylabel(ax(3), 'correlation', 'fontsize', fontsize)
ylabel(ax(4), 'KL div.', 'fontsize', fontsize)
ylabel(ax(5), 'KS stat', 'fontsize', fontsize)

xlabel(ax(5), 'time(sec)', 'fontsize', fontsize)

if iUnit < 10
    fileName = sprintf('%s_unit%d%d.pdf', sessionName, 0, iUnit);
else
    fileName = sprintf('%s_unit%d.pdf', sessionName, iUnit);
end
    
print(gcf, fileName, '-dpdf', '-painters')

close all
end



end

%%
function output = nanprctile(x, percent)

output = prctile(x(~isnan(x)), percent);

end


function drawpatch(x,y, color, faceAlpha)

idx = isnan(y);
y(idx) = 0;

patch(x, y, color, 'EdgeColor', 'none', 'FaceAlpha', faceAlpha)

end

