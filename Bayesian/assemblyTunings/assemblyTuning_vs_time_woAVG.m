function assemblyTuning_vs_time_woAVG(sessionNumber)


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
storagePath = fullfile(basePath, 'assemblyTunings');
mkdir(storagePath)


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


ifshuffle = 0;
[AT_prefLocation.data, AT_placeFieldCorr.data, AT_placeFieldKL.data, unit_firings, concatPBEbinCenters] = calculateAssemblyTuning_woAVG(PBEInfo, spikes, clusterQuality, ifshuffle); % 


% 
% ifshuffle = 1;
% [AT_prefLocation.ui, AT_placeFieldCorr.ui, AT_placeFieldKL.ui, AT_placeFieldKL2.ui] = calculateAssemblyTuning_woAVG(PBEInfo, spikes, clusterQuality, ifshuffle); % 



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
%         asTuning_PF_KLdistance2.(currEpoch).(curDataType) = zeros(nUnits, nTimeBins);
        
        asTuning_PFcorrelation.(currEpoch).(curDataType) = zeros(nUnits, nTimeBins);

        for iUnit = 1:nUnits

            learnedTuning = squeeze(assemblyTunings_time.(currEpoch).(curDataType)(iUnit, :, :, 1));
            learnedTuning = learnedTuning ./ repmat(sum(learnedTuning, 1), [size(learnedTuning, 1) 1]);
            learnedTuning = learnedTuning + eps;

            asTuning_PF_KLdistance.(currEpoch).(curDataType)(iUnit, :)  = sum(learnedTuning .* (log(learnedTuning) - log(repmat(spatialTuning(:, iUnit), [1 size(learnedTuning, 2)]))), 1);
%             asTuning_PF_KLdistance2.(currEpoch).(curDataType)(iUnit, :) = sum(repmat(spatialTuning(:, iUnit), [1 size(learnedTuning, 2)]) .* (log(repmat(spatialTuning(:, iUnit), [1 size(learnedTuning, 2)])) - log(learnedTuning)), 1);
            
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
gw  = gausswindow(3, 9);

posBins  = 1:nPosBins;
corrBins = -1:0.02:1;



for iEpoch = 1:3
    
%     iEpoch = 3;
    currEpoch = epochNames{iEpoch};

    fprintf(['\nProcessing ' currEpoch ' ..'])

    
    totalT = endT.(currEpoch) - startT.(currEpoch);
    nBins  = floor((totalT - binDur)/stepDur) + 1;

    binStarts   = startT.(currEpoch) + (0:nBins-1)*stepDur;
    binEnds     = binStarts + binDur;
    binCenters.(currEpoch)  = binStarts + binDur/2;


    unit_fr.(currEpoch) = nan(nUnits, nBins);
    
    prefLocation_med.(currEpoch)  = nan(nUnits, nBins);
    prefLocation_lq.(currEpoch)   = nan(nUnits, nBins);
    prefLocation_hq.(currEpoch)   = nan(nUnits, nBins);
    prefLocation_mode.(currEpoch) = nan(nUnits, nBins);
    KL_prefLocation.(currEpoch)   = nan(nUnits, nBins);
    
    
    placeFieldCorr_med.(currEpoch) = nan(nUnits, nBins);
    placeFieldCorr_lq.(currEpoch) = nan(nUnits, nBins);
    placeFieldCorr_hq.(currEpoch) = nan(nUnits, nBins);
    KL_placeFieldCorr.(currEpoch) = nan(nUnits, nBins);
    corrSign_PFcorrelation.(currEpoch) = nan(nUnits, nBins);   
    
    
    placeFieldKL_med.(currEpoch) = nan(nUnits, nBins);
    placeFieldKL_lq.(currEpoch) = nan(nUnits, nBins);
    placeFieldKL_hq.(currEpoch) = nan(nUnits, nBins);
    KL_currPlaceFieldKL.(currEpoch) = nan(nUnits, nBins);
    corrSign_KLdistance.(currEpoch) = nan(nUnits, nBins);

    
    % %
    currPrefLocation_ui = AT_prefLocation.ui; 
    currPrefLocation_ui = currPrefLocation_ui(:); 
    
    pVect_PrefLocation_ui = hist(currPrefLocation_ui, posBins);
    pVect_PrefLocation_ui = pVect_PrefLocation_ui ./sum(pVect_PrefLocation_ui);
    pVect_PrefLocation_ui = conv(pVect_PrefLocation_ui, gw, 'same');
    pVect_PrefLocation_ui = pVect_PrefLocation_ui + eps;

    
    % %
    currPlaceFieldCorr_ui = AT_placeFieldCorr.ui;
    currPlaceFieldCorr_ui = currPlaceFieldCorr_ui(:);
    
    pVect_PlaceFieldCorr_ui = hist(currPlaceFieldCorr_ui, corrBins);
    pVect_PlaceFieldCorr_ui = pVect_PlaceFieldCorr_ui ./sum(pVect_PlaceFieldCorr_ui);
    pVect_PlaceFieldCorr_ui = conv(pVect_PlaceFieldCorr_ui, gw, 'same');
    pVect_PlaceFieldCorr_ui = pVect_PlaceFieldCorr_ui + eps;
    
   
    
    % %
    currPlaceFieldKL_ui = AT_placeFieldKL.ui;
    currPlaceFieldKL_ui = currPlaceFieldKL_ui(:);
    
    pooled = [AT_placeFieldKL.data(:); currPlaceFieldKL_ui];
    
    klBins = linspace(min(pooled), max(pooled), 100);
    
    pVect_PlaceFieldKL_ui = hist(currPlaceFieldKL_ui, klBins);
    pVect_PlaceFieldKL_ui = pVect_PlaceFieldKL_ui ./sum(pVect_PlaceFieldKL_ui);
    pVect_PlaceFieldKL_ui = conv(pVect_PlaceFieldKL_ui, gw, 'same');
    pVect_PlaceFieldKL_ui = pVect_PlaceFieldKL_ui + eps;
    
    
    
%     placeFieldCorr_mode.(currEpoch) = nan(nUnits, nBins);
    
    for iBin = 1:nBins
        
%         iBin = 6;
        idx = concatPBEbinCenters > binStarts(iBin) & concatPBEbinCenters < binEnds(iBin);
        
        for iUnit = 1:nUnits
%             iUnit = 55;

%             idx2 = AT_prefLocation(iUnit).tvec > binStarts(iBin) & AT_prefLocation(iUnit).tvec < binEnds(iBin);
            
            currPrefLocation_data = AT_prefLocation.data(iUnit, idx);
            
            unit_fr.(currEpoch)(iUnit, iBin) = numel(find(unit_firings(iUnit, idx)))/(numel(idx) * 0.02);

            prefLocation_med.(currEpoch)(iUnit, iBin) = nanmedian(currPrefLocation_data);
            prefLocation_lq.(currEpoch)(iUnit, iBin)  = nanprctile(currPrefLocation_data, 25);
            prefLocation_hq.(currEpoch)(iUnit, iBin)  = nanprctile(currPrefLocation_data, 75);
            
            
            pVect_data = hist(currPrefLocation_data, posBins);
            pVect_data = pVect_data ./sum(pVect_data);            
            pVect_data = conv(pVect_data, gw, 'same');
            pVect_data = pVect_data + eps;
            
            
            if ~isempty(pVect_data)
                [~, prefLocation_mode.(currEpoch)(iUnit, iBin)] = max(pVect_data);
            end
            
            % the kl divergence
            KL_prefLocation.(currEpoch)(iUnit, iBin) = sum(pVect_data .* (log2(pVect_data) - log2(pVect_PrefLocation_ui)));
           
            
            % %
            currPlaceFieldCorr_data = AT_placeFieldCorr.data(iUnit, idx);
            

            placeFieldCorr_med.(currEpoch)(iUnit, iBin) = nanmedian(currPlaceFieldCorr_data);
            placeFieldCorr_lq.(currEpoch)(iUnit, iBin)  = nanprctile(currPlaceFieldCorr_data, 25);
            placeFieldCorr_hq.(currEpoch)(iUnit, iBin)  = nanprctile(currPlaceFieldCorr_data, 75);
            
            try
                [~, temp] = ranksum(currPlaceFieldCorr_data, AT_placeFieldCorr.ui(iUnit, :), 'tail', 'right', 'alpha', 0.05);

                if temp == 1
                    corrSign_PFcorrelation.(currEpoch)(iUnit, iBin) = 1;
                end
            catch
            end
            
            corrBins = -1:0.02:1;
            pVect_data = hist(currPlaceFieldCorr_data, corrBins);
            pVect_data = pVect_data ./sum(pVect_data);
            pVect_data = conv(pVect_data, gw, 'same');
            pVect_data = pVect_data + eps;
            
            % the kl divergence
            KL_placeFieldCorr.(currEpoch)(iUnit, iBin) = sum(pVect_data .* (log2(pVect_data) - log2(pVect_PlaceFieldCorr_ui)));
           
            
            
            % %
            currPlaceFieldKL_data = AT_placeFieldKL.data(iUnit, idx);
            
            placeFieldKL_med.(currEpoch)(iUnit, iBin) = nanmedian(currPlaceFieldKL_data);
            placeFieldKL_lq.(currEpoch)(iUnit, iBin)  = nanprctile(currPlaceFieldKL_data, 25);
            placeFieldKL_hq.(currEpoch)(iUnit, iBin)  = nanprctile(currPlaceFieldKL_data, 75);
            
            try
                [~, temp] = ranksum(currPlaceFieldKL_data, AT_placeFieldKL.ui(iUnit, :), 'tail', 'left', 'alpha', 0.05);

                if temp == 1
                    corrSign_KLdistance.(currEpoch)(iUnit, iBin) = 1;
                end
            catch
            end
            
            
            
            pVect_data = hist(currPlaceFieldKL_data, klBins);
            pVect_data = pVect_data ./sum(pVect_data);
            pVect_data = conv(pVect_data, gw, 'same');
            pVect_data = pVect_data + eps;
            
            
            % the kl divergence
            KL_currPlaceFieldKL.(currEpoch)(iUnit, iBin) = sum(pVect_data .* (log2(pVect_data) - log2(pVect_PlaceFieldKL_ui)));

        end
    end
end


%%

placeCells = find(peakPFfiring >= 3);

for iUnit = 67%nUnits

    
if ~ismember(iUnit, placeCells)
    continue
end

plotheight = 10;
plotwidth  = 15;
fontsize   = 8;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);



ax(1) = subplot(21,1,1);
hold on
for iEpoch = 1:3
    currEpoch = epochNames{iEpoch};

    plot(binCenters.(currEpoch), unit_fr.(currEpoch)(iUnit, :), 'linewidth', 2, 'color', [0 0 0 0.7]);
end
xlim([0 binCenters.post(end)])
xticklabels([])


yl = ylim(ax(1));

text(startT.run, yl(1)+2*diff(yl), 'run', 'fontsize', 10)
text(mean([startT.pre endT.pre]), yl(1)+2*diff(yl), 'pre', 'fontsize', 10)
text(mean([startT.post endT.post]), yl(1)+2*diff(yl), 'post', 'fontsize', 10)

text(startT.pre, yl(1)+2*diff(yl),  sprintf('unit-%d', iUnit), 'fontweight', 'bold')

% %
ax(5) = subplot(21,1, 2:6);
hold on

for iEpoch = 1:3
    
    currEpoch = epochNames{iEpoch};
    
    learnedTuning = squeeze(assemblyTunings_time.(currEpoch).data(iUnit, :, :, 1));
    learnedTuning = learnedTuning ./ repmat(max(learnedTuning, [], 1), [size(learnedTuning, 1) 1]);
    
    imagesc(binCenters_avg.(currEpoch), 1:nPosBins, learnedTuning); colormap('jet')

end


% %

ax(2) = subplot(21,1, 7:11);
hold on

for iEpoch = 1:3
    
    currEpoch = epochNames{iEpoch};
    

    h1 = plot(binCenters.(currEpoch), prefLocation_med.(currEpoch)(iUnit, :), 'linewidth', 2, 'color', [0 0 0 0.7], 'DisplayName', '20msbin-resolved median');
    
    
    drawpatch([binCenters.(currEpoch) fliplr(binCenters.(currEpoch))], [prefLocation_lq.(currEpoch)(iUnit, :) fliplr(prefLocation_hq.(currEpoch)(iUnit, :))], [.5 .5 .5], 0.2)
    
    
    h2 = plot(binCenters_avg.(currEpoch), AT_prefLocation_avg.(currEpoch).data(iUnit, :), 'linewidth', 2, 'color', [[255 140 0]/255 0.7], 'DisplayName', 'avg learned tuning');
    
    h3 = plot(binCenters.(currEpoch), prefLocation_mode.(currEpoch)(iUnit, :), 'linewidth', 2, 'color', [0 1 0 0.5], 'DisplayName', '20msbin-resolved mode');
%     h4 = plot(binCenters.(currEpoch), KL_prefLocation.(currEpoch)(iUnit, :)*60, 'linewidth', 2, 'color', [0 1 0 0.5], 'DisplayName', 'KL distance');

end
xlim([0 binCenters.post(end)])
% xl = xlim;
% line(xl, [nanmedian(AT_prefLocation.data(:)) nanmedian(AT_prefLocation.data(:))], 'color', [1 0 0 0.3], 'linewidth', 1)
% patch([xl fliplr(xl)], [nanprctile(AT_prefLocation.data(:), 25) nanprctile(AT_prefLocation.data(:), 25) nanprctile(AT_prefLocation.data(:), 75) nanprctile(AT_prefLocation.data(:), 75)], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1)

% line(patch(AT_prefLocation.data(iUnit, idx);))

h4 = plot(ax(2), 500*spikes(iUnit).spatialTuning_smoothed.uni + startT.run - diff([startT.run endT.run])/2 , 1:numel(spikes(iUnit).spatialTuning_smoothed.uni), 'linewidth', 2, 'color', [0 0 1 0.5], 'DisplayName', 'track spatial tuning');

xticklabels([])


yl = ylim(ax(2));

line(ax(2), 500*[0 1] + startT.run - diff([startT.run endT.run])/2 , [yl(1)-0.2*diff(yl) yl(1)-0.2*diff(yl)], 'linewidth', 2, 'color', [0 0 1 0.5])
text(startT.run - diff([startT.run endT.run])/2, yl(1)-0.1*diff(yl), '1Hz', 'fontsize', 6, 'color', 'b')


oldPos = get(gca, 'position');

tt = legend([h1, h3, h2, h4], 'location', 'northout', 'box', 'off', 'fontsize', 6);

oldPos2 = get(tt, 'position');

set(tt, 'position', [oldPos2(1)+ 0.3 oldPos2(2)+0.12 oldPos2(3:4)])

set(gca, 'position', oldPos)




ax(3) = subplot(21, 1, 12:16);
hold on

for iEpoch = 1:3
    
    currEpoch = epochNames{iEpoch};

    plot(binCenters.(currEpoch), placeFieldCorr_med.(currEpoch)(iUnit, :), 'linewidth', 2, 'color', [0 0 0 0.7], 'DisplayName', 'indiv time bins median');
    
    drawpatch([binCenters.(currEpoch) fliplr(binCenters.(currEpoch))], [placeFieldCorr_lq.(currEpoch)(iUnit, :) fliplr(placeFieldCorr_hq.(currEpoch)(iUnit, :))], [.5 .5 .5],  0.2)
    
    plot(binCenters_avg.(currEpoch), asTuning_PFcorrelation.(currEpoch).data(iUnit, :), 'linewidth', 2, 'color', [[255 140 0]/255 0.7], 'DisplayName', 'avg of multiple time bins');
    
%     plot(binCenters.(currEpoch), KL_placeFieldCorr.(currEpoch)(iUnit, :), 'linewidth', 2, 'color', [0 1 0 0.5], 'DisplayName', 'KL distance');

%     h3 = plot(binCenters.(currEpoch), placeFieldCorr_mode.(currEpoch)(iUnit, :), 'linewidth', 2, 'color', [1 0 0 0.5], 'DisplayName', 'indiv time bins mode');

    plot(binCenters.(currEpoch), 0.95*corrSign_PFcorrelation.(currEpoch)(iUnit, :), '.', 'markersize', 5, 'color', [.7 .7 .7 .7])

end
xticklabels([])
xlim([0 binCenters.post(end)])
xl = xlim;
line(xl, [nanmedian(AT_placeFieldCorr.ui(iUnit, :)) nanmedian(AT_placeFieldCorr.ui(iUnit, :))], 'color', [0 0 0 0.3], 'linewidth', 1)
line(xl, [nanprctile(AT_placeFieldCorr.ui(iUnit, :), 25) nanprctile(AT_placeFieldCorr.ui(iUnit, :), 25)], 'color', [0 0 0 0.3], 'linewidth', 0.5, 'linestyle', '--')
line(xl, [nanprctile(AT_placeFieldCorr.ui(iUnit, :), 75) nanprctile(AT_placeFieldCorr.ui(iUnit, :), 75)], 'color', [0 0 0 0.3], 'linewidth', 0.5, 'linestyle', '--')
% patch([xl fliplr(xl)], [nanprctile(AT_placeFieldCorr.ui(iUnit, :), 25) nanprctile(AT_placeFieldCorr.ui(iUnit, :), 25) nanprctile(AT_placeFieldCorr.ui(iUnit, :), 75) nanprctile(AT_placeFieldCorr.ui(iUnit, :), 75)], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1)



line(xl, [nanmedian(asTuning_PFcorrelation.(currEpoch).ui(iUnit, :)) nanmedian(asTuning_PFcorrelation.(currEpoch).ui(iUnit, :))], 'color', [[255 140 0]/255 0.7], 'linewidth', 1, 'DisplayName', 'pooled ui median_indiv bins');
line(xl, [nanprctile(asTuning_PFcorrelation.(currEpoch).ui(iUnit, :), 25) nanprctile(asTuning_PFcorrelation.(currEpoch).ui(iUnit, :), 25)], 'color', [[255 140 0]/255 0.7], 'linewidth', 0.5, 'linestyle', '--')
line(xl, [nanprctile(asTuning_PFcorrelation.(currEpoch).ui(iUnit, :), 75) nanprctile(asTuning_PFcorrelation.(currEpoch).ui(iUnit, :), 75)], 'color', [[255 140 0]/255 0.7], 'linewidth', 0.5, 'linestyle', '--')

% patch([xl fliplr(xl)], [nanprctile(asTuning_PFcorrelation.(currEpoch).ui(iUnit, :), 25) nanprctile(asTuning_PFcorrelation.(currEpoch).ui(iUnit, :), 25) nanprctile(asTuning_PFcorrelation.(currEpoch).ui(iUnit, :), 75) nanprctile(asTuning_PFcorrelation.(currEpoch).ui(iUnit, :), 75)], 'b' , 'EdgeColor', 'none', 'FaceAlpha', 0.1)


% 
% xl = xlim;
% line(xl, [0 0], 'color', [0 0 0 0.5], 'linewidth', 0.6)



ax(4) = subplot(21, 1, 17:21);
hold on

for iEpoch = 1:3
    
    currEpoch = epochNames{iEpoch};

    plot(binCenters.(currEpoch), placeFieldKL_med.(currEpoch)(iUnit, :), 'linewidth', 2, 'color', [0 0 0 0.7], 'DisplayName', 'indiv time bins median');
    
    drawpatch([binCenters.(currEpoch) fliplr(binCenters.(currEpoch))], [placeFieldKL_lq.(currEpoch)(iUnit, :) fliplr(placeFieldKL_hq.(currEpoch)(iUnit, :))], [.5 .5 .5], 0.2)
    
    plot(binCenters_avg.(currEpoch), asTuning_PF_KLdistance.(currEpoch).data(iUnit, :), 'linewidth', 2, 'color', [[255 140 0]/255 0.7], 'DisplayName', 'avg of multiple time bins');
    plot(binCenters_avg.(currEpoch), asTuning_PF_KSstat.(currEpoch).data(iUnit, :)*10, 'linewidth', 2, 'color', [0 0 1 0.7], 'DisplayName', 'avg of multiple time bins');

  
%     plot(binCenters.(currEpoch), KL_placeFieldCorr.(currEpoch)(iUnit, :), 'linewidth', 2, 'color', [0 1 0 0.5], 'DisplayName', 'KL distance');

end

yl = ylim;

for iEpoch = 1:3
    currEpoch = epochNames{iEpoch};
    plot(binCenters.(currEpoch), 0.95*yl(2)*corrSign_KLdistance.(currEpoch)(iUnit, :), '.', 'markersize', 5, 'color', [.7 .7 .7 .7])
end



xlim([0 binCenters.post(end)])
xl = xlim;

h1 = line(xl, [nanmedian(AT_placeFieldKL.ui(iUnit, :)) nanmedian(AT_placeFieldKL.ui(iUnit, :))], 'color', [0 0 0 0.3], 'linewidth', 1, 'DisplayName', 'pooled ui median: avg learned tuning');
line(xl, [nanprctile(AT_placeFieldKL.ui(iUnit, :), 25) nanprctile(AT_placeFieldKL.ui(iUnit, :), 25)], 'color', [0 0 0 0.3], 'linewidth', 0.5, 'linestyle', '--')
line(xl, [nanprctile(AT_placeFieldKL.ui(iUnit, :), 75) nanprctile(AT_placeFieldKL.ui(iUnit, :), 75)], 'color', [0 0 0 0.3], 'linewidth', 0.5, 'linestyle', '--')

% patch([xl fliplr(xl)], [nanprctile(AT_placeFieldKL.ui(iUnit, :), 25) nanprctile(AT_placeFieldKL.ui(iUnit, :), 25) nanprctile(AT_placeFieldKL.ui(iUnit, :), 75) nanprctile(AT_placeFieldKL.ui(iUnit, :), 75)], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1)

h2 = line(xl, [nanmedian(asTuning_PF_KLdistance.(currEpoch).ui(iUnit, :)) nanmedian(asTuning_PF_KLdistance.(currEpoch).ui(iUnit, :))], 'color', [[255 140 0]/255 0.7], 'linewidth', 1, 'DisplayName', 'pooled ui median: 20msbin-resolved');
line(xl, [nanprctile(asTuning_PF_KLdistance.(currEpoch).ui(iUnit, :), 25) nanprctile(asTuning_PF_KLdistance.(currEpoch).ui(iUnit, :), 25)], 'color', [[255 140 0]/255 0.7], 'linewidth', 0.5, 'linestyle', '--')
line(xl, [nanprctile(asTuning_PF_KLdistance.(currEpoch).ui(iUnit, :), 75) nanprctile(asTuning_PF_KLdistance.(currEpoch).ui(iUnit, :), 75)], 'color', [[255 140 0]/255 0.7], 'linewidth', 0.5, 'linestyle', '--')


h3 = line(xl, [nanmedian(asTuning_PF_KSstat.(currEpoch).ui(iUnit, :)) nanmedian(asTuning_PF_KSstat.(currEpoch).ui(iUnit, :))]*10, 'color', [0 0 1 0.7], 'linewidth', 1, 'DisplayName', 'pooled ui median: 20msbin-resolved');
line(xl, [nanprctile(asTuning_PF_KSstat.(currEpoch).ui(iUnit, :), 25) nanprctile(asTuning_PF_KSstat.(currEpoch).ui(iUnit, :), 25)]*10, 'color', [0 0 1 0.7], 'linewidth', 0.5, 'linestyle', '--')
line(xl, [nanprctile(asTuning_PF_KSstat.(currEpoch).ui(iUnit, :), 75) nanprctile(asTuning_PF_KSstat.(currEpoch).ui(iUnit, :), 75)]*10, 'color', [0 0 1 0.7], 'linewidth', 0.5, 'linestyle', '--')

% patch([xl fliplr(xl)], [nanprctile(asTuning_PF_KLdistance.(currEpoch).ui(iUnit, :), 25) nanprctile(asTuning_PF_KLdistance.(currEpoch).ui(iUnit, :), 25) nanprctile(asTuning_PF_KLdistance.(currEpoch).ui(iUnit, :), 75) nanprctile(asTuning_PF_KLdistance.(currEpoch).ui(iUnit, :), 75)], 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.1)


oldPos = get(gca, 'position');
legend([h1, h2], 'location', 'northeast', 'box', 'off', 'fontsize', 6);
set(gca, 'position', oldPos)




linkaxes(ax, 'x')

ylabel(ax(1), 'PBE fr(Hz)', 'fontsize', fontsize)
ylabel(ax(2), 'pref. location', 'fontsize', fontsize)
ylabel(ax(3), 'PF correlation', 'fontsize', fontsize)
ylabel(ax(4), 'KL distance', 'fontsize', fontsize)
ylabel(ax(5), 'position bin', 'fontsize', fontsize)

xlabel(ax(3), 'time(sec)', 'fontsize', fontsize)

if iUnit < 10
    fileName = sprintf('%s_unit%d%d.pdf', sessionName, 0, iUnit);
else
    fileName = sprintf('%s_unit%d.pdf', sessionName, iUnit);
end
    
print(gcf, fileName, '-dpdf', '-painters')

% close all
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

