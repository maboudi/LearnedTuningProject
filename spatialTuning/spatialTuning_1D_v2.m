function [spatialTunings, PF_sorted, activeUnits_sorted] = spatialTuning_1D_v2(spikeStruct, qclus, fileinfo, behavior, speed, direction, posBinSize, runSpeedThresh, flaggedUnits, timeUnit, FileBase)

% EXPLAIN THE FUNCTION

totoalNumofUnits = max(spikeStruct.unit);

% including only the spikes meeting certain criteria

if strcmp(direction, 'LR')
    desiredMod = 0; % even traversals
elseif strcmp(direction, 'RL')
    desiredMod = 1; % odd traversals
elseif strcmp(direction, 'uni')
    desiredMod = [0 1];
end


totalNumofLaps   = max(spikeStruct.lap);
currDirlaps      = find(ismember(mod(1:totalNumofLaps, 2), desiredMod));
currDirNumoflaps = numel(currDirlaps);


linearPos = fileinfo.xyt2(:, 1);
nPosBins  = floor((max(linearPos) - min(linearPos))/posBinSize);

nIter      = 1000;
subsetPrct = 90;

subsetSize = floor(numel(currDirlaps) * subsetPrct/100);

spatialTunings_iters = zeros(totoalNumofUnits, nPosBins, subsetSize);


for ii = 1:nIter
    printf('.')
    % subset of all laps that is used in the current calculation
    
    currLapSubset = currDirlaps(randperm(currDirNumoflaps, subsetSize));
    [spatialTunings_iters(:,:,ii), linearPoscenters] = calSpatialTuning(spikeStruct, qclus, fileinfo, behavior, currLapSubset, speed, posBinSize, runSpeedThresh, flaggedUnits, timeUnit); % non-smoothed tunings
    
end

nonSmoothed_spatialTunings = min(spatialTunings_iters, [], 3);

spatialTunings = zeros(size(nonSmoothed_spatialTunings));
peakPosBin     = zeros(totoalNumofUnits, 1);
peakRates      = zeros(totoalNumofUnits, 1);

win = gausswindow(3,7);

for unit = 1:totoalNumofUnits
    
    spatialTunings(unit , :) = conv(nonSmoothed_spatialTunings(unit, :), win, 'same');    
    [peakRates(unit), peakPosBin(unit)] = max(spatialTunings(unit , :));
  
end



% plot the place fields

activeUnits = find(sum(spatialTunings, 2) > 0); % units with non-zero tunings

peakPosBin   = peakPosBin(activeUnits); 

[~, sortIdx] = sort(peakPosBin, 'ascend');
peakRates_sorted = peakRates(activeUnits(sortIdx)); % maximum firing rate

PF_sorted = spatialTunings(activeUnits(sortIdx), :); % sort the place fields based on the peaks


activeUnits_sorted = activeUnits(sortIdx); % keep track of the unit labels

PF_sorted_norm = PF_sorted ./ repmat(max(PF_sorted, [], 2), [1 size(PF_sorted, 2)]); % Normalize the peaks to one for visulaization



% units with peak firing rates below the threshold (2 Hz) will be shown in black

figure;

x0=0;
y0=0;
width=400;
height=400* size(PF_sorted_norm, 1)/20;
difpos = linearPoscenters(2)- linearPoscenters(1);
set(gcf,'units','points','position',[x0,y0,width,height])

tt = 0;

for jj = 1 : size(PF_sorted_norm, 1)
    
    if peakRates_sorted(jj) > 1
        
        tt = tt + 1;
        
        if ismember(activeUnits_sorted(jj), flaggedUnits)
            cl = [238,130,238]/255; % violet
        else 
            cl = 'r';
        end
        
        fill([linearPoscenters fliplr(linearPoscenters)], [0.06*tt+PF_sorted_norm(jj, :)/20 fliplr(0.06*tt*ones(size(PF_sorted_norm(jj, :))))], cl,'LineStyle','none')
        hold on

        plot(linearPoscenters, 0.06*tt+PF_sorted_norm(jj, :)/20,'color', 'k','linewidth', 0.5);

        alpha(0.5)
        
%         sparsity(tt) = mean(PF_sorted_norm(jj, :))^2 / mean(PF_sorted_norm(jj, :).^2);
%         peakBins(tt) = peakPosBin_sorted(jj) * posBinSize;
        
    end
    
    

    set(gca, 'YTick', [], 'YTickLabel', [], 'color', 'none', 'YColor', 'none', 'box', 'off')
%     text(linearPoscenters(1)-5*difpos, 0.06*jj, num2str(activeUnits_sorted(jj)), 'fontsize', 7, 'HorizontalAlignment', 'center');
%     text(linearPoscenters(end)+5*difpos, 0.06*jj, sprintf('%.1f Hz', peakRates_sorted(jj)), 'fontsize', 7, 'HorizontalAlignment', 'center');
end


% tt = tt + 3;
% cl = 'k';
% 
% fill([linearPoscenters fliplr(linearPoscenters)], [0.06*tt+sum_PF/2 fliplr(0.06*tt*ones(size(sum_PF)))], cl,'LineStyle','none')
% 
% hold on
% plot(linearPoscenters, 0.06*tt+sum_PF/2,'color', 'k','linewidth', 0.5);
% 
% alpha(0.5)


xlim([linearPoscenters(1)-difpos/2 linearPoscenters(end)+difpos/2])

xlabel('Position on track (cm)', 'fontsize', 10)

h = text(linearPoscenters(1)-5*difpos, 0.06*(tt)/2, 'Unit', 'fontsize', 10, 'HorizontalAlignment', 'center');
set(h, 'rotation', 90)
% 
% h = text(linearPoscenters(1)-5*difpos, 0.06*tt, 'Norm pooled', 'fontsize', 10, 'HorizontalAlignment', 'center');
% set(h, 'rotation', 90)



if ~strcmp(direction, 'uni')

    ha = annotation('arrow');  
    ha.Parent = gca;  

    if strcmp(direction, 'LR')
        ha.X = [linearPoscenters(1) linearPoscenters(15)]; 
    elseif strcmp(direction, 'RL')
        ha.X = [linearPoscenters(end) linearPoscenters(end-15)]; 
    end

    ha.Y = [0.06*(tt+2) 0.06*(tt+2)];   

    ha.LineWidth  = 2;          % make the arrow bolder for the picture
    ha.HeadWidth  = 10;
    ha.HeadLength = 10;
end



if ~isempty(FileBase)

    if ~strcmp(direction, 'uni')

        filename = [fileinfo.name '_spatialTunings1D_' direction];

    else
        filename = [fileinfo.name '_spatialTunings1D_uni'];
    end


    savepdf(gcf, fullfile(FileBase, filename), '-dpng')
    savefig(gcf, fullfile(FileBase, filename))
end


end




function [spatialTunings, linearPoscenters] = calSpatialTuning(spikeStruct, qclus, fileinfo, behavior, currLapSubset, speed, posBinSize, runSpeedThresh, flaggedUnits, timeUnit)

totoalNumofUnits = max(spikeStruct.unit);

spikeInd = find(spikeStruct.t >= behavior.time(2,1) & spikeStruct.t < behavior.time(2,2) ... % within the RUN period
                & ismember(spikeStruct.qclu, qclus) ... % only stable pyramidal units are included
                & spikeStruct.speed > runSpeedThresh ... % only the spikes happening when the velocity of the animal is higher than the threshold (usually 10 cm/sec)
                & spikeStruct.theta == 1 ... % only spike during high theta power
                & spikeStruct.lap > 0 & ismember(spikeStruct.lap, currLapSubset)); % limiting to the travels in specific direction
                
               
spikePositions = spikeStruct.linearPos(spikeInd);
spikeUnit      = spikeStruct.unit(spikeInd);

activeUnits    = unique(spikeUnit); 


% finding position samples pertaining to the maze period and certain direction of travel 

speedatPos = interp1(speed.t, speed.v, fileinfo.xyt(:, 3));

positionIdx    = find(fileinfo.xyt(:, 3) > behavior.time(2,1) & fileinfo.xyt(:, 3) < behavior.time(2,2) ... % within the run period
                    & speedatPos > runSpeedThresh ... % position bins when animal's velocity is higher than threshold
                    & fileinfo.xyt2(:, 3) > 0 & ismember(fileinfo.xyt2(:, 3), currLapSubset)); % limiting to the travels in specific direction
           

linearPos = fileinfo.xyt2(:, 1);                    
directionLinearPos = nan(size(linearPos));          
directionLinearPos(positionIdx) = linearPos(positionIdx);


% defining the position bins

nPosBins = floor((max(linearPos) - min(linearPos))/posBinSize);
posBinEdges = min(linearPos): posBinSize: max(linearPos); % center of the position bins
linearPoscenters = posBinEdges(1:end-1) + posBinSize/2;
posSamplingPeriod = median(diff(fileinfo.xyt(:, 3)))/timeUnit; % 1/sampling frequency - Hiro's dataset: timeunit is microsecond


% spatialTunings = zeros(length(activeUnits), nPosBins);
% peakPosBin = zeros(length(activeUnits), 1);

totalNumSpikes = zeros(totoalNumofUnits, 1);
% combinedflag   = zeros(length(activeUnits), 1);
spatialTunings = zeros(totoalNumofUnits, nPosBins);

for ii = 1: length(activeUnits)
    
    unit                  = activeUnits(ii);
    
    unitSpikePositions    = spikePositions(spikeUnit == unit);
    
    totalNumSpikes(unit) = length(unitSpikePositions);
    
%     if totalNumSpikes(unit) <= 10
%         continue
%     end
    
    posBinnSpikeCnts      = histc(unitSpikePositions, posBinEdges);
    posBinnSpikeCnts(end) = [];
    
    posBinDwelltime       = histc(directionLinearPos, posBinEdges) * posSamplingPeriod;
    posBinDwelltime(end)  = [];
    
    
    
    if size(posBinnSpikeCnts, 2) > 1
        posBinnSpikeCnts = posBinnSpikeCnts';
    end
    
    nonSmoothed_tuning = posBinnSpikeCnts ./ posBinDwelltime;

    nonSmoothed_tuning(isnan(nonSmoothed_tuning)) = 0;
    nonSmoothed_tuning(isinf(nonSmoothed_tuning)) = 0;
    
    spatialTunings(activeUnits(ii), :) = nonSmoothed_tuning;
    
end


end