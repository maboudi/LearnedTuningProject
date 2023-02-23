function     PCtuning(spike, qclus, fileinfo, whichPart, speedThresh, posBinSize)

% PCtuning calculates the firing rate of each unit within position bins of
% certain size. To do this, number of spikes within each bin is divided by
% the bin's occupancy time. 

Fs = fileinfo.Fs; % sampling frequency

%%% Load position data

tbegin = fileinfo.tbegin;
tend = fileinfo.tend;


if ~isfield(fileinfo, 'pix2cm')
    pix2cm = 1;
else
    pix2cm = fileinfo.pix2cm;
end


xyt = fileinfo.xyt;


withinRange = find(xyt(:,3) >= tbegin & xyt(:,3) <= tend);
xyt = xyt(withinRange, :);

timepnts = (xyt(:,3) - tbegin)/10^6;

xpos = xyt(:,1)*pix2cm;
ypos = xyt(:,2)*pix2cm;


%%%% calculation of speed

runWin = ones(15,1)/15;

diffx = [0; abs(diff(filtfilt(runWin,1,xpos)))];
diffy = [0; abs(diff(filtfilt(runWin,1,ypos)))];


difftt = [1; diff(timepnts)];

speed = sqrt(diffx.^2 + diffy.^2)./difftt; %% speed before smoothing


%%% defining the boundaries of the position bins

binEdges = min(xpos): posBinSize: max(xpos);

noPosBins = floor((max(xpos) - min(xpos))/posBinSize); %% posBinSize: the length of each position bin in cm
positionBins = binEdges(1:end-1) + posBinSize/2; %% center of the position bins


%%% to have a higher resolution in calculating position bin occupancy time

resolution = 0.001; %% in sec
timepnts_hr = 0.001:resolution:timepnts(end); %%   hr:high resolution
xpos_hr = interp1(timepnts, xpos, timepnts_hr);


%%%% spike position and convert to cm. Here I used only the X position of
%%%% the animal to linearize the animal's position. In other datasets we
%%%% may do a different conversion if animal's position in changing along
%%%% other axes

spikePos = spike.x * pix2cm;

%%% in case we are going to apply a restriction on the speed of the animal
%%% at the time of spiking

spkSpeed = interp1(timepnts, speed, spike.t/Fs);



%% Calcuation of the occupancy times and spike count of each unit per each position bin


totNoLaps = max(spike.lap); %%% total number of laps in the session
noUnits = length(unique(spike.unit));


lapDirection = zeros(1, totNoLaps);%%% initialization
   
occupTime = zeros(totNoLaps, noPosBins);

lapUnits = cell(1,totNoLaps);
noLapUnits = zeros(1,totNoLaps);

binFirRate = zeros(noUnits, noPosBins, totNoLaps);

for lap = 1: totNoLaps
    
    %%% direction and boundaries of the lap (Here we used spike.lap to
    %%% to define the boundaries)
    
    begSpikeInd = find(spike.lap == lap , 1, 'first');
    endSpikeInd = find(spike.lap == lap , 1, 'last');
    
    begTime = spike.t(begSpikeInd)/Fs; 
    endTime = spike.t(endSpikeInd)/Fs;
    
    
    if spikePos(begSpikeInd) > spikePos(endSpikeInd) 
       lapDirection(lap) = 1; %% right to left (leftward) (RL = higher X position to a lower X position)
    elseif spikePos(begSpikeInd) < spikePos(endSpikeInd) 
       lapDirection(lap) = 2; %% left to right (rightward) (LR)
    end
    

    %%% Calculation of the occupancy time of each position bin
    lapXpos = xpos_hr(find(timepnts_hr >= begTime & timepnts_hr <= endTime)); 

    for bin = 1 : noPosBins %%% calculating the occupancy duration for each position bin

        occupiedPosInd = find(lapXpos > binEdges(bin) & lapXpos <= binEdges(bin+1)); %% calculation of the number of position samples within the bin's position boundaries

        if isempty(occupiedPosInd)
            occupTime(lap,bin) = 1e3; %%% if the animal was not on specific position bin (probably if the animal do not enter some positions bins in either ends of the track)
        else
            occupTime(lap,bin) = length(occupiedPosInd) * resolution;
        end

    end

    %%% spike count per bin
    
    lapSpikeInd = find(spike.lap == lap & spkSpeed > speedThresh & ismember(spike.qclu , qclus));

    lapSpikeUnit = spike.unit(lapSpikeInd); 

    lapSpikeX = spikePos(lapSpikeInd);

    lapUnits{lap} = unique(lapSpikeUnit);
    noLapUnits(lap) = length(lapUnits{lap});

    for unit = 1 : noLapUnits(lap)
        
        unitSpkX = lapSpikeX(find(lapSpikeUnit == lapUnits{lap}(unit)));
        temp = histc(unitSpkX, binEdges);

        binFiring = temp(1: end-1); %%% the very end element is the number of spikes falls on the last edge, allways zero
        
        if size(binFiring, 1) > size(binFiring, 2)
            binFiring = binFiring';
        end
        
        binFirRate(lapUnits{lap}(unit), :, lap) = binFiring ./ occupTime(lap,:);

    end   

end

    
% alternating = unique(abs(diff(lapDirection))); %% if this is equal one then the data is alternating between LR and RL laps


directionName = {'left';'right'};


allActUnits = unique(cell2mat(lapUnits')); %%% set of firing units during the task



%%%% Smoothed firing rates

%%% smoothing kernel

sigma = 2; %%% 5cm
halfwidth = 3 * sigma; 
win = gausswindow(sigma, halfwidth);

binFirRate_sm = zeros(size(binFirRate)); 

for lap = 1 : totNoLaps
    for unit = 1 : noUnits
        temp = conv(binFirRate(unit,:,lap), win);
        binFirRate_sm(unit, :, lap) = temp(halfwidth+1 : end - halfwidth); 
    end
end




%% calculation of tuning and template for each trajectory direction

for direction = 1:2
    
    directionLaps = find(lapDirection == direction); %%% considering all the laps that have the same direction
%     dirTotNoLaps = length(directionLaps); %%% total number of laps with the desired direction
    
    
    %%% Calculating the mean over the laps
    
    MeanoverLaps = mean(binFirRate(:, :, directionLaps), 3); 
    %%% this is the unsmoothed tuning of the units

    
    %%% calculate the mean of smoothed firing rates over the laps
    MeanoverLaps_sm = mean(binFirRate_sm(:, :, directionLaps), 3);
    
    
    %%%% to check whether each unit is firing steadily within each bin
    
    firingThresh_low = 2; %% in Hz, to see if the distribution of the firing rate of each unit across laps passes a threshold significantly
    firingThresh_high = 5; %% threshold on peak firing to determine whether the unit is representing a place cell
    
    considerOnlySigBins = 0;
    
    binSigMat = zeros(noUnits, noPosBins);
    for bin = 1 : noPosBins
        for unit = 1 : noUnits
            binSigMat(unit, bin) = ttest(binFirRate_sm(unit, bin, :), firingThresh_low, 0.05,'right');
        end
    end
    
    
    %% Calculation of the trajectory template
    
    
    %%% Calculate the peak firing bin for each unit 


    peakRateBin = zeros(noUnits,1);
    isPlaceCell = zeros(noUnits,1); %% an index of whether the peak firing pas the higher threshold
    
    
    for ii = 1 : length(allActUnits) %%% we are just considering the active units here, not all the units
        
        unit = allActUnits(ii);
        sigbins = find(binSigMat(unit, :)); %%% bins with the unit's firing significantly above the threshold
        
%         if isempty(sigbins)
%            isPlaceCell(unit) = 0;
%             continue
%         end

        if considerOnlySigBins == 1

            [peakRate, maxBin] = max(MeanoverLaps_sm(unit, sigbins));
            peakRateBin(unit) = sigbins(maxBin);
            
        elseif considerOnlySigBins == 0

            [peakRate, maxBin] = max(MeanoverLaps_sm(unit, :));
            peakRateBin(unit) = maxBin;

        end
        
        if peakRate > firingThresh_high
            isPlaceCell(unit) = 1;
        end 

    end

    
    plcUnits = find(isPlaceCell);
    
    [sortedMaxBins, sortInd] = sort(peakRateBin(plcUnits), 'descend');
    sortedUnits = plcUnits(sortInd); %%% sorted units based on the corresponding peak rate bins
    
    
    sortedTunings = MeanoverLaps(sortedUnits, :);
    sortedTunings_sm = MeanoverLaps_sm(sortedUnits, :); %% used for plotting

    trajTemplate = sortedUnits;
    
    placeinfo = [sortedUnits (min(xpos)+sortedMaxBins*posBinSize)]

    
    currDir = pwd;
    FolderBase = [currDir '/' fileinfo.name '/data part' num2str(whichPart) '/PlaceFields/' directionName{direction} 'ward'];
    mkdir(FolderBase)
    
    save([FolderBase '/' 'PeakRateSortedTuning_' directionName{direction} 'ward.mat'], 'sortedTunings');
    save([FolderBase '/' 'AllCellsTuning_' directionName{direction} 'ward.mat'], 'MeanoverLaps');
    save([FolderBase '/' 'sortedUnits_' directionName{direction} 'ward.mat'], 'trajTemplate');
    save([FolderBase '/' 'activeUnits.mat'], 'allActUnits');
    save([FolderBase '/' 'PositionBinCenters.mat'], 'positionBins');

    
    %%%% ploting the sorted Tunings
    
    figure('Visible','off');
    for ii = 1 : length(sortedUnits)
        
        unit = sortedUnits(ii);
        
        subplot(length(sortedUnits), 1, ii) 

        area(positionBins, sortedTunings_sm(ii,:),'linewidth',1)

        if ii == 1
            title(['Place Fields for ' directionName{direction} 'ward trajectories'], 'fontsize',20)
        end
        
        ylabel(num2str(unit))
        
        set(gca, 'ytick', [])

        if ii < length(sortedUnits)
            set(gca,'xtick',[])
        end

    end
    
    saveas(gcf, [FolderBase '/' directionName{direction} 'ward_PlaceFields.fig'])
    
    
    %%% making clu and res files 
    
    rest = [];
    cluorder = [];
    shank = zeros(length(sortedUnits), 1);
    cluster = zeros(length(sortedUnits), 1);
    
    for ii = 1 : length(sortedUnits)
        spkInd = find(spike.unit == sortedUnits(ii));
        
        currSpiketimes = spike.t(spkInd);
        rest = [rest; currSpiketimes];
        cluorder = [cluorder; ii*ones(length(currSpiketimes),1)];
        
        shank(ii) = unique(spike.shank(spkInd));
        cluster(ii) = unique(spike.cluster(spkInd));
    end

    [sortedRest, ind] = sort(rest);
    cluorder = cluorder(ind);
    
    
    numclu = [shank cluster];
    numclu = numclu';

    Saveres([FolderBase '/' directionName{direction} 'Active.10' num2str(direction) '.res.' num2str(direction)],sortedRest);
    SaveClu([FolderBase '/' directionName{direction} 'Active.10' num2str(direction) '.clu.' num2str(direction)],cluorder);
    dlmwrite([FolderBase '/' directionName{direction} 'Active.10' num2str(direction) '.numclu.' num2str(direction)],numclu);
end

end
    
