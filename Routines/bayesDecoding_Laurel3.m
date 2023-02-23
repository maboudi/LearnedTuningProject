clear; clc

cd('/home/kouroshmaboudi/Documents/LaurelData/ObjectLocationMemory/Charlie')

% Dataset directory
currDir = pwd;

%% session info

animal = 'charlie';
sessionName = 'PlacePreference-Charlie-20150623-01-nlx';

fileinfo = struct('name', sessionName, 'animal', animal, 'xyt', [], 'tbegin', [], 'tend', ...
                    [], 'Fs', [], 'nCh', [], 'lfpSampleRate', [], 'CA1thetach', [], 'CA', [1 1 1 1], 'pos2cm', 240); 


% (1)  Sleep 1: ~20 minute baseline sleep in homecage.
% (2)  Exposure 1: exposure to empty box (~6 min).
% (3)  Homecage 1: return to home cage (~3 min).
% (4)  Exposure 2: exposure to box with objects (~6 min).
% (5)  Homecage 2: return to homecage (~3 min). 
% (6)  Exposure 3: exposure to box with objects (~6 min).
% (7)  Homecage 3: return to homecage (~3 min). 
% (8)  Exposure 4: exposure to box with objects (~6 min). 
% (9)  Sleep 2: extended sleep (~4 hours). 
% (10) Memory recall test: exposure to box with displaced object (~6 min).
% (11) rest of the reocring for a short period (~20 min)


FileBase = [currDir '/' fileinfo.name '/' fileinfo.name];

% Loading the mat fiels
% varList = {'basics'; 'speed'; 'spikes'; 'hiAmpSpindles'};
varList = {'-basics'; '-speed'; '-spikes'; '.SleepState.states'};

for var = 1:length(varList)
    load([FileBase  varList{var} '.mat'])
end


% fileinfo.Fs = basics.SampleRate;

fileinfo.Fs = 1e6; % we don't use the original sample rate, because the time stamps were converted to micro sec. 
fileinfo.lfpSampleRate = basics.lfpSampleRate; 
fileinfo.nCh = basics.nChannels;

fileinfo.xyt = [speed.x*fileinfo.pos2cm speed.y*fileinfo.pos2cm speed.t]; 


% periods
periods = basics.period.time;


%% processing spikes

nUnits = length(spikes);

for unit = 1:length(spikes)
    
    spikes(unit).rate = length(spikes(unit).time)/((periods(end, 2) - periods(1,1))/1e6);

end



spike = struct('t', [], 'unit', [], 'shank', [], 'qclu', [], 'x', [], 'y', [], 'speed', [], 'rate', []);

for unit = 1: nUnits
    
    if spikes(unit).quality == 0
        continue
    end
    
    unitSpikes = spikes(unit).time;
    
    spike.t = [spike.t; unitSpikes];
    
    spike.unit = [spike.unit; repelem(unit, length(unitSpikes))'];
    spike.shank = [spike.shank; repelem(spikes(unit).id(1), length(unitSpikes))'];
    spike.qclu = [spike.qclu; repelem(spikes(unit).quality, length(unitSpikes))'];
    spike.rate = [spike.rate; repelem(spikes(unit).rate, length(unitSpikes))'];
    
    spike.x = [spike.x; spikes(unit).x*fileinfo.pos2cm];
    spike.y = [spike.y; spikes(unit).y*fileinfo.pos2cm];
    spike.speed = [spike.speed; spikes(unit).speed];
    
    
end


% sorting the spike time stamps
[~, spikeSortIdx] = sort(spike.t, 'ascend');

spike = structfun(@(x)(x(spikeSortIdx)), spike,'UniformOutput',false);


% spike.qclu(spike.rate > 2 & ismember(spike.qclu, [1 3])) = 2; % characterizing high firing units, to exclude them later




%% Active-run periods


velocity = speed.v; 


% smoothing the speed; the sampling rate is 100 Hz (after interpolation).

sigma = 100; % the sampling is every 0.01 sec, so the sigma is equal to 1 sec
halfwidth = 3*sigma;
win = gausswindow(sigma, halfwidth); % smoothing kernel

velocity = conv(velocity, win, 'same');  



% active run periods

velocityIdx = velocity > 2; % velocity above x cm/s

tran2HiVelPnts = find([0; diff(velocityIdx)] == 1);
tran2LoVelPnts = find([0; diff(velocityIdx)] == -1);



if length(tran2HiVelPnts) == length(tran2LoVelPnts)+1
    
    tran2LoVelPnts = [tran2LoVelPnts; find(velocityIdx, 1, 'last')];
    
elseif length(tran2LoVelPnts) == length(tran2HiVelPnts)+1

    tran2HiVelPnts = [find(velocityIdx, 1, 'first'); tran2HiVelPnts];
    
elseif tran2LoVelPnts(1) < tran2HiVelPnts(1)

    tran2LoVelPnts = [tran2LoVelPnts; find(velocityIdx, 1, 'last')];
    tran2HiVelPnts = [find(velocityIdx, 1, 'first'); tran2HiVelPnts];
    
end


activePeriods = [tran2HiVelPnts tran2LoVelPnts];
activePeriods = activePeriods/100*1e6 + periods(1,1);



% Active periods that happen during the exposure of the animal to the
% objects


boxActivePeriods = activePeriods(activePeriods(:,1) > periods(4,1) & activePeriods(:,2) < periods(8,2), :);


% boxActivePeriods = [periods(4,1) periods(8,2)];


%% place field calculation: 2D


qclus = [1 3];

posBinSize = 2; % in cm

[rateMaps, xposcenters, yposcenters, occupancy, spikeCount] = placeField_2D(spike, qclus, boxActivePeriods, posBinSize, fileinfo);



noUnits = size(rateMaps, 3);
noUnitsperRow = 5;


peakRate = zeros(noUnits, 1);


figure; 

for ii = 1:noUnits
    
    
subplot(ceil(noUnits/noUnitsperRow), noUnitsperRow, ii)


peakRate(ii) = max(max(rateMaps(:, :, ii)));


plot(fileinfo.xyt(:, 1), fileinfo.xyt(:, 2), '.', 'markersize', 1, 'MarkerEdgeColor', [0.7 0.7 0.7])


hold on

h1 = imagesc(xposcenters, yposcenters, spikeCount(:,:,ii)'); %
set(h1, 'AlphaData', 0.7)
text(min(xposcenters), max(yposcenters)+8, sprintf('%.1f', peakRate(ii)))

xlim([min(xposcenters) max(xposcenters)])
ylim([min(yposcenters) max(yposcenters)])


cornerPlot = (ceil(noUnits/noUnitsperRow)-1) * noUnitsperRow + 1;

if ii == cornerPlot
     
    set(gca, 'XTick', [min(xposcenters) max(xposcenters)], 'XTickLabel', {sprintf('%.1f', min(xposcenters)); sprintf('%.1f', max(xposcenters))}, ...
        'YTick', [min(yposcenters) max(yposcenters)], 'YTickLabel', {sprintf('%.1f', min(yposcenters)); sprintf('%.1f', max(yposcenters))})

else
    set(gca, 'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', [])
end


end



%% Functions

function [rateMaps, xposCenters, yposCenters, occupancy, spikeCount] = placeField_2D(spike, qclus, activePeriods, posBinSize, fileinfo)


minX = min(fileinfo.xyt(:,1));
minY = min(fileinfo.xyt(:,2));


noXBins = ceil(range(fileinfo.xyt(:,1))/posBinSize);
noYBins = ceil(range(fileinfo.xyt(:,2))/posBinSize);

xposCenters = minX + (1:noXBins)*posBinSize - posBinSize/2 ;
yposCenters = minY + (1:noYBins)*posBinSize - posBinSize/2 ;


% only the spikes from (low firing) pyramidal units which ocurrs during the active periods 

keepIdx = find(ismember(spike.qclu, qclus)); 
pyrSpikes = structfun(@(x)(x(keepIdx)), spike,'UniformOutput',false);

uniqUnits = unique(pyrSpikes.unit); % all the units


% smoothing kernel for spikeCount and occupancy maps 

sigma = 1;
halfwidth = 2*sigma;
smoothwin = gausswindow(sigma, halfwidth);
smoothwin2D = smoothwin' * smoothwin;


spikeCount = zeros(noXBins, noYBins, length(uniqUnits));
occupancy = zeros(noXBins, noYBins);

for ii = 1: size(activePeriods, 1)
    
    perSpikeCount = zeros(noXBins, noYBins, length(uniqUnits));
    
    for unit = 1: length(uniqUnits)
        
        unitID = uniqUnits(unit);
        
        unitSpkIdx = find(pyrSpikes.t > activePeriods(ii, 1) & pyrSpikes.t < activePeriods(ii,2) & pyrSpikes.unit == unitID);
        unitSpikes = structfun(@(x)(x(unitSpkIdx)), pyrSpikes, 'UniformOutput', false);

        perSpikeCount(:,:,unit) = full(sparse(ceil((unitSpikes.x - minX + eps)/posBinSize), ceil((unitSpikes.y - minY + eps)/posBinSize), 1, noXBins, noYBins));
               
    end
    
    
    periodPosIdx = find(fileinfo.xyt(:,3) > activePeriods(ii, 1) & fileinfo.xyt(:,3) < activePeriods(ii, 2));
    
    spikeCount = spikeCount + perSpikeCount;
    occupancy = occupancy + full(sparse(ceil((fileinfo.xyt(periodPosIdx,1) - minX + eps)/posBinSize), ceil((fileinfo.xyt(periodPosIdx,2) - minY + eps)/posBinSize), 1, noXBins, noYBins))*0.01 ; % in second
    
end


rateMaps = spikeCount ./ repmat(occupancy, [1,1, length(uniqUnits)]);

rateMaps(isnan(rateMaps)) = 0;
rateMaps(isnan(rateMaps)) = 0;


for ii = 1: length(uniqUnits)
    rateMaps(:,:, ii) = conv2(rateMaps(:,:, ii), smoothwin2D, 'same');   
end


for ii = 1: length(uniqUnits)
    spikeCount(:,:, ii) = conv2(spikeCount(:,:, ii), smoothwin2D, 'same');   
end

occupancy = conv2(occupancy, smoothwin2D, 'same');  


end

