function [PBEs, sdfZ] = DetectPopBursts(spike, binDur, threshZ, exclude, fileinfo)


% (1) spike density function(sdf) is calculated by counting the number of multiunit
%     spikes in one milisecond time bins.
% (2) Peaks of sdf passing a certain threshold (mean(sdf) + threshZ*std(sdf)) are chalculated. 
% (3) The boundary of each event is determined as the mean-crossing points at each side of the peak.


minInterEvent = 0;



% Loading dataset info

% nCh = fileinfo.nCh;
Fs = fileinfo.Fs;
tbegin = fileinfo.tbegin; 
tend = fileinfo.tend;


spike.t = spike.t(spike.t > tbegin & spike.t < tend); % including the spikes from all units
spikeTimes = (spike.t - tbegin)/Fs ; 


%% Spike Density function 


noBins = floor(spikeTimes(end)/binDur);
binEdges = 0:binDur:spikeTimes(end);

[spikeBinned, ~] = histc(spikeTimes, binEdges);
spikeBinned(end) = []; 
spikeDensity = spikeBinned/binDur; % spike rate at each time bin


% smoothing the spike rate

sigma = 20; %% in ms
halfwidth = 3 * sigma;

smoothwin = gausswindow(sigma, halfwidth);
spikeDensity = conv(spikeDensity, smoothwin, 'same');



%%%% Mean and Standard deviation of the Spike Density 

if ~isempty(exclude)
    exclude = exclude((exclude(:, 2) > fileinfo.tbegin & exclude(:, 1) < fileinfo.tend), :);
end


ii_ok = ones(length(spikeDensity),1);
if ~isempty(exclude)
    
    
    if exclude(1,1) <= tbegin
        exclude(1,1) = tbegin + 0.001 * Fs; 
    end

    if exclude(end, 2) > tend
        exclude(end, 2) = tend;
    end

    exclude = exclude - tbegin; 

    
    exclude = floor(exclude/Fs * 1000); % in ms
    
    
    for ii = 1:size(exclude, 1)

        ii_ok(exclude(ii,1):exclude(ii,2)) = 0;

    end

end



rateMean = mean(spikeDensity(ii_ok>0)); 
rateSD = std(spikeDensity(ii_ok>0)); 


sdfZ = zeros(length(spikeDensity), 1); 
sdfZ(ii_ok>0) = (spikeDensity(ii_ok>0) - rateMean)/rateSD; % the sdfZ is zero over the excluded periods




% PRIMARY PERIODS (periods with MUA above the threshold)


highRateBins = find(sdfZ > threshZ); 

primPeriodStr = [highRateBins(1); highRateBins([0; diff(highRateBins)] > 1)];
primPeriodEnd = [highRateBins(find([0; diff(highRateBins)] > 1)-1); highRateBins(end)];


% find the CENTER (bin with the peak rate) within each of the periods

noPrimPeriods = length(primPeriodStr);
primPeriodCntr = zeros(noPrimPeriods, 1); % the time bin with the peak rate

primPeakRate = zeros(noPrimPeriods, 1);

for ii = 1 : noPrimPeriods
    
    [primPeakRate(ii), idx] = max(sdfZ(primPeriodStr(ii) : primPeriodEnd(ii))); 
    primPeriodCntr(ii) = idx + primPeriodStr(ii) - 1; 
end

primPeriods = [primPeriodStr primPeriodEnd primPeriodCntr primPeakRate];



% SECONDARY PERIODS (calculating the boundaries of the periods)

halfSearchLen = 1000; % in ms the search period from either sides of the peaks for the mean-crossing points (if there was any)


% Initializing

risingEdge = zeros(noPrimPeriods, 1);
fallingEdge = zeros(noPrimPeriods, 1);
eventCntr = zeros(noPrimPeriods, 1);
eventPeakFiring = zeros(noPrimPeriods, 1);



secEvtIdx = 0;

for evt = 1 : noPrimPeriods
    
    
    peakTime = primPeriods(evt, 3); % peak time of the current event
    
    startsearchPnt = max(1, peakTime-halfSearchLen); 
    endsearchPnt = min(noBins, peakTime + halfSearchLen);
    
    % Find the mean-passing points prior to the peak (RISING EDGE)

    meanPassPnt_prior = find(sdfZ(startsearchPnt : peakTime) < 0, 1, 'last');
    
    
    % find the mean-passing points post to the peak (FALLING EDGE)
    
    meanPassPnt_post = find(sdfZ(peakTime : endsearchPnt) < 0, 1, 'first');
        
    if isempty(meanPassPnt_prior) || isempty(meanPassPnt_post) 

       continue
    else
        
       secEvtIdx = secEvtIdx + 1;
       
       risingEdge(secEvtIdx) = peakTime - (halfSearchLen-meanPassPnt_prior)-1;
       fallingEdge(secEvtIdx) = peakTime + meanPassPnt_post - 1;
       eventCntr(secEvtIdx) = peakTime;
       eventPeakFiring(secEvtIdx) = primPeriods(evt, 4);
       
    end
    
end


% Create the secondary periods

secondaryPeriods = [risingEdge fallingEdge eventCntr eventPeakFiring]; 
secondaryPeriods(risingEdge == 0, :) = [];


% remove the redundancies

risefall_temp = secondaryPeriods(:, 1:2)'; 
edges = risefall_temp(:); 
repeated_edges = find([edges(1); diff(edges)] < minInterEvent); % 0 can be replaced with a minimum inter-event interval, if we want to merge the events that are close together

secondaryPeriods(floor(repeated_edges/2)+1 , :) = [];


% converting from milisecond to the original sampling frames

PBEs = (secondaryPeriods(:, 1:3) - 1) * binDur * Fs + fileinfo.tbegin; % the fourth column is just the peak rate

PBEs(:,4) = secondaryPeriods(:, 4); 



% save the PBEs

currDir = pwd;
basePath = [currDir '/' fileinfo.name '/PopulationBurstEvents'];
mkdir(basePath)


save(fullfile(basePath, [fileinfo.name '-PBEs.mat']), 'PBEs') 

MakeEvtFile(PBEs(:, 1:3), fullfile(basePath, [fileinfo.name '.pbe.evt']), {'beg', 'end', 'peak'}, Fs, 1)


end

function win = gausswindow(sigma, halfwidth)

mu = 0;
x = mu - halfwidth : mu + halfwidth;

for i = 1 : length(x)
    y(i) = (1/sigma*sqrt(2*pi)) * exp(-(x(i) - mu)^2 / 2/sigma^2);
end

win = y./sum(y);

end
