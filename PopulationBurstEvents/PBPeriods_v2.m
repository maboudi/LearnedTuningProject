function [PBEs, sdat] = PBPeriods_v2(spike, fileinfo, whichPart, binDur, ripplePeriods, exclude, VelocityFilter)

% Function PBPeriods characterizes periods which then will be used for determining
% population burst events (PBEs). To do this, first, spike density function
% (sdf or called spike histogram) is calculated by counting the multiunit spikes from all units (pyramidal
% and interneurons). Second, periods with sdf above x standard deviation
% above the mean in the whole period (with excluding theta or run periods)
% are characterized (pramary periods). Finally, the boundaries of the
% events will be the mean-crossing bins in either sides of the priamry
% periods (more specifically, the peaks within the primary periods are calculated
% and the mean-crossing points are dearched for within a window around the peak)

% the nonThetaPBEs will be in the format of [beginnig end center peakrate position speed]


% dataset info

% nCh = fileinfo.nCh;
tbegin = fileinfo.tbegin;
tend   = fileinfo.tend;

% Let's change the zero time and convert anything to second to make the next calculations easier
spike.t = spike.t(spike.t > tbegin & spike.t < tend);
spikeTimes = (spike.t - tbegin)/fileinfo.Fs ; %% (find(ismember(spike.qclu , qclus))) %% in case we want to see the population burst (synchrony) among just pyramidal units


%% Position Data processing

% in case we are going to do any filteration based on the speed or we need
% to find the locations where the PBEs are happening


% Loading POSITION data 

xyt = fileinfo.xyt;

if ~isfield(fileinfo, 'pix2cm')
    pix2cm = 1/3;
else
    pix2cm = fileinfo.pix2cm;
end


withinRange = find(xyt(:,3) >= tbegin & xyt(:,3) <= tend); % selecting only the position samples occured during neural recording
xyt = xyt(withinRange, :);

timepnts = (xyt(:,3) - tbegin)/fileinfo.Fs; % u seconds to seconds

xpos = xyt(:,1) * pix2cm; % Converting the position data (from pixel number to centimeter)
ypos = xyt(:,2) * pix2cm;


%%%% SPEED %%%%%%%%%%%%

diffx = [0; abs(diff(xpos))];
diffy = [0; abs(diff(ypos))];

difftt = [1; diff(timepnts)];

% calculation of speed in just x direction (we then can use the speed or theta for filtering the events)

% speed = abs(diffx)./difftt; %% speed before smoothing
speed = sqrt(diffx.^2 + diffy.^2)./difftt; 

% smoothing the speed

sigma = 15; %%% smoothing the speed, the positon sampling rate is around 30 Hz, so duration of the sigma is about 500 ms
halfwidth = 3 * sigma;

smoothwin = gausswindow(sigma, halfwidth);
speed = conv(speed, smoothwin, 'same');



%% Calculation of Spike Density function 

% Time binnig of the spike train, bin size = 1 ms and calculating the
% firing rate

noBins = floor(spikeTimes(end) / binDur);
binEdges = 0:binDur:spikeTimes(end);

[spikeBinned, ~] = histc(spikeTimes, binEdges);
spikeBinned(end) = []; %% removing the spikes that matches the last edge (we need the counts between the edges)

spikeDensity = spikeBinned / binDur;

% Final spike density function after smoothing the binned spikes

sigma = 20; %% in milisecond since we are using 1 milisecond time bins
halfwidth = 3 * sigma;

smoothwin = gausswindow(sigma, halfwidth);
spikeDensity = conv(spikeDensity, smoothwin, 'same');



% EXCLUDE 


ii_ok = ones(length(spikeDensity),1);


% exclude = floor(exclude .* 1000/lfpSampleRate);

% for analyzing post sleep

if ~isempty(exclude)

exclude = exclude((exclude(:, 2) > fileinfo.tbegin & exclude(:, 1) < fileinfo.tend), :);

end



if ~isempty(exclude)
    
% if exclude(1,1) <= tbegin
%     exclude(1,1) = tbegin + 0.001 * fileinfo.Fs; % 1 ms after the zero time (fileinfo.tbegin)
% end
% 
% if exclude(end, 2) > tend
%     exclude(end, 2) = tend;
% end


exclude(find(exclude(:,1) <= tbegin), 1) = tbegin + 0.001 * fileinfo.Fs;
exclude(find(exclude(:, 2) > tend), 2) = tend;


exclude = (exclude - tbegin)/fileinfo.lfpSampleRate; % Check whether we need to change this for Hiro's data


exclude = floor(exclude *  1000); % in ms

for ii = 1:size(exclude, 1)

    try
        ii_ok(exclude(ii,1):exclude(ii,2)) = 0;
    catch
        exclude(ii,1)
        exclude(ii,2)
    end

end

end

%%%% Mean and Standard deviation of the Spike Density in the non_theta periods

rateMean = mean(spikeDensity(ii_ok>0)); 
rateSD = std(spikeDensity(ii_ok>0)); 


sdat = zeros(length(spikeDensity), 1); 
sdat(ii_ok>0) = (spikeDensity(ii_ok>0) - rateMean)/rateSD; % note: the sdat is zero over the excluded theta periods


% PRIMARY PERIODS (periods with MUA above x SD)

%%% My strategy here was to find the periods with rate above the threshold
%%% and then find the maximum rate within those periods as the center for
%%% the rest of the analysis (e.g., defining the boundaries and so on)

% 
% % The epochs
% 
% highRateBins = find(sdat > threshZ); 
% 
% primPeriodStr = [highRateBins(1); highRateBins([0; diff(highRateBins)] > 1)];
% primPeriodEnd = [highRateBins(find([0; diff(highRateBins)] > 1)-1); highRateBins(end)];
% 
% 
% % Find the CENTER (bin with peak rate) in each of the epoches
% 
% noPrimPeriods = length(primPeriodStr);
% primPeriodCntr = zeros(noPrimPeriods, 1); % the center or bin with peak rate
% 
% primpeak = zeros(noPrimPeriods, 1);
% 
% for ii = 1 : noPrimPeriods
%     
%     [primpeak(ii), Ind] = max(sdat(primPeriodStr(ii) : primPeriodEnd(ii))); %% Ind signifiy the center bin in respect to the current period
%     primPeriodCntr(ii) = Ind + primPeriodStr(ii) - 1; 
% end
% 
% primPeriods = [primPeriodStr primPeriodEnd primPeriodCntr primpeak];

ripplePeriods(:,1:3) = round((ripplePeriods(:, 1:3) - tbegin)*1e3);
primPeriods = ripplePeriods;
noPrimPeriods = length(primPeriods);

% SECONDARY PERIODS (finging the edges around the peaks)

% Having determined the peaks, calculate the epochs' RISING and FALLING edges 
% first mean-passing point or local minimum afterward in the sides of the peak
% will be used as the edges. If niether found in either side within a
% search window, that epoch will be ignored.


halfSearchLen = 1000; % in ms the search period from either sides of the peaks for mean-crossing point, if any.

% Initializing

risingEdge = zeros(noPrimPeriods, 1);
fallingEdge = zeros(noPrimPeriods, 1);
eventCntr = zeros(noPrimPeriods, 1);
eventPeakFiring = zeros(noPrimPeriods, 1);
eventpeakRipplePower = zeros(noPrimPeriods, 1);

% Lets do the calculations for all the PRIMARY periods and correspnding
% peaks

epochInd = 0; % indices of the epochs with side mean-passing points within search windows of 1000 ms toward each side of the peak 

for epoch = 1 : noPrimPeriods
    
    
    peakTime = primPeriods(epoch, 3); % peak time of the current event
    
    startsearchPnt = max(1, peakTime-halfSearchLen); 
    endsearchPnt = min(noBins, peakTime + halfSearchLen);
    
    % Find the mean-passing points prior to the peak (RISING EDGE)

    meanPassPnt_prior = find(sdat(startsearchPnt : peakTime) < 0, 1, 'last');
    
    
    % Find the mean-passing points post to the peak (FALLING EDGE)
    
    meanPassPnt_post = find(sdat(peakTime : endsearchPnt) < 0, 1, 'first');
        
    if isempty(meanPassPnt_prior) || isempty(meanPassPnt_post) || peakTime < halfSearchLen

       continue
    else
        
       epochInd = epochInd + 1;
       
       risingEdge(epochInd) = peakTime - (halfSearchLen-meanPassPnt_prior)-1;
       fallingEdge(epochInd) = peakTime + meanPassPnt_post - 1;
       eventCntr(epochInd) = peakTime;
       
       % % modified from last version
       eventPeakFiring(epochInd) = max(sdat(risingEdge(epochInd):fallingEdge(epochInd)));
       eventpeakRipplePower(epochInd) = primPeriods(epoch, 4);
       % %
    end
    
end


% Create the secondary periods

secondaryPeriods = [risingEdge fallingEdge eventCntr eventPeakFiring eventpeakRipplePower]; %% the third and fourth columns are the center and the peak
secondaryPeriods(risingEdge == 0, :) = []; %% removing the extra zeros


% Detect overlapping epochs, if there is any

risefall_temp = secondaryPeriods(:, 1:2)'; % An array made of all rise and fall edges 
edges = risefall_temp(:); 
repeated_edges = find([edges(1); diff(edges)] < 0);

secondaryPeriods(floor(repeated_edges/2)+1 , :) = [];


% Calcualte the events' speed and position and add them to the 

if VelocityFilter
    
eventSpeed = interp1(timepnts, speed, secondaryPeriods(:, 3)/1e3); % calculating animal's speed at the peak of events
% % eventPos = interp1(timepnts, xpos, PBEs(:,3) / 1e3);

secondaryPeriods(eventSpeed > 5, :) = [];

end


% converting from milisecond resolution back to resolution of sampling frequency

% PBEs = (secondaryPeriods(:, 1:3) - 1) * binDur * fileinfo.Fs + fileinfo.tbegin + 1; %%% the fourth column is just the peak rate
PBEs = (secondaryPeriods(:, 1:3) - 1) * binDur * fileinfo.Fs + fileinfo.tbegin; %%% the fourth column is just the peak rate


PBEs(:,4:5) = secondaryPeriods(:, 4:5); % the peak population firing



% Generate the OUTPUT files realted to the nonThetaPBEs. Although these are not the final PBEs,
% visualizing the result before any further filteration could be needed as
% well.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% last commented
% 
% currDir = pwd;
% % Folderbase = [currDir '/' fileinfo.name '/data part' num2str(whichPart) '/Population Burst Events'];
% Folderbase = [currDir '/' fileinfo.name '/PopulationBurstEvents'];
% mkdir(Folderbase)
% 
% % los = {'long'; 'short'};
% 
% % save([currDir '/' fileinfo.name '-' fileinfo.animal '-' los{whichPart} '-' 'sdfEvents.mat'], 'PBEs') % The name signifies that the events are from processing sdf file without further binning and filtering, etc.
% save([currDir '/' fileinfo.name  '-' 'sdfEvents.mat'], 'PBEs') 
% 
% 
% % Filename = [currDir '/' fileinfo.name '-' fileinfo.animal '-' los{whichPart} '.sdf.evt'];
% Filename = [currDir '/' fileinfo.name '.sdf.evt'];
% 
% MakeEvtFile(PBEs(:, 1:3), Filename, {'beg', 'end', 'peak'}, fileinfo.Fs, 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% 
% % Add the speed and sdf to eeg files 
% 
% FileBase = [currDir '/' fileinfo.name '/' fileinfo.name];
% eegInfo = dir([FileBase '.eeg']);
% noTimePnts = eegInfo.bytes/nCh/2;
% 
% t = (1 : 1 : noTimePnts)* 1/lfpSampleRate;
% 
% 
% speedAdd2eeg = interp1(timepnts, speed, t);
% 
% sdfAdd2eeg = interp1((1:noBins)*binDur, sdat, t);
% 
% 
% speedAdd2eeg = (2 * (speedAdd2eeg - min(speedAdd2eeg)) / (max(speedAdd2eeg) - min(speedAdd2eeg)) - 1)*12000;
% 
% sdfAdd2eeg = (2 * (sdfAdd2eeg - min(sdfAdd2eeg)) / (max(sdfAdd2eeg) - min(sdfAdd2eeg)) - 1)*12000;
% 
% 
% inputFileHandle = fopen([FileBase,'.eeg']);
% outputFileHandle=fopen([FileBase,'-2.eeg'],'w');
% 
% 
% bufferSize=4096;
% 
% doneFrame=0;
% while ~feof(inputFileHandle)
%     
%     data = fread(inputFileHandle,[nCh,bufferSize],'int16')';
%     for frame=1:size(data,1)
%         fwrite(outputFileHandle,[data(frame,:),sdfAdd2eeg(frame+doneFrame),speedAdd2eeg(frame+doneFrame)]','int16');
%     end
%     doneFrame=doneFrame+frame;
%     
% end
% 
% fclose(inputFileHandle);
% fclose(outputFileHandle);


end

