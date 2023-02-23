function [PBEs, sdat] = PBPeriods_MUA_PYR2(MUA, PYR, detectParams, fileInfo)


% Function PBPeriods characterizes periods which later will be used for determining
% population burst events (PBEs). To do this, first, spike density function
% (sdf or called spike histogram) is calculated by counting the multiunit spikes from all units (pyramidal
% and maybe interneurons). Second, periods with sdf above x standard deviation
% above the mean in the whole period (with excluding theta or run periods)
% are characterized (pramary periods). Finally, the boundaries of the
% events will be the mean/zero-crossing points in either sides of the priamry
% periods (more specifically, the peaks within the primary periods are calculated
% and the mean-crossing points are dearched within a window around the peak)


binDur  = detectParams.time_resolution;
threshZ = detectParams.threshold;
exclude = detectParams.exclude;

spike_set = {'MUA'; 'PYR'};


% dataset info
tbegin = fileInfo.tbegin;
tend   = fileInfo.tend;


% Let's change the zero time and convert anything to second to make the next calculations easier

spikeTimes.MUA = MUA.time;
spikeTimes.MUA = spikeTimes.MUA(spikeTimes.MUA > tbegin & spikeTimes.MUA < tend);
spikeTimes.MUA = (spikeTimes.MUA - tbegin)/fileInfo.timeUnit;

spikeTimes.PYR = PYR.time;
spikeTimes.PYR = spikeTimes.PYR(spikeTimes.PYR > tbegin & spikeTimes.PYR < tend);
spikeTimes.PYR = (spikeTimes.PYR - tbegin)/fileInfo.timeUnit;


%% Calculation of Spike Density function 

% Time binnig of the spike train, bin size = 1 ms and calculating the
% firing rate

nBins = floor(spikeTimes.MUA(end)/binDur);
binEdges = 0:binDur:spikeTimes.MUA(end);


% Exclude the REM/Active wake periods if we want to limit the z-scoring to
% the NREM/quite wake periods

ii_ok = ones(nBins,1);

if ~isempty(exclude)
    exclude = exclude((exclude(:, 2) > fileInfo.tbegin & exclude(:, 1) < fileInfo.tend), :);
end


if ~isempty(exclude)

    exclude(exclude(:,1) <= tbegin, 1) = tbegin + 0.001 * fileInfo.Fs;
    exclude(exclude(:, 2) > tend, 2)   = tend;

    exclude = (exclude - tbegin)/fileInfo.timeUnit; 

    exclude = floor(exclude*1e3); % in ms

    for ii = 1:size(exclude, 1)

        try
            ii_ok(exclude(ii,1):exclude(ii,2)) = 0;
        catch
            exclude(ii,1)
            exclude(ii,2)
        end
    end
end



for is = 1:2

    curr_spike_set = spike_set{is};

    spikeBinned.(curr_spike_set) = histc(spikeTimes.(curr_spike_set), binEdges);
    spikeBinned.(curr_spike_set)(end) = []; %% removing the spikes that matches the last edge (we need the counts between the edges)

    spikeDensity.(curr_spike_set) = spikeBinned.(curr_spike_set)/binDur;

    % Final spike density function after smoothing the binned spikes
    
    
    sigma     = detectParams.smoothingSigma.(curr_spike_set)*1e3; %% in milisecond since we are using 1 milisecond time bins
    halfwidth = 3 * sigma;
    
    smoothwin    = gausswindow(sigma, halfwidth);


    spikeDensity.(curr_spike_set) = conv(spikeDensity.(curr_spike_set), smoothwin, 'same');



    %%%% Mean and Standard deviation of the Spike Density in the non_theta periods

    rateMean = mean(spikeDensity.(curr_spike_set)(ii_ok>0)); 
    rateSD = std(spikeDensity.(curr_spike_set)(ii_ok>0)); 

    sdat.(curr_spike_set) = (spikeDensity.(curr_spike_set) - rateMean)/rateSD;

    
    %% PRIMARY PERIODS (periods with MUA above x SD)

    % only MUA spikes are used to detect peak of population bursts

    % The epochs

    highRateBins = find(sdat.(curr_spike_set) > threshZ.(curr_spike_set)); 

    primPeriodStr = [highRateBins(1); highRateBins([0; diff(highRateBins)] > 1)];
    primPeriodEnd = [highRateBins(find([0; diff(highRateBins)] > 1)-1); highRateBins(end)];


    % Find the CENTER (bin with peak rate) in each of the epoches

    nPrimPeriods = length(primPeriodStr);
    primPeriodCntr = zeros(nPrimPeriods, 1); % the center or bin with peak rate

    primpeak = zeros(nPrimPeriods, 1);

    for ii = 1 : nPrimPeriods

        [primpeak(ii), Ind] = max(sdat.(curr_spike_set)(primPeriodStr(ii) : primPeriodEnd(ii))); %% Ind signifiy the center bin in respect to the current period
        primPeriodCntr(ii) = Ind + primPeriodStr(ii) - 1; 
    end

    primPeriods = [primPeriodStr primPeriodEnd primPeriodCntr primpeak];



    %% SECONDARY PERIODS (detect edges around the peak MUAs)

    % Having determined the peaks, calculate the events RISING and FALLING edges 
    % as first mean-passing points immediately before and following the peak. If niether found in either side within a
    % search window, that event will be omitted.


    halfSearchLen = 1000; % in ms the search period from either sides of the peaks for mean-crossing point.


    curr_spike_set = spike_set{is};

    % Initializing

    risingEdge      = zeros(nPrimPeriods, 1);
    fallingEdge     = zeros(nPrimPeriods, 1);
    eventCntr       = zeros(nPrimPeriods, 1);
    eventPeakFiring = zeros(nPrimPeriods, 1);



    eventNum = 0; 
    for epoch = 1 : nPrimPeriods


        peakTime = primPeriods(epoch, 3); % peak time of the current event

        startsearchPnt = max(1, peakTime-halfSearchLen); 
        endsearchPnt   = min(nBins, peakTime+halfSearchLen);



        % Find the mean-passing points prior to the peak (RISING EDGE)

        zeroCrossing_prior = find(sdat.(curr_spike_set)(startsearchPnt:peakTime) < 0, 1, 'last');



        % Find the mean-passing points post to the peak (FALLING EDGE)

        zeroCrossing_post = find(sdat.(curr_spike_set)(peakTime : endsearchPnt) < 0, 1, 'first');




        if isempty(zeroCrossing_prior) || isempty(zeroCrossing_post) || startsearchPnt == 1 || endsearchPnt == nBins
           continue
        else

           eventNum = eventNum + 1;

           risingEdge(eventNum)  = peakTime - (halfSearchLen-zeroCrossing_prior)-1;
           fallingEdge(eventNum) = peakTime + zeroCrossing_post - 1;

           eventCntr(eventNum) = peakTime;
           eventPeakFiring(eventNum) = primPeriods(epoch, 4);

        end

    end


    % Create the secondary periods

    secondaryPeriods = [risingEdge fallingEdge eventCntr eventPeakFiring]; %% the third and fourth columns are the center and the peak
    secondaryPeriods(risingEdge == 0, :) = []; %% removing the extra zeros



    % get rid of repeated PBEs
    min_distance = 0;
    tmp_pbe = secondaryPeriods(1,:);
    third = [];

    for ii = 2:size(secondaryPeriods,1)
        if (secondaryPeriods(ii,1)-tmp_pbe(2)) < min_distance

            if secondaryPeriods(ii,4) > tmp_pbe(4)
                currpeakT   = secondaryPeriods(ii, 3);
                currpeakMUA = secondaryPeriods(ii, 4);
            else
                currpeakT   = tmp_pbe(3);
                currpeakMUA = tmp_pbe(4);
            end

            tmp_pbe = [tmp_pbe(1) max(tmp_pbe(2), secondaryPeriods(ii,2)) currpeakT currpeakMUA]; % merge two PBEs
        else
            third = [third; tmp_pbe];
            tmp_pbe = secondaryPeriods(ii,:);
        end
    end

    third = [third;tmp_pbe]; 

    PBEs.(curr_spike_set)      = third(:, 1:3)*binDur + fileInfo.tbegin; %%% the fourth column is just the peak rate
    PBEs.(curr_spike_set)(:,4) = third(:, 4);

    
    duration = PBEs.(curr_spike_set)(:, 2) - PBEs.(curr_spike_set)(:, 1);
    PBEs.(curr_spike_set) = PBEs.(curr_spike_set)(duration >= detectParams.minDur & duration <= detectParams.maxDur, :);


end


% Accepting PYR-based detected PBEs that overlapped with the MUA-based
% ones (the rest are most probably false positives)

n_PYR_PBEs = size(PBEs.PYR, 1);
PYR_PBE_OverlappingIDX = false(n_PYR_PBEs, 1);

for ipbe = 1:n_PYR_PBEs
    
   idx = find(PBEs.MUA(:, 1) <= PBEs.PYR(ipbe, 3) & PBEs.MUA(:, 2) >= PBEs.PYR(ipbe, 3), 1, 'first'); 
   
   if ~isempty(idx)
       PYR_PBE_OverlappingIDX(ipbe) = true;
   end
   
end

PBEs.PYR = PBEs.PYR(PYR_PBE_OverlappingIDX, :);



end

