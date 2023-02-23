function thetaPeriods = ThetaRatio_KM(fileBase, sampleRate, nCh, thetaChannels, period)


max_thresholdf = 3;

% filter configuration

forder = 256;  % filter order has to be even. The longer the more selective, but the operation
% will be linearly slower to the filter order. 100-125 for 1.25Hz (500 for 5 KH
forder = ceil(forder/2)*2; % to make sure filter order is even

hTheta = 12; % bandpass filter range
lTheta = 6; %


hDelta = 5;
lDelta = 1;

thresholdf = 0  ; % power SD threshold for high theta detection

min_theta_period = 500 ; % minimum theta period, 50ms ~ 4 cycles
max_theta_period = inf;

min_itheta_period = 100; % minimum inter-theta period


firfiltb_th = fir1(forder,[lTheta/sampleRate*2, hTheta/sampleRate*2]); % calculate convolution func
firfiltb_d = fir1(forder,[lDelta/sampleRate*2, hDelta/sampleRate*2]); % calculate convolution func

% smoothing kernel

sigma = 0.5 * sampleRate; %% in lfp sampling period (equal to 500 ms)
halfwidth = 3 * sigma;
smoothwin = gausswindow(sigma, halfwidth);



% load the data (.eeg from the select channels)
data = readmulti([fileBase '.eeg'], nCh, thetaChannels); % load .eeg trace


if ~isempty(period)
    data = data(period(1)*sampleRate+1 : period(2)*sampleRate, :);
end


% filter the data in theta and delta frequency bands

filtered_data_th = Filter0(firfiltb_th, data); % the actual filtering
filtered_data_d = Filter0(firfiltb_d, data); % the actual filtering



% calculate the theta power and smooth is using the 

% filtered_data2 = abs(Shilbert(filtered_data)); % to power trace


% calculating amplitude in each of the frequency bands and smmothing
% afterwards

% theta

filtered_data2_th = filtered_data_th.^2;
filtered_data2_th = sum(filtered_data2_th, 2); % taking average over the channels
filtered_data2_th = conv(sum(filtered_data2_th,2), smoothwin, 'same');


% delta

filtered_data2_d = filtered_data_d.^2;
filtered_data2_d = sum(filtered_data2_d, 2); % taking average over the channels
filtered_data2_d = conv(sum(filtered_data2_d,2), smoothwin, 'same');


% calcualte the theta to delta ratio
thetaDeltaRatio = filtered_data2_th./filtered_data2_d;


% calculate the z scored amplitude 
sdat = (thetaDeltaRatio - repmat(mean(thetaDeltaRatio), [size(thetaDeltaRatio, 1), 1])) ./ repmat(std(thetaDeltaRatio), [size(thetaDeltaRatio, 1), 1]);



% [sdat, stdA] = unity(conv(sum(filtered_data2,2), smoothwin, 'same'),exclude); % averaging & standardizing
% sdat = unity(sum(filtered_data2,2),exclude); % averaging & standardizing


% (1) primary detection of theta periods, based on thresholding


thresholded = sdat > thresholdf;
primary_detection_b = find(diff(thresholded) > 0);
primary_detection_e = find(diff(thresholded) < 0);

if primary_detection_b(1) > primary_detection_e(1)
    primary_detection_e = primary_detection_e(2:end);
end

sum_detection =[primary_detection_b(1:size(primary_detection_e,1),1) primary_detection_e];
% sum_detection_e = [sum_detection_e primary_detection_e];

sum_detection = sortrows(sum_detection);
primary_detection_b = sum_detection(:,1);
primary_detection_e = sum_detection(:,2);

if (length(primary_detection_e) == length(primary_detection_b)-1) % exclude ranged-out (1st or
    primary_detection_b = primary_detection_b(1:end-1);
end

if (length(primary_detection_e)-1 == length(primary_detection_b))
    primary_detection_e = primary_detection_e(2:end);
end

primary = [primary_detection_b,primary_detection_e]; % [start, end]


% (2) merging ripples, if inter-ripples period is less than min_isw_period

min_itheta_period = min_itheta_period/1000*sampleRate; % in eeg
min_theta_period = min_theta_period/1000*sampleRate; % in eeg
max_theta_period = max_theta_period/1000*sampleRate;
secondary=[];
tmp_rip = primary(1,:);

for ii=2:size(primary,1)
    if (primary(ii,1)-tmp_rip(2)) < min_itheta_period
        tmp_rip = [tmp_rip(1),primary(ii,2)]; % merge two ripples
    elseif abs(primary(ii,1)-tmp_rip(1)) < min_itheta_period
        tmp_rip = [min([tmp_rip(1) primary(ii,1)]),max([tmp_rip(2) primary(ii,2)])]; % merge two ripples
    else
        secondary = [secondary;tmp_rip];
        tmp_rip = primary(ii,:);
    end
end

secondary = [secondary;tmp_rip]; % [start, end]

keeper = find((secondary(:,2)-secondary(:,1)) > min_theta_period & (secondary(:,2)-secondary(:,1)) < max_theta_period);
secondary = secondary(keeper,:);



% (3) Theta periods must have its peak power of > max_thresholdf

third = [];
% SDmax = [];

for ii=1:size(secondary,1)
    [max_val,max_idx] = max(sdat([secondary(ii,1):secondary(ii,2)]));
%     mean_val = mean(sdat([secondary(ii,1):secondary(ii,2)]));
    if max_val > max_thresholdf % & mean_val>thresholdf
        third = [third; secondary(ii,:)]; % [start, end]
%         SDmax = [SDmax; max_val];
    end
end


if ~isempty(period)
    thetaPeriods = third + period(1)*sampleRate;
else
    thetaPeriods = third;
end


MakeEvtFile(thetaPeriods,[fileBase '.tht.evt'], {'beg', 'end'}, sampleRate); % make evt file for neuroscope browsing


end