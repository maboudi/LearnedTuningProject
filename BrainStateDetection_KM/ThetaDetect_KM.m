function [thetaPeriods, sdat] = ThetaDetect_KM(fileBase, sampleRate, nCh, thetaChannels, period)


max_thresholdf = 2;

% filter configuration

forder = 1024;  % filter order has to be even. The longer the more selective, but the operation
% will be linearly slower to the filter order. 100-125 for 1.25Hz (500 for 5 KH
forder = ceil(forder/2)*2; % to make sure filter order is even

hTheta = 12; % bandpass filter range
lTheta = 6; %

thresholdf = 0  ; % power SD threshold for high theta detection

min_theta_period = 500 ; % minimum theta period, 500ms ~ 4 cycles
max_theta_period = inf;

min_itheta_period = 100; % minimum inter-theta period


firfiltb = fir1(forder,[lTheta/sampleRate*2, hTheta/sampleRate*2]); % calculate convolution func


% smoothing kernel

sigma = 0.2 * sampleRate; %% in lfp sampling period (equal to 500 ms)
halfwidth = 3 * sigma;
smoothwin = gausswindow(sigma, halfwidth);



% load the data (.eeg from the select channels)
filtered_data = readmulti([fileBase '.eeg'], nCh, thetaChannels); % load .eeg trace


if ~isempty(period)
    filtered_data = filtered_data(period(1)*sampleRate+1 : period(2)*sampleRate, :);
end


% filter the data in theta frequency band
filtered_data2 = Filter0(firfiltb, filtered_data); % the actual filtering



% calculate the theta power and smooth is using the 

% filtered_data2 = abs(Shilbert(filtered_data)); % to power trace



% filtered_data2 = filtered_data2.^2; 
filtered_data2 = abs(hilbert(filtered_data2)); % to power trace


% filtered_data2 = sum(filtered_data2, 2); % taking average over the channels

filtered_data2 = conv(sum(filtered_data2,2), smoothwin, 'same');



% calculate the z scored amplitude 
sdat = (filtered_data2 - repmat(mean(filtered_data2), [size(filtered_data2, 1), 1])) ./ repmat(std(filtered_data2), [size(filtered_data2, 1), 1]);


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


% (2) merging theta periods, if inter-theta interval is less than min_itheta_period

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




%% Add the speed and sdf to eeg files 


eegInfo = dir([fileBase '.eeg']);
noTimePnts = eegInfo.bytes/nCh/2;

t = (1 : 1 : noTimePnts)* 1/sampleRate;


% sdtAdd2eeg = interp1((1:noBins)*binDur, sdat, t);
sdtAdd2eeg = sdat;

sdtAdd2eeg = (2 * (sdtAdd2eeg - min(sdtAdd2eeg)) / (max(sdtAdd2eeg) - min(sdtAdd2eeg)) - 1)*12000;


inputFileHandle = fopen([fileBase,'.eeg']);
outputFileHandle=fopen([fileBase,'-2.eeg'],'w');


bufferSize=4096;

doneFrame=0;
while ~feof(inputFileHandle)
    
    data = fread(inputFileHandle,[nCh,bufferSize],'int16')';
    for frame=1:size(data,1)
        fwrite(outputFileHandle,[data(frame,:), sdtAdd2eeg(frame+doneFrame)]','int16');
    end
    doneFrame=doneFrame+frame;
    
end

fclose(inputFileHandle);
fclose(outputFileHandle);


end