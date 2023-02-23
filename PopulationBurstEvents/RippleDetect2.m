function [ripple_result, rippleLFP, sdat, sdat_summed] = RippleDetect2(FileBase, fileinfo, varargin)
% function [ripple_result] = RippleDetect(fileinfo,without,max_thresholdf,homefilter,runfilter,thetafilter)
%
% function is intended to detect ripple events, from 130-230 Hz, in the eeg file, given
% fileinfo, an array which contains the name and selected channels of the
% directory in question.

[max_thresholdf,homefilter,runfilter,thetafilter, period] = DefaultArgs(varargin,{5,0,0,0,[]});

% homefilter: whether or not the rat should be in the home cage aread
%  runfilter: whether or not there should be a limit on the rat's running speed
% thetafilter whether or not we require robust theta oscillations
% hi!slacker ;)
% max_thresholdf = 5; peak SD threshold

videorate = 120;

% currentdir = pwd;
% FileBase = [currentdir '/' fileinfo.name '/' fileinfo.name];

fprintf('%s  detecting ripples...',fileinfo.sessionName);

highband = 250; % bandpass filter range
lowband = 150; % (250Hz to 100Hz)

thresholdf = 2  ; % power SD threshold for ripple detection (data mean in case of zero)

% min_sw_period = 50 ; % minimum ripple period, 50ms ~ 5-10 cycles
min_sw_period = 15 ; % minimum ripple period, 50ms ~ 5-10 cycles\
max_sw_period = 300;
% min_isw_period = 100; % minimum inter-ripple period
min_isw_period = 0; % minimum inter-ripple period



%%%%% Configuration %%%%%

Par = LoadPar([FileBase '.xml']); % either par or xml file needed

numchannel = Par.nChannels; % total channel number, 6, 10 or 18
EegRate = Par.lfpSampleRate; % sampling rate of eeg file
% SampleRate = Par.SampleRate; % sampling rate of dat file


forder = 100;  % filter order has to be even. The longer the more selective, but the operation
% will be linearly slower to the filter order. 100-125 for 1.25Hz (500 for 5 KH)

avgfilorder = round(min_sw_period/1000*EegRate/2)*2+1 ; % should not change this. length of averaging filter
% avgfilorder = 101; % should not change this. length of averaging filter
forder = ceil(forder/2)*2; % to make sure filter order is even

%%%%% Ripple detection %%%%%
%       
% if isempty(fileinfo.ripch)
%     fprintf(1,'No CA1 ripple channels available\n')
%     ripple_result = [];
%     save([FileBase '-rip.mat'],'ripple_result'); %
%     MakeEvtFile([],[FileBase '.rip.evt'],[],EegRate); % make evt file for neuroscope browsing
%     return
% end

% (0) Loading eeg file for ripple detection
% rip_channels = fileinfo.ripch(1:end-1) + 1; % tetrode/silicon probe 1-4,8,16ch
% subtract_channel = fileinfo.ripch(end) + 1; % need a relative channel to subtract

rip_channels = fileinfo.RippleChannels; % tetrode/silicon probe 1-4,8,16ch
% subtract_channel = fileinfo.bestch(end) + 1; % need a relative channel to subtract



% [min_median,ix] = max(chm(1,:));
% subtract_channel = chm(2,ix);
% GammaCh = fileinfo.GammaCh(find(fileinfo.CA==1));
% 
% if ismember(subtract_channel,GammaCh)
%     GammaCh = GammaCh(find(GammaCh~=subtract_channel));
% end
    
% [sum_detection] = [];
% for shank = 1:size(bc,2);
% select_channel = [rip_channels GammaCh];  %  not yet finalized--how to select channels! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
select_channel = [rip_channels]; 
nChs = numel(select_channel);
% select_channel = bc(shank,:)';



firfiltb = fir1(forder,[lowband/EegRate*2, highband/EegRate*2]); % calculate convolution func
% avgfiltb = ones(avgfilorder,1)/avgfilorder; % passbands are normalized to Nyquist frequency.


sigma = 4/1000*EegRate; %% 4 ms s.d. of the Gaussian kernel 
halfwidth = 3 * sigma;
smoothwin = gausswindow(sigma, halfwidth);



rippleLFP = readmulti([FileBase '.eeg'],numchannel,select_channel); % load .eeg trace
%     thresholdbuffer = readlocalch([FileBase
%     '.eeg'],numchannel,select_channel,subtract_channel);

if ~isempty(period)
    rippleLFP = rippleLFP(period(1)*EegRate+1 : period(2)*EegRate, :);
end

rippleLFP = Filter0(firfiltb, rippleLFP); % the actual filtering


if thetafilter
    exclude = dlmread([FileBase '.theta.1']);
else
    exclude = [];
end


% 
% if FileExists([FileBase '.ripExclude.mat']) & isempty(without)
%     fprintf(1,'Warning!  Special ripple detection in place \n')
%     load([FileBase '.ripExclude.mat'])
% end
% exclude = sortrows([exclude;without]);
  
envelope = abs(hilbert(rippleLFP)); % to power trace
for ich = 1:nChs
    envelope(:, ich) = conv(envelope(:, ich), smoothwin, 'same');
end


% summed ripple power or envelope (this can be combined with the couple of lines above)
summed_envelope = abs(hilbert(rippleLFP)); % to power trace
summed_envelope = conv(sum(summed_envelope, 2), smoothwin, 'same');
summed_envelope = sum(summed_envelope, 2); % taking average over the channels


sdat_summed = (summed_envelope - repmat(mean(summed_envelope), [size(summed_envelope, 1), 1])) ./ repmat(std(summed_envelope), [size(summed_envelope, 1), 1]);
%%%



min_isw_period = min_isw_period/1000*EegRate; % in eeg
min_sw_period = min_sw_period/1000*EegRate; % in eeg
max_sw_period = max_sw_period/1000*EegRate;



primaryList = cell(nChs, 1);
sdat        = zeros(size(envelope, 1), nChs);

for ich = 1:nChs

    currRipplePow = envelope(:, ich);

    sdat(:, ich) = (currRipplePow - repmat(mean(currRipplePow), [size(currRipplePow, 1), 1])) ./ repmat(std(currRipplePow), [size(currRipplePow, 1), 1]);

    % [sdat, stdA] = unity(conv(sum(currRipplePow,2), smoothwin, 'same'),exclude); % averaging & standardizing
    % sdat = unity(sum(currRipplePow,2),exclude); % averaging & standardizing


    % (1) primary detection of ripple periods, based on thresholding
    thresholded = sdat(:, ich) > thresholdf;
    primary_detection_b = find([0; diff(thresholded)]>0);
    primary_detection_e = find([0; diff(thresholded)]<0);

    if primary_detection_b(1) > primary_detection_e(1)
        primary_detection_e = primary_detection_e(2:end);
    end

    currEvents =[primary_detection_b(1:size(primary_detection_e,1),1) primary_detection_e];
    
    
    %% only primary detected ripples with peak ripple amplitude passing the threshold are kept
   
    nEvents = size(currEvents, 1);
    
    keeper = zeros(nEvents, 1);
    for ii = 1:nEvents
        
        peakk = max(sdat(currEvents(ii,1):currEvents(ii,2), ich));
%         keeper(ii) = peakk >= max_thresholdf;
        keeper(ii) = peakk >= max_thresholdf & (currEvents(ii, 2)- currEvents(ii,1)) >= min_isw_period;

    end
    
    primaryList{ich} = currEvents(find(keeper), :);
    
end

primaryList = cell2mat(primaryList);


primaryList = sortrows(primaryList);
primary_detection_b = primaryList(:,1);
primary_detection_e = primaryList(:,2);

% if (length(primary_detection_e) == length(primary_detection_b)-1) % exclude ranged-out (1st or
%     primary_detection_b = primary_detection_b(1:end-1);
% end
% 
% if (length(primary_detection_e)-1 == length(primary_detection_b))
%     primary_detection_e = primary_detection_e(2:end);
% end

primary = [primary_detection_b,primary_detection_e]; % [start, end]

% 
% if isempty(primary)
%     disp(sprintf('ERROR1; no more ripple %s\n',FileBase));
%     ripple_result = [];
%     if FileExists([FileBase '-rip.mat'])
%         fprintf('save old ripple data...');
%         com = ['mv ' FileBase '-rip.mat '  FileBase '-rip.bak.mat'];
%         system(com);
%         com = ['mv ' FileBase '.rip.evt '  FileBase '.rip.bak.evt'];
%         system(com);
%     end
%     save([FileBase '-rip.mat'],'ripple_result'); %
%     MakeEvtFile([],[FileBase '.rip.evt'],[],EegRate); % make evt file for neuroscope browsing
%     return
% end

% (2) merging ripples, if inter-ripples period is less than min_isw_period;

secondary = [];


tmp_rip = primary(1,:);

for ii = 2:size(primary,1)
    if (primary(ii,1)-tmp_rip(2)) < min_isw_period
        tmp_rip = [tmp_rip(1), max(tmp_rip(2), primary(ii,2))]; % merge two ripples
%     elseif abs(primary(ii,1)-tmp_rip(1)) < min_isw_period
%         tmp_rip = [min([tmp_rip(1) primary(ii,1)]),max([tmp_rip(2)
%         primary(ii,2)])]; % I guess it's for the case if we want to merge ripples across channels 
    else
        secondary = [secondary;tmp_rip];
        tmp_rip = primary(ii,:);
    end
end

secondary = [secondary;tmp_rip]; % [start, end]

% keeper = find((secondary(:,2)-secondary(:,1)) > min_sw_period & (secondary(:,2)-secondary(:,1)) < max_sw_period);
% keeper = find((secondary(:,2)-secondary(:,1)) > min_sw_period);
% secondary = secondary(keeper,:);

% if isempty(secondary)
%     disp(sprintf('ERROR2; no more ripple %s\n',FileBase));
%     ripple_result = [];
%     if FileExists([FileBase '-rip.mat'])
%         fprintf('save old ripple data...');
%         com = ['mv ' FileBase '-rip.mat '  FileBase '-rip.bak.mat'];
%         system(com);
%         com = ['mv ' FileBase '.rip.evt '  FileBase '.rip.bak.evt'];
%         system(com);
%     end
%     save([FileBase '-rip.mat'],'ripple_result'); %
%     MakeEvtFile([],[FileBase '.rip.evt'],[],EegRate); % make evt file for neuroscope browsing
%     return
% end

% (3) ripples must have its peak power of > max_thresholdf
third = zeros(size(secondary, 1), 5);
% SDmax = [];


startT = secondary(:,1);
endT   = secondary(:,2);

for ii=1:numel(startT)
    
%     max_val = cellfun(@(x) max(x(startT(ii):endT(ii))), sdat, 'UniformOutput',false);
%     [max_val, max_ID] = max(cell2mat(max_val));
%     [~, peakID] = max(sdat{max_ID}(startT(ii):endT(ii)));
   
    currMat = sdat(startT(ii):endT(ii), :);
    currDur = endT(ii) - startT(ii) + 1;
    
    [max_val, max_ID] = max(currMat(:));
    
    peakID = mod(max_ID, currDur);
    peakCh = floor(max_ID/currDur)+1;
    
    third(ii, :) = [startT(ii) endT(ii) startT(ii)+peakID-1 max_val peakCh]; % [start, end]
    
end

third = third(third(:,4) > max_thresholdf, :);

% 
% if isempty(third)
%     disp(sprintf('ERROR3; no more ripple %s\n',FileBase));
%     ripple_result = [];
%     save([FileBase '-rip.mat'],'ripple_result'); %
%     MakeEvtFile([],[FileBase '.rip.evt'],[],EegRate); % make evt file for neuroscope browsing
%     return
% end

% (4) detection of negative peak position of each ripple
% peak_idx = zeros(size(third,1),1);
% 
% for ii=1:size(third,1)
%     [min_val,min_idx] = min(rippleLFP(third(ii,1):third(ii,2)));
%     peak_idx(ii,1) = min_idx+third(ii,1)-1;
% end
% 
% Fourth = [third(:,1),peak_idx,third(:,2),SDmax]; % [start, peak, end, SD]


% (7) save detected ripples
ripple_result = third; % [start, peak, end, SD] (in eeg)

ripple_result(:, 1:3) = ripple_result(:, 1:3)/EegRate;
% %

if homefilter & ~isempty(ripple_result)
    chts = fileinfo.chts;
    tbegin = fileinfo.tbegin;
    tend = fileinfo.ftend;
    xyt = fileinfo.xyt;
    xr2 = 0.5 + [-.13 .13];
    xr1 = 0.5 + [-.2 .2];
    tchg = (chts-tbegin)/1e6*SampleRate;

    before = find(ripple_result(:,1) < tchg);
    xbin = interp1(tt,xyt(:,1),ripple_result(before,1));
    cageout = find(xbin > max(xr1));
    cagein = find(xbin < min(xr1));

    resultb = [ripple_result(cageout,:);ripple_result(cagein,:)];

    after = find(ripple_result(:,1) > tchg);
    xbin = interp1(tt,xyt(:,1),ripple_result(after,1));

    tempvec = ripple_result(after,:);
    cageout = find(xbin > max(xr2));
    cagein = find(xbin < min(xr2));

    resulta = [tempvec(cageout,:);tempvec(cagein,:)];

    ripple_result = sortrows([resultb;resulta]);
end

if runfilter & ~isempty(ripple_result)
    
    
    % Loading POSITION data 

    xyt = fileinfo.xyt;

    if ~isfield(fileinfo, 'pix2cm')
        pix2cm = 1;
    else
        pix2cm = fileinfo.pix2cm;
    end


    withinRange = find(xyt(:,3) >= fileinfo.tbegin & xyt(:,3) <= fileinfo.tend); % selecting only the position samples occured during neural recording
    xyt = xyt(withinRange, :);

    timepnts = (xyt(:,3) - fileinfo.tbegin)/fileinfo.Fs; % u seconds to seconds

    xpos = xyt(:,1) * pix2cm; % Converting the position data (from pixel number to centimeter)
    ypos = xyt(:,2) * pix2cm;


    %%%% SPEED %%%%%%%%%%%%

    diffx = [0; abs(diff(xpos))];
    diffy = [0; abs(diff(ypos))];

    difftt = [1; diff(timepnts)];
    speed = sqrt(diffx.^2 + diffy.^2)./difftt; 

    sigma = 5; %%% smoothing the speed, the positon sampling rate is around 30 Hz, so duration of the sigma is about 500 ms
    halfwidth = 3 * sigma;
    smoothwin = gausswindow(sigma, halfwidth);
    
    speed = conv(speed, smoothwin, 'same');
    
    eventSpeed = interp1(timepnts, speed, ripple_result(:, 3)); % calculating animal's speed at the peak of events
    
    eventDur = ripple_result(:,2) - ripple_result(:,1);
   ripple_result(eventSpeed > 5, :) = [];
    
end

if thetafilter & ~isempty(ripple_result) & FileExists([FileBase 'theta.1']);
    sts = dlmread([FileBase '.theta.1'])';
    rip_ind = STSfilt(sts,ripple_result(:,2),~thetafilter);
    ripple_result = ripple_result(rip_ind,:);
end
% 
% 
% if fileinfo.NumberOfRipples~=size(ripple_result,1)
%     disp(['Previous = ' num2str(fileinfo.NumberOfRipples) ', New = ' num2str(size(ripple_result,1))])
% % ex1 = ['neuroscope ' FileBase '.eeg &']; 
% % system(ex1)
% 
% end
% 
% if FileExists([FileBase '-rip.mat'])
%     fprintf('save old ripple data...');
%     com = ['mv ' FileBase '-rip.mat '  FileBase '-rip.bak.mat'];
%     system(com);
%     com = ['mv ' FileBase '.rpl.evt '  FileBase '.rpl.bak.evt'];
%     system(com);
% end

Filename = [FileBase '.rp3.evt'];
MakeEvtFile(third(:, [1 3 2]), Filename, {'beg','peak', 'end'}, EegRate, 1)

% if ~isempty(period)
%     ripple_result(:,1:3) = ripple_result(:,1:3) + period(1)*EegRate;
% end
% 
% ripple_number = size(ripple_result,1);
% fprintf('%s  %d ripples detected successfully!\n',fileinfo.name,ripple_number);

end


