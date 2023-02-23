dataDirectory = '/home/kouroshmaboudi/Documents/NCMLproject/UtkuData';


sessionName = 'AG_2019-12-23_NSD';
basePath = fullfile(dataDirectory, sessionName);



%% recording info
par = LoadXml(fullfile(basePath, sessionName));


datInfo = dir([fullfile(basePath, sessionName) '.lfp']);
nlfpTimePnts = datInfo.bytes/2/par.nChannels; 


fileInfo.sessionName   = sessionName;
fileInfo.Fs            = par.SampleRate; 
fileInfo.lfpSampleRate = par.lfpSampleRate;
fileInfo.timeUnit      = 1; % in sec
fileInfo.tbegin        = [];
fileInfo.tend          = [];
fileInfo.nCh           = par.nChannels;
fileInfo.xyt           = [];
fileInfo.xyt2          = []; 



%% epochs (PRE1, PRE2 (NSD), MAZE, POST)

epochs = csvread(fullfile(basePath, [sessionName '.TimeIntervalCombined.csv']), 1, 1);
epochs_n_lfp_pnts = epochs(:, 1)/par.lfpSampleRate;

epochsTimeTable = readtable(fullfile(basePath, [sessionName '.TimeIntervalCombined.csv'])); % e-phys recordings
epochsTimeTable = table2struct(epochsTimeTable);



for ii = 1:numel(epochsTimeTable)
    currStartTime = datevec(epochsTimeTable(ii).StartTime);
    epochsTimeTable(ii).StartTime = ((currStartTime(4)*60)+currStartTime(5))*60 + currStartTime(6);
    epochsTimeTable(ii).Duration  = epochsTimeTable(ii).NumberOfPoints/epochsTimeTable(ii).SampleRate;

    epochsTimeTable(ii).EndTime   = epochsTimeTable(ii).StartTime + epochsTimeTable(ii).Duration;
end


blocks_startTimes_concatenatedData = [0 cumsum([epochsTimeTable(1:end-1).NumberOfPoints])/epochsTimeTable(1).SampleRate];
gaps_bw_blocks = [epochsTimeTable(2:end).StartTime] - ([epochsTimeTable(1:end-1).StartTime] + [epochsTimeTable(1:end-1).NumberOfPoints]/epochsTimeTable(1).SampleRate);



for ii = 1:numel(epochsTimeTable)
    if ii == 1
        behavior.time(ii,1) = blocks_startTimes_concatenatedData(ii);
    else
        behavior.time(ii,1) = blocks_startTimes_concatenatedData(ii) + sum(gaps_bw_blocks(1:ii-1));
    end
    
    behavior.time(ii,2) = behavior.time(ii,1) + epochsTimeTable(ii).Duration;
end

behavior2.time(1,1) = behavior.time(1,1);
behavior2.time(1,2) = behavior.time(2,2);
behavior2.time(2,1) = behavior.time(3,1);
behavior2.time(2,2) = behavior.time(3,2);
behavior2.time(3,1) = behavior.time(4,1);
behavior2.time(3,2) = behavior.time(4,2);

behavior = behavior2;

fileInfo.behavior      = behavior;

fileInfo.tbegin = behavior.time(1,1);
fileInfo.tend   = behavior.time(3,2);



%% zero padding the lfp/eeg file

% lpfFileName = fullfile(basePath, [sessionName '.lfp']);
% inputFileHandle = fopen(lpfFileName);
% 
% lfp_correced_filename = fullfile(basePath, [sessionName '_corrected.lfp']);
% 
% 
% nEpochs = numel(epochsTimeTable);
% 
% bufferSize = 4096;% the number of frames read at each time
% outputFileHandle = fopen(lfp_correced_filename,'w');
% 
% 
% for ii = 1:nEpochs
%         
%     if ii == 1
%         offset = 0;
%     else
%         offset = offset + epochsTimeTable(ii-1).NumberOfPoints-1;
%     end
%     
%     offset2 = offset*par.nChannels*2; % number of bytes
%     
%     frewind(inputFileHandle)
%     fseek(inputFileHandle, offset2, 'bof'); % move the file pointer to a specific point
%     
%     doneFrame = 0;
%     remainingFrames = Inf;
%     
%     while remainingFrames > 0
%         
%         nFrames = min(bufferSize, remainingFrames);
%         
%         data = (fread(inputFileHandle,[par.nChannels, nFrames],'int16'))';
%             
%         for frame=1:size(data, 1)
%             fwrite(outputFileHandle, data(frame,:)','int16');
%         end
% 
%         
%         doneFrame=doneFrame + nFrames;
%         remainingFrames = epochsTimeTable(ii).NumberOfPoints - doneFrame;
%     end
%         
%     
%     if ii < nEpochs
%         currGap = gaps_bw_blocks(ii);
%         nZeroFrames = floor(currGap * par.lfpSampleRate);
%         
%         for frame = 1:nZeroFrames
%             fwrite(outputFileHandle, zeros(par.nChannels, 1), 'int16');
%         end
%         
%     end
%    
% end
% 
% fclose(inputFileHandle);
% fclose(outputFileHandle);



%% position data

positionValues  = readtable(fullfile(basePath, [sessionName '.location.points.csv']));
positionValues  = table2struct(positionValues);

posTimeTable = readtable(fullfile(basePath, [sessionName '.location.time.xlsx']));
posTimeTable = table2struct(posTimeTable);

position_samplingRate = posTimeTable(1).SampleRate;

posTimePnts = cell(numel(posTimeTable), 1);
for ii = 1:numel(posTimeTable)   
    
    currStartTime = posTimeTable(ii).StartTime;
    currStartTime = datevec(currStartTime);
   
    startTime = ((currStartTime(4)*60)+currStartTime(5))*60 + currStartTime(6);
    posTimePnts{ii} = startTime + (1:posTimeTable(ii).NumberOfPoints)/position_samplingRate - epochsTimeTable(1).StartTime;
end

posTimePnts = cell2mat(posTimePnts');

idx = posTimePnts>0;

fileInfo.xyt = [[positionValues.x]' [positionValues.z]' posTimePnts'];
fileInfo.xyt = fileInfo.xyt(idx, :);   


speed = calspeed(fileInfo, position_samplingRate); % animal's velocity
    
fileInfo.speed = speed;
    


linearPos = linearizePosition2(fileInfo, behavior, 'L-shape');

fileInfo.xyt2 = [linearPos fileInfo.xyt(:, 3)];


direction = 'bi';
lapsStruct = calculateLapTimings(fileInfo, direction, []); 

if length(lapsStruct.RL) > length(lapsStruct.LR)
   lapsStruct.RL(1,:) = [];
end

totNumLaps = size(lapsStruct.RL, 1) + size(lapsStruct.LR, 1);
laps = zeros(totNumLaps, 2);
laps(1:2:totNumLaps, :)  = lapsStruct.LR;
laps(2:2:totNumLaps, :)  = lapsStruct.RL;


laps(:, 3) = 1:size(laps, 1); 


fileInfo.xyt2(:, 3) = zeros(size(linearPos)); % labeling the postion samples with the calculated laps (if not part of any lap the label is zero)

for ii = 1: length(laps)
   idx =  fileInfo.xyt2(:, 2) > laps(ii, 1) & fileInfo.xyt2(:, 2) < laps(ii, 2);
   fileInfo.xyt2(idx, 3) = laps(ii, 3);
end

runSpeedThresh = 10; % cm/s

%%


clu_info = readtable(fullfile(basePath, [sessionName '.clu_info.csv']));


spikeTimes  = load(fullfile(basePath, [sessionName '.res.0']));
cluContent  = load(fullfile(basePath, [sessionName '.clu.0']));

clusterIDs = cluContent(2:end);

uniqClusters = unique(clusterIDs); 
nClusters    = numel(uniqClusters);


pre_sub_edges = linspace(behavior.time(1,1), behavior.time(1,2), 4); 
stb_pre1 = [pre_sub_edges(1) pre_sub_edges(2)]; % each sub epoch is ~ 3 hours
stb_pre2 = [pre_sub_edges(2) pre_sub_edges(3)];
stb_pre3 = [pre_sub_edges(3) pre_sub_edges(4)];

stb_maze = behavior.time(2,:);

post_sub_edges = linspace(behavior.time(3,1), behavior.time(3,2), 3);
stb_post1 = [post_sub_edges(1) post_sub_edges(2)];
stb_post2 = [post_sub_edges(2) post_sub_edges(3)];

stb = [stb_pre1; stb_pre2; stb_pre3; stb_maze; stb_post1; stb_post2];



spikes = struct('time', [], 'id', [], 'amp', [], 'depth', [], 'fr', [], 'fr_by_epoch', [],...
    'fr_pre_post_ratio', [], 'fr_pre_post_bias', [], 'fr_isOK', [], 'group', [], 'n_spikes', [], 'sh', [], 'location', []);

for ic = 1:nClusters
    
    
    currCluster = uniqClusters(ic);   
    
    spikes(ic).time = spikeTimes(clusterIDs == currCluster)/fileInfo.Fs;
    
    % we need to correct the spike times, because the spike times are not
    % absolute, it was generated based on the concatenated data.
    
    % spike times, respecting the gaps 
    % |  epoch1   |%%gap1%%|  epoch2   |%%gap2%%|  epoch3   |%%gap3%%|  epoch4    
    
    % spike times in concaenated data, diregarding the gaps
    %|  epoch1   |  epoch2   |  epoch3   |  epoch4
    
    
    idx = spikes(ic).time >= blocks_startTimes_concatenatedData(4);
    spikes(ic).time(idx) = spikes(ic).time(idx) + sum(gaps_bw_blocks);
    
    idx = spikes(ic).time >= blocks_startTimes_concatenatedData(3) & spikes(ic).time < blocks_startTimes_concatenatedData(4);
    spikes(ic).time(idx) = spikes(ic).time(idx) + sum(gaps_bw_blocks(1:2));
    
    idx = spikes(ic).time >= blocks_startTimes_concatenatedData(2) & spikes(ic).time < blocks_startTimes_concatenatedData(3);
    spikes(ic).time(idx) = spikes(ic).time(idx) + gaps_bw_blocks(1);
    
    
    
    
    spikes(ic).id   = currCluster;
    
    % get info in the cluster file
    
    tempIdx = find(clu_info.id == currCluster);
    
    spikes(ic).amp      = clu_info.amp(tempIdx);
    spikes(ic).depth    = clu_info.depth(tempIdx);
    spikes(ic).fr       = clu_info.fr(tempIdx);
    spikes(ic).group    = clu_info.group(tempIdx);
    spikes(ic).n_spikes = clu_info.n_spikes(tempIdx);
    
    spikes(ic).sh       = clu_info.sh(tempIdx);
    spikes(ic).id       = [spikes(ic).sh spikes(ic).id];
    
    spikes(ic).location = clu_info.location(tempIdx);
    
    spikes(ic).fr_by_epoch = zeros(1,6); 
    for ip = 1:size(stb, 1)
       spikes(ic).fr_by_epoch(ip) = numel(find(spikes(ic).time > stb(ip, 1) & spikes(ic).time <= stb(ip, 2)))/diff(stb(ip, :));
    end
    
    spikes(ic).fr_pre_post_ratio = spikes(ic).fr_by_epoch(5)/ spikes(ic).fr_by_epoch(3); % the difference between the pre and post epoch next the the maze is enough
    
    if spikes(ic).fr_pre_post_ratio > 1
        spikes(ic).fr_pre_post_bias = 1; % smaller as a percent of the greater firing rate
        spikes(ic).fr_pre_post_ratio = 1/spikes(ic).fr_pre_post_ratio;
    else
        spikes(ic).fr_pre_post_bias = -1;
    end
    
    
    if spikes(ic).fr_pre_post_ratio > 0.3 && spikes(ic).fr_by_epoch(4) > 0.01 && strcmp(spikes(ic).group, 'good')
        spikes(ic).fr_isOK = 1;
    else
        spikes(ic).fr_isOK = 0;
    end
    
    spikes(ic).linearPos = interp1(fileInfo.xyt2(:,2), fileInfo.xyt2(:,1), spikes(ic).time);
    spikes(ic).speed = interp1(speed.t, speed.v, spikes(ic).time);
    
    spikes(ic).lap = zeros(numel(spikes(ic).time), 1);
    for spk = 1:numel(spikes(ic).time)
        index = find(laps(:, 1) < spikes(ic).time(spk) & laps(:, 2) > spikes(ic).time(spk));
        if ~isempty(index)
            spikes(ic).lap(spk) = laps(index, 3);
        end
    end

end


%% create a separate .clu and .res file for each shank




%% 

posBinSize = 2;
subfolder = fullfile(fullfile(basePath, sessionName), 'placeFields');
if ~exist(subfolder, 'dir') 
    mkdir(subfolder)
end


spikes = spatialTuning_1D_June21(spikes, behavior, [], [], speed, runSpeedThresh, 'LR', posBinSize, fileInfo, subfolder);
spikes = spatialTuning_1D_June21(spikes, behavior, [], [], speed, runSpeedThresh, 'RL', posBinSize, fileInfo, subfolder);




%% ripple detection 


fileInfo.CA = [1 1 1 1 3 3];

best_channels = BigRippleChannels(fullfile(basePath, [sessionName '_corrected']), fileInfo);
fileInfo.RippleChannels = best_channels; 





threshSD = 2;
[ripplePeriods, rippleLFP, ripplePower_summed] = RippleDetect(fullfile(basePath, [sessionName '_corrected']), fileInfo, threshSD); % the results in seconds


% downsample the ripple band lfp or power to 1kHz (1ms sampling period) % maybe we don't
% need to downsample at all?

ts    = 1/fileInfo.lfpSampleRate;

tPnts = (1:length(rippleLFP))*ts;

tq    = 1e-3:1e-3:fileInfo.tend;


rippleLFP_adjust = zeros(numel(tq), size(rippleLFP, 2));
for ich = 1:size(rippleLFP, 2)
    rippleLFP_adjust(:, ich) = interp1(tPnts, rippleLFP(:, ich), tq);
end

ripplePower_summed_adjust = interp1(tPnts, ripplePower_summed, tq);



%% adding variables to the EEG file

nCh = fileInfo.nCh;

baseFile = fullfile(basePath, [sessionName '_corrected']);
EEGdata = readmulti([baseFile '.lfp'], nCh, 1); % load .eeg trace
totalT  = length(EEGdata)/fileInfo.lfpSampleRate;


samplingDur  = 1/fileInfo.lfpSampleRate;
samplingPnts = samplingDur:samplingDur:totalT;

% sdat_mua_2add = interp1((1:length(sdat.MUA))*1e-3, sdat.MUA, samplingPnts);
% sdat_pyr_2add = interp1((1:length(sdat.PYR))*1e-3, sdat.PYR, samplingPnts);
% speed_2add    = interp1(fileInfo.speed.t, fileInfo.speed.v, samplingPnts);
ripple_2add   = interp1(tq, ripplePower_summed_adjust, samplingPnts);
% theratio_2add = interp1(timePnts, theratio, samplingPnts);


% sdat_mua_2add(isnan(sdat_mua_2add)) = 0;
% sdat_mua_2add(sdat_mua_2add < 0) = 0;
% sdat_mua_2add(sdat_mua_2add > 2) = 2;
% sdat_mua_2add = (2*(sdat_mua_2add-min(sdat_mua_2add))/(max(sdat_mua_2add)-min(sdat_mua_2add))-1)*2^10;
% 
% 
% sdat_pyr_2add(isnan(sdat_pyr_2add)) = 0;
% sdat_pyr_2add(sdat_pyr_2add < 0) = 0;
% sdat_pyr_2add(sdat_pyr_2add > 2) = 2;
% sdat_pyr_2add = (2*(sdat_pyr_2add-min(sdat_pyr_2add))/(max(sdat_pyr_2add)-min(sdat_pyr_2add))-1)*2^10;
% 
% 
% 
% speed_2add(speed_2add < 10) = 0;
% speed_2add = (2*(speed_2add-min(speed_2add))/(max(speed_2add)-min(speed_2add))-1)*2^11;


ripple_2add(isnan(ripple_2add)) = 0;
ripple_2add(ripple_2add > 3) = 3;
ripple_2add(ripple_2add < 0) = 0;
ripple_2add = (2*(ripple_2add - min(ripple_2add))/(max(ripple_2add)-min(ripple_2add))-1)*2^10;


% theratio_2add(isnan(theratio_2add)) = 0;
% theratio_2add(theratio_2add < 1) = 0;
% theratio_2add = (2*(theratio_2add - min(theratio_2add))/(max(theratio_2add)-min(theratio_2add))-1)*2^11;


inputFileHandle  = fopen([baseFile,'.lfp']);
outputFileHandle = fopen([baseFile,'-2.lfp'],'w');
doneFrame = 0;
bufferSize = 4096;
while ~feof(inputFileHandle)
    data = [fread(inputFileHandle,[nCh,bufferSize],'int16')]';
    
    for frame=1:size(data,1)
        fwrite(outputFileHandle,[data(frame,:), ...
            ripple_2add(frame+doneFrame)]','int16');
        %             theratio_2add(frame+doneFrame), ...
        %             sdat_mua_2add(frame+doneFrame), ...
        %             sdat_pyr_2add(frame+doneFrame), ...
        %             speed_2add(frame+doneFrame), ...

    end
    doneFrame=doneFrame+frame;

end

fclose(inputFileHandle);
fclose(outputFileHandle);


%%

function speed = calspeed(fileInfo, position_samplingRate)


xpos = fileInfo.xyt(:, 1);
ypos = fileInfo.xyt(:, 2);

timepnts = fileInfo.xyt(:, 3);


diffx = [0; abs(diff(xpos))];
diffy = [0; abs(diff(ypos))];

difft = [1; diff(timepnts)];


velocity = sqrt(diffx.^2 + diffy.^2)./difft; 

% interpolate the velocty at time point with missing position data
velocity(isnan(velocity)) = interp1(timepnts(~isnan(velocity)), velocity(~isnan(velocity)), timepnts(isnan(velocity)));


% smoothing the speed
sigma = 0.25/(1/position_samplingRate); 
halfwidth = 3 * sigma;


smoothwin = gausswindow(sigma, halfwidth); % smoothing kernel
speed.v = conv(velocity, smoothwin, 'same'); 

% speed.v = velocity;
speed.t = timepnts;


end
