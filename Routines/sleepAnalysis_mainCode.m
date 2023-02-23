%% Loading Grosmark's datasets


currDir = '/home/kouroshmaboudi/Documents/HMM project/Grossmark_BapunReclustered';
cd(currDir)



%% loading session data

VarList = {'spikes','behavior','position'};

% note that the spikes within the .mat file are from just the pyramidal
% units, for MUA I need to load the clu and res files again


% Although no difference regarding the behavior and position data


for var = 1 : length(VarList)
    load([currDir '/CRCNSReclustered-' VarList{var} '.mat'])
end


sessionNames = fieldnames(spikes);
noSessions = numel(sessionNames);


cov2cmfac = 100; % the position are in m in the dataset, hence positions must be multiplied with 100


for currSession = 5%: noSessions


sessionName = sessionNames{currSession};


spikes = eval(['spikes.' sessionName]);
behavior = eval(['behavior.' sessionName]);
position = eval(['Position.' sessionName]);

behavior.time = [behavior.PREEpoch; behavior.MazeEpoch; behavior.POSTEpoch];

speed = calspeed(position); % animal's velocity



%% sessioninfo

fileBase = [currDir '/' sessionName '/' sessionName];

fileinfo = struct('name', sessionName, 'animal', sessionName(1:strfind(sessionName, '_')-1), 'xyt', [], 'tbegin', [], 'tend', [], 'Fs', [], 'lfpSampleRate', [], 'nCh', []); 


% Load the recording configurations from the xml file

Par = LoadXml(fileBase);
fileinfo.lfpSampleRate = Par.lfpSampleRate; 
fileinfo.nCh = Par.nChannels;
fileinfo.Fs = 1; % Since the spike timestamps are already in seconds we don't need to use the original sampling rate of recording




% Position data

xpos = position.TwoDLocation(:,1)* cov2cmfac;
ypos = position.TwoDLocation(:,2)* cov2cmfac;


fileinfo.xyt = [xpos ypos position.TimeStamps']; 

figure; plot(xpos, ypos, '.', 'markersize', 3)



% fileinfo.xyt(find(ypos > 13 | xpos > 130 | xpos < -20), 1:2) = NaN;

% fileinfo.xyt(find(ypos > 13), 1:2) = NaN;


%% Behavioral states

% nrem = 1, Drowsy = 3.5; rem = 2, Intermediate = 1.5, wake = 4


% Based on the description of the dataset on CRCNS, the Intermediate could
% be considered a different type of NREM, characterized by high spindle (12-20 Hz) power and low movement/EMG (but considerbale change in respect to NREM).
% So, we can include them as NREM in our analyses. 
% Drowsy states happens in transition from wake to NREM, REM to NREM or within the NREM periods,
% charachterized by low overall spectral power. Can we consider them as NREM???



bvrTimeList = [behavior.Wake; behavior.Drowsy; behavior.NREM; behavior.Intermediate; behavior.REM];
bvrState = [4   *   ones(length(behavior.Wake), 1); ...
            1 *   ones(length(behavior.Drowsy), 1); ... % considering for now Drowsy as NREM
            1   *   ones(length(behavior.NREM), 1); ... 
            1   *   ones(length(behavior.Intermediate), 1); ... % considering for now Intermediate as NREM
            2   *   ones(length(behavior.REM), 1)];



%% making a new spike data just to conform with the format of Kamran's data


qual2consider = 'all';

% mySpikes = spikeBehaviorAnalysis(spikes, laps, ripple.time, speed, qual2consider, fileinfo, Fs);
mySpikes = spikeBehaviorAnalysis2(spikes, speed, qual2consider, fileinfo);




% We need the multiunit (any detected spike) as well for defining the population burst events

numberofShanks = length(dir([fileBase '.res.*']));

MUA.t = [];
for shank = 1:numberofShanks
    
    currSpikes = load([fileBase '.res.' num2str(shank)]);
    MUA.t = [MUA.t; currSpikes]; % the first is the total number of clusters detected for the shank
    
end


MUA.t = sort(MUA.t, 'ascend');
% The originla smapling frequency is 20000

MUA.t = MUA.t./20000;

