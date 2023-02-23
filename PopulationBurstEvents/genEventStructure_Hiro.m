function varargout = genEventStructure_Hiro(EVTs, spikes_pyr, sdat, detectParams, fileInfo, folderName)


sessionName     = fileInfo.sessionName;
behavior        = fileInfo.behavior;
speed           = fileInfo.speed;
position        = fileInfo.linearPos;
brainStates     = fileInfo.brainStates;

states = {'NREM'; 'REM'; 'QW'; 'WAKE'};

EVTInfo = struct('sessionName', [], 'n', [], ...
                    'startT', [], 'endT', [], 'peakT', [], ...
                    'peakMUA', [], 'MUAWave', [], ...
                    'nFiringUnits', [], 'duration', [], ...
                    'nTimeBins', [], ...
                    'peakRippleA', [], 'rippleWave', [], ...
                    'thetaRatio', [], 'SWA', [], 'emg', [], ...
                    'epoch', [], 'brainState', [], ...
                    'fr_1msbin', [], 'fr_20msbin', []);


for evt = 1:numel(EVTs)

    currEVT = EVTs(evt);


    EVTInfo(evt).sessionName = sessionName;

    if isfield(currEVT, 'n')
        EVTInfo(evt).n       = currEVT.n;
    else
        EVTInfo(evt).n       = evt;
    end

    EVTInfo(evt).startT      = currEVT.startT;
    EVTInfo(evt).endT        = currEVT.endT;
    EVTInfo(evt).peakT       = currEVT.peakT;

    EVTInfo(evt).duration    = EVTInfo(evt).endT - EVTInfo(evt).startT; 

    EVTInfo(evt).peakMUA     = currEVT.peakMUA;


    try
        MUAWave = sdat(floor((EVTInfo(evt).peakT*1e3)-700) : (floor(EVTInfo(evt).peakT*1e3)+700), 1);
        temp    = nan(size(MUAWave));
        temp((floor((EVTInfo(evt).startT*1e3))  - (floor(EVTInfo(evt).peakT*1e3) - 700)) + (1:floor(EVTInfo(evt).duration*1e3))) = 0;

        EVTInfo(evt).MUAWave(:, 1) = MUAWave;
        EVTInfo(evt).MUAWave(:, 2) = temp;

    catch
        EVTInfo(evt).MUAWave = [];
    end


    EVTInfo(evt).linearPos   = interp1(position(:, 2), position(:, 1), EVTInfo(evt).peakT);

    EVTInfo(evt).velocity    = interp1(speed.t, speed.v, EVTInfo(evt).peakT);
    EVTInfo(evt).velocity(isnan(EVTInfo(evt).velocity)) = -1;


    % divding the event to 20 ms time bins
    [binnedEVT, nFiringUnits, EVTlength] = timeBinning(currEVT, spikes_pyr, detectParams.binDur, fileInfo);

    EVTInfo(evt).nTimeBins   = EVTlength;
    EVTInfo(evt).nFiringUnits = nFiringUnits;
    EVTInfo(evt).fr_1msbin   = binnedEVT{1};
    EVTInfo(evt).fr_20msbin  = binnedEVT{2};


    if EVTInfo(evt).startT > behavior.time(1,1) && EVTInfo(evt).startT < behavior.time(1,2)
       EVTInfo(evt).epoch = 'pre';
    elseif EVTInfo(evt).startT > behavior.time(2,1) && EVTInfo(evt).startT < behavior.time(2,2)
       EVTInfo(evt).epoch = 'run';
    elseif EVTInfo(evt).startT > behavior.time(3,1) && EVTInfo(evt).startT < behavior.time(3,2)
       EVTInfo(evt).epoch = 'post';
    end



    % estimating brain states of each EVT

    idx = find(EVTInfo(evt).peakT > brainStates(:, 1) & EVTInfo(evt).peakT < brainStates(:, 2), 1, 'first');
    EVTInfo(evt).brainState = states{brainStates(idx, 3)};

    if evt == 1
        fprintf('\n dividing the EVTs into 20 ms time bins ..')
    end

    if mod(evt, 1000) == 0
        fprintf('.')
    end

end


%     PBEInfo = EVTInfo(ismember({EVTInfo.brainState}, {'NREM'; 'QW'}) & ...
%                             [EVTInfo.velocity] < detectParams.maxVelocity);               
PBEInfo = EVTInfo;


% generate event files of vidualizing on Neuroscope (all the EVTs not only the accepted ones for Bayesian decoding)
nPBEs = numel(PBEInfo);

event_neuroscope       = zeros(nPBEs, 3);
event_neuroscope(:, 1) = [PBEInfo.startT];
event_neuroscope(:, 2) = [PBEInfo.endT];
event_neuroscope(:, 3) = [PBEInfo.peakT];

% make an .evt file for the EVTs
fileName = fullfile(folderName, [fileInfo.sessionName '.mua.evt']);
MakeEvtFile(event_neuroscope(:, 1:3), fileName, {'beg', 'end', 'peak'}, 1, 1)



% A smaller subset of PBEs with stricter criteria, on which we are going to
% apply Bayesian decoding
idx_Bayesian = [PBEInfo.nFiringUnits] >= detectParams.minFiringUnits & ...
                [PBEInfo.nTimeBins] >= [detectParams.minNTimeBins];
PBEInfo_Bayesian = PBEInfo(idx_Bayesian);



% save a .mat file contatining PBEs that passed the criteria
fileName = fullfile(folderName, [fileInfo.sessionName '.PBEInfo.mat']);
save(fileName, 'PBEInfo_Bayesian', 'detectParams', 'sdat', '-v7.3')


 % store all the PBEs before doing further filterings based on duration
% or unit participation
fileName = fullfile(folderName, [fileInfo.sessionName '.PBEInfo_all.mat']);
save(fileName, 'PBEInfo', 'idx_Bayesian', '-v7.3')


if nargout  > 0 ; varargout{1} = PBEInfo; end
if nargout  > 1 ; varargout{2} = PBEInfo_Bayesian; end
if nargout  > 2 ; varargout{3} = idx_Bayesian;end 


end