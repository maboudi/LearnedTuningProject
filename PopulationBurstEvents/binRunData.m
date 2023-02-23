function [runBinnedfiring, posbinIdx, xposcenters, yposcenters] = binRunData(mySpikes, activePeriods, behavior, binDur, posBinSize, fileinfo)


% This function is intended to divide the RUN neuronal activity to activity bouts within time bins with
% duration determined by binDur. Also, it calculates index of the position
% bin corresponding to each time bin.

% 2-D maze positions are used here


%% binning the run data


qclus = [1 2 3]; % pyramidal and interneurons

runBinnedfiring2 = timeBinning(activePeriods , mySpikes, qclus, binDur, fileinfo); % the theta periods



% Calculate the track positions corresponding to each time bin

periods = activePeriods * 1/fileinfo.lfpSampleRate;
runBinPos = cell(size(periods, 1), 1);

for ii = 1 : size(periods, 1)
    
    runBinPos{ii} = interp1(fileinfo.xyt(:, 3)*1/fileinfo.lfpSampleRate, fileinfo.xyt(:,[1 2]), (periods(ii, 1):binDur:periods(ii, 2))+binDur/2); %% at the center of each time bin
    runBinPos{ii}(end, :) = [];
end
    

% runBinPos = cell2mat(runBinPos');


% Go from the continuous position to position bins

% Defining the position bins; the same we had in calculating the place
% fields


% posBinSize = 0.5; % in cm3
runidx = find(fileinfo.xyt(:, 3) > behavior.time(2,1) & fileinfo.xyt(:, 3) < behavior.time(2,2));

xpos = fileinfo.xyt(runidx,1);
ypos = fileinfo.xyt(runidx,2);


% noXPosBins = floor((max(xpos) - min(xpos))/posBinSize); % posBinSize: the length of each position bin in cm
% noYPosBins = floor((max(ypos) - min(ypos))/posBinSize);

xposBins = min(xpos): posBinSize: max(xpos); % center of the position bins
xposBins(end) = max(xpos);

xposcenters = xposBins + posBinSize/2;
xposcenters(end) = [];


yposBins = min(ypos): posBinSize: max(ypos); 
yposBins(end) = max(ypos);

yposcenters = yposBins + posBinSize/2;
yposcenters(end) = [];


posbinIdx2 = cell(size(periods, 1), 1);
for ii = 1 : size(periods, 1)
    
    [~, posbinIdx2{ii}(:,1)] = histc(runBinPos{ii}(:,1), xposBins);
    [~, posbinIdx2{ii}(:,2)] = histc(runBinPos{ii}(:,2), yposBins);
end


runBinnedfiring{1,1} = cell2mat(runBinnedfiring2(:,1)'); % concatenating all the run binned spikes
runBinnedfiring{1,2} = cell2mat(runBinnedfiring2(:,2)'); % the same as the first column but with milisecond bins

posbinIdx = cell2mat(posbinIdx2);


end

