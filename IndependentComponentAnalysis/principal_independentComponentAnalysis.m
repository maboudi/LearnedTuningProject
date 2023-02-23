% This script in intended to calcualte ICs and their corresponding weights
% from RUN period and their activation strengths during sleep periods: PRE and POST. 
% This analysis is based on two time scales: 250ms to detect reactivations of coarse
% time-scale assemblies and 20ms to detect fine time-scale firing regimes,
% each possiblty corresponding to a track position or in general a neuronal
% population state. 


%% RUN principle-independent components (corresponding weights)

% binning RUN data

binDur = 0.25; % in sec
posBinSize = 2; % in cm

runSpeedThresh = 0;
[runBinnedfiring, posbinIdx, linposcenters] = binRunData_1D(spikeStruct, behavior.time(2, :), behavior, speed, runSpeedThresh, binDur, posBinSize, fileinfo);
runBinnedfiring = runBinnedfiring{:};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the RUN ICs and IC weights

[ICs, ICweights] = pca_ica(runBinnedfiring);
nICs = size(ICs, 1);
randomPosbinIdx  = posbinIdx(randperm(length(posbinIdx(:, 1))), 1);


for ii = 1:nICs; ICweights(ii, :) = ICweights(ii, :)/norm(ICweights(ii, :)); end

% spatial tunings

uniqPositions = unique(posbinIdx(:, 1));

ICtunings         = zeros(size(ICs, 1), numel(uniqPositions));
ICtunings_shuffle = zeros(size(ICs, 1), numel(uniqPositions));

for ii = 1: numel(uniqPositions)
    ICtunings(:, ii)         = sum(ICs(:, posbinIdx(:, 1) == uniqPositions(ii)), 2);
    ICtunings_shuffle(:, ii) = sum(ICs(:, randomPosbinIdx == uniqPositions(ii)), 2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% recalculate the ICs by dividing the whole set of RUN time bins to two
% parts, train and test sets, and do the cross-validation to evaluate their
% activation strengths and spatial tuning

nRUNtimebins = size(runBinnedfiring, 2);
trainSet     = randperm(nRUNtimebins, floor(nRUNtimebins/2));
testSet      = setdiff(1:nRUNtimebins, trainSet);

[ICs_CV, ICweights_CV] = pca_ica(runBinnedfiring(:, trainSet));
nICs_CV = size(ICs_CV, 1);

for ii = 1:nICs_CV; ICweights_CV(ii, :) = ICweights_CV(ii, :)/norm(ICweights_CV(ii, :)); end

activationStrength = struct('PRE', [], 'RUN', [], 'POST', [], 'PREn', [], 'POSTn', []);

activationStrength.RUN = calActStrength(ICweights_CV, runBinnedfiring(:, testSet));

% spatial tunings CV

posbinIdx_CV       = posbinIdx(trainSet, 1);
randomPosbinIdx_CV = posbinIdx_CV(randperm(length(posbinIdx_CV)), 1); 

ICtunings_CV         = zeros(nICs_CV, numel(uniqPositions));
ICtunings_CV_shuffle = zeros(nICs_CV, numel(uniqPositions));

for ii = 1: numel(uniqPositions)
    ICtunings_CV(:, ii)         = sum(activationStrength.RUN(:, posbinIdx_CV == uniqPositions(ii)), 2);
    ICtunings_CV_shuffle(:, ii) = sum(activationStrength.RUN(:, randomPosbinIdx_CV == uniqPositions(ii)), 2);
end



%% Sleep reactivations

qclus = [1 2 3]; % pyramidal and interneurons
binOverlapRatio = 0.5;
stepSize = (1-binOverlapRatio)*binDur;
avgWinSize = 15*60; % in seconds the duration of each sleep bout for averaging the reactivation strengths
avgWinOverlapRatio = 2/3;
avgWinStepSize = (1-avgWinOverlapRatio)*avgWinSize;

% PRE
binCentersPRE = (behavior.time(1, 1):stepSize:(behavior.time(1, 2)-binDur)) + binDur/2;
preBinnedfiring = timeBinning_withOverlap(behavior.time(1, :) , spikeStruct, qclus, binDur, binOverlapRatio, fileinfo); % the theta periods
activationStrength.PRE = calActStrength(ICweights, preBinnedfiring{1,2}); 
activationStrength.PREn = (activationStrength.PRE - repmat(mean(activationStrength.PRE, 2), [1 size(activationStrength.PRE, 2)]))./repmat(std(activationStrength.PRE, [], 2), [1 size(activationStrength.PRE, 2)]);


PREavgWins(:, 1) = (behavior.time(1, 1):avgWinStepSize:(behavior.time(1, 2)-avgWinSize));
PREavgWins(:, 2) = PREavgWins(:, 1) + avgWinSize;
avgWinCentersPRE = PREavgWins(:, 1) + avgWinSize/2;

navgWins = size(PREavgWins, 1); 

% navgWins   = ceil(range(behavior.time(1,:))/avgWinSize);
% [~, idx] = histc(binCentersPRE, PREavgWins);

PREactTimeCourse = zeros(nICs, navgWins);

for ii = 1: navgWins
    binsIncluded = binCentersPRE > PREavgWins(ii, 1) & binCentersPRE < PREavgWins(ii, 2);  
    PREactTimeCourse(:, ii) = mean(activationStrength.PREn(:, binsIncluded), 2); % should I average over absolute strengths or the actual ones
end



% POST
binCentersPOST = (behavior.time(3, 1):stepSize:(behavior.time(3, 2)-binDur)) + binDur/2;
postBinnedfiring = timeBinning_withOverlap(behavior.time(3, :) , spikeStruct, qclus, binDur, binOverlapRatio, fileinfo); % the theta periods
activationStrength.POST = calActStrength(ICweights, postBinnedfiring{1,2}); 
activationStrength.POSTn = (activationStrength.POST - repmat(mean(activationStrength.PRE, 2), [1 size(activationStrength.POST, 2)]))./repmat(std(activationStrength.PRE, [], 2), [1 size(activationStrength.POST, 2)]);


POSTavgWins(:, 1) = (behavior.time(3, 1):avgWinStepSize:(behavior.time(3, 2)-avgWinSize));
POSTavgWins(:, 2) = POSTavgWins(:, 1) + avgWinSize;
avgWinCentersPOST = POSTavgWins(:, 1) + avgWinSize/2;

navgWins = size(POSTavgWins, 1); 

POSTactTimeCourse = zeros(nICs, navgWins);

for ii = 1: navgWins
    binsIncluded = binCentersPOST > POSTavgWins(ii, 1) & binCentersPOST < POSTavgWins(ii, 2);  
    POSTactTimeCourse(:, ii) = mean(activationStrength.POSTn(:, binsIncluded), 2); % should I average over absolute strengths or the actual ones
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot IC weights(or bases), and their spatial tunings
% recalculate the ICs and plot average activation strength and sort the ICs accordingly

figure;
set(gcf, 'position', [2500 700 620 750])

cm = redblue;
fontsize = 10;

subplot(311)
imagesc(ICweights', [-max(abs(ICweights(:))) max(abs(ICweights(:)))]); colormap(cm); colorbar; 
set(gca, 'ydir', 'normal', 'fontsize', fontsize)
xlabel('independent components', 'fontsize', 12)
ylabel('units', 'fontsize', 12)

subplot(312)
ICtunings2 = zscore(ICtunings(:, 2:end), [], 2); % excluding the zero position bin and z scoring
imagesc(1:size(ICtunings2, 2)*2, 1:nICs,  ICtunings2, [-max(abs(ICtunings2(:))) max(abs(ICtunings2(:)))]); colormap(cm); colorbar;
set(gca, 'ydir', 'normal', 'fontsize', fontsize)
xlabel('track position(cm)', 'fontsize', 12)
ylabel('independent components', 'fontsize', 12)

subplot(313)
ICtunings2_shuffle = zscore(ICtunings_shuffle(:, 2:end), [], 2); % excluding the zero position bin and z scoring
imagesc(0:size(ICtunings2_shuffle, 2)*2, 1:nICs,  ICtunings2_shuffle, [-max(abs(ICtunings2(:))) max(abs(ICtunings2(:)))]); colormap(cm); colorbar;
set(gca, 'ydir', 'normal', 'fontsize', fontsize)
xlabel('track position(cm)', 'fontsize', 12)
ylabel('independent components', 'fontsize', 12)

% savepdf(gcf, 'ICactivationStrength_AchillesLinear_20msICs', '-dpdf')
%%%%%%%%%%%%%%%%%%
% cross-validation

figure;
set(gcf, 'units', 'points', 'position', [100 200 1200 700])

subplot(4,5,1)
imagesc(ICweights_CV', [-max(abs(ICweights_CV(:))) max(abs(ICweights_CV(:)))]); colormap(cm); colorbar
set(gca, 'ydir', 'normal', 'fontsize', fontsize)
xlabel('independent components', 'fontsize', 8)
ylabel('units', 'fontsize', 8)

subplot(4,5,6)
ICtunings2_CV = zscore(ICtunings_CV(:, 2:end), [], 2); % excluding the zero position bin and z scoring
imagesc(1:size(ICtunings2_CV, 2)*2, 1:nICs,  ICtunings2_CV, [-max(abs(ICtunings2_CV(:))) max(abs(ICtunings2_CV(:)))]); colormap(cm); colorbar;
set(gca, 'ydir', 'normal', 'fontsize', fontsize)
xlabel('track position(cm)', 'fontsize', 8)
ylabel('independent components', 'fontsize', 8)
title('activation strength (CV)', 'fontsize', 8)

subplot(4,5,11)
ICtunings2_CV_shuffle = zscore(ICtunings_CV_shuffle(:, 2:end), [], 2); % excluding the zero position bin and z scoring
imagesc(1:size(ICtunings2_CV_shuffle, 2)*2, 1:nICs,  ICtunings2_CV_shuffle, [-max(abs(ICtunings2_CV(:))) max(abs(ICtunings2_CV(:)))]); colormap(cm); colorbar;
set(gca, 'ydir', 'normal', 'fontsize', fontsize)
xlabel('track position(cm)', 'fontsize', 8)
ylabel('independent components', 'fontsize', 8)
title({'activation strength (CV)';'shuffled position'}, 'fontsize', 8)

%%%%%%%%%%%%%%%%
NREMperiods = bvrTimeList(ismember(bvrState, 1),:);

% PRE
ax1 = subplot(4,5, 2:5);
imagesc(binCentersPRE, 1:nICs, activationStrength.PRE, [0 10]); colormap(gca, flipud(gray)); 
PRENREMperiods = NREMperiods(NREMperiods(:, 1) > behavior.time(1, 1) & NREMperiods(:, 2) < behavior.time(1, 2), :)-behavior.time(1,1);

hold on
for ii = 1:size(PRENREMperiods,1)
patch([PRENREMperiods(ii, 1) PRENREMperiods(ii, 2) PRENREMperiods(ii, 2) PRENREMperiods(ii, 1)], [0 0 nICs+1 nICs+1], 'b', 'facealpha', 0.1, 'edgeColor', 'none')
end

xlim(behavior.time(1,:))

set(gca, 'ydir', 'normal', 'fontsize', fontsize)
xlabel('time(sec)', 'fontsize', 8)
ylabel('independent components', 'fontsize', 8)
title({'PRE';'absolute activation strength'})

ax2 = subplot(4,5, 7:10);
h1 = plot(avgWinCentersPRE , PREactTimeCourse', 'linewidth', 2);

hold on
plot(avgWinCentersPRE , mean(PREactTimeCourse, 1)', 'linewidth', 3, 'color', 'k');

xlim(behavior.time(1,:))

xlabel('time(sec)', 'fontsize', 8)
ylabel('independent components', 'fontsize', 8)
title('average activation strength')

linkaxes([ax1 ax2], 'x')

% POST
ax3 = subplot(4,5, 12:15);
imagesc(binCentersPOST-behavior.time(3,1), 1:nICs, activationStrength.POST, [0 10]); colormap(gca, flipud(gray)); 
POSTNREMperiods = NREMperiods(NREMperiods(:, 1) > behavior.time(3, 1) & NREMperiods(:, 2) < behavior.time(3, 2), :)-behavior.time(3,1);

hold on
for ii = 1:size(POSTNREMperiods,1)
patch([POSTNREMperiods(ii, 1) POSTNREMperiods(ii, 2) POSTNREMperiods(ii, 2) POSTNREMperiods(ii, 1)], [0 0 nICs+1 nICs+1], 'b', 'facealpha', 0.1, 'edgeColor', 'none')
end

xlim(behavior.time(3,:)-behavior.time(3,1))

set(gca, 'ydir', 'normal', 'fontsize', fontsize)
xlabel('time(sec)', 'fontsize', 8)
ylabel('independent components', 'fontsize', 8)
title({'POST';'absolute activation strength'})

oldPosition = get(gca, 'position');
set(gca, 'position', [oldPosition(1) oldPosition(2) oldPosition(3)*length(postBinnedfiring{1,2})/length(preBinnedfiring{1,2}) oldPosition(4)])

ax4 = subplot(4,5, 17:20);
h2 = plot(avgWinCentersPOST-behavior.time(3,1) , POSTactTimeCourse', 'linewidth', 2);

hold on
plot(avgWinCentersPOST-behavior.time(3,1) , mean(POSTactTimeCourse, 1)', 'linewidth', 3, 'color', 'k');


xlim(behavior.time(3,:)-behavior.time(3,1))

xlabel('time(sec)', 'fontsize', 8)
ylabel('independent components', 'fontsize', 8)
title('average activation strength')

oldPosition = get(gca, 'position');
set(gca, 'position', [oldPosition(1) oldPosition(2) oldPosition(3)*length(postBinnedfiring{1,2})/length(preBinnedfiring{1,2}) oldPosition(4)])

linkaxes([ax3 ax4], 'x')
linkaxes([ax2 ax4], 'y')

for ii = 1:nICs
    h1(ii).Color = [h1(ii).Color 0.5];
    h2(ii).Color = [h1(ii).Color 0.5];
end

% savepdf(gcf, 'ICactivationStrength_AchillesLinear_20msICs_CV', '-dpdf')

function activationStrength = calActStrength(ICweights, x)

x = zscore(x, [], 2);
nTimebins = size(x,2);
nICs = size(ICweights, 1); % rows correspond to weights for each IC

activationStrength = zeros(nICs, nTimebins);

for ii = 1:nICs
    currentICbase = ICweights(ii, :)';
    ProjectorOp = currentICbase * currentICbase';
    ProjectorOp = ProjectorOp - ProjectorOp.*eye(size(ProjectorOp)); % removing the diagonal to prevent high activation strength casued by the isolated activity of a single neuron with high weight in the pattern
    
    for jj = 1:nTimebins
        activationStrength(ii, jj) = x(:, jj)' * ProjectorOp * x(:, jj);
    end
end

end

