% Analyzing aRUN model quality as a function of number of states, duration of the time bins
% to control for non-specific factores, we can plot the same for shuffled
% (time-bin shuffled data)


% time bin duration
timebins = 20;
% timebins = [10 15 20 25 30 40 50 75 100 125 150 200 250 300]; % in miliseconda
nStates  = 40;
nEvents  = length(runBouts);

noActiveUnits = size(RUNbinnedPBEs.data{1,1}, 1);
nShuffle      = 500;

nPBEs = length(RUNbinnedPBEs.data);

% 
% trainLikelihood = zeros(numel(timebins), nEvents);
% cvLikelihood    = zeros(numel(timebins), nEvents);
% 
% for ii = 1: numel(timebins)
%     
%     [trainLikelihood(ii, :), cvLikelihood(ii, :)] = HMMaRUNquality(runBouts, spikeStruct, qclus, timebins(ii), nStates, fileinfo, 0);
% 
% end

[transmat0, lambda0, prior0] = initialModel(nStates, noActiveUnits);

% nPBEs = 385;

transmatActiveRUN = cell(numel(timebins), 1);
lambdaActiveRUN  = cell(numel(timebins), 1);

HMMSeqScores = zeros(nPBEs, numel(timebins), 2);
dataLL       = zeros(nPBEs, numel(timebins), 2); %length(RUNbinnedPBEs.data)


for ii = 1: numel(timebins)
    
    fprintf('\n running the calculations using active RUN time bins with duration of %d ms\n', timebins(ii))
    % time-binning the events (here the theta bouts/high speed bouts instead of PBEs)

    runBoutsBinnedfiring  = timeBinning(runBouts, spikeStruct, qclus, timebins(ii)/1000, fileinfo);
    
    % exclude the time bins overlaping with ripples
    
    % ADD THE CODE HERE FOR EXCLUDING RIPPLE PERIODS
    
    runBoutsBinnedfiring  = removeSideSilentBins(runBoutsBinnedfiring, runBouts(:, 1), timebins(ii)/1000, fileinfo);

    [~, prior, transmatActiveRUN{ii}, lambdaActiveRUN{ii}] = phmm_em(runBoutsBinnedfiring(:,2), prior0, transmat0, lambda0, 'verbose', 1, 'max_iter', 100);
    
    [transmatActiveRUN{ii}, lambdaActiveRUN{ii}] = sortStates(transmatActiveRUN{ii}, lambdaActiveRUN{ii}, prior);
    
%     
%     gammaDist = cell(nPBEs, 1);
%     for pbe = 1:nPBEs
%    
%         currEvent = RUNbinnedPBEs.data{pbe, 2};
%         B = poisson_prob(currEvent, lambdaActiveRUN{ii},1);
%         prior = 1/nStates * ones(nStates,1); %%% a uniform prior probability distribution over the states
% 
%         %%% observation likelihood for the actual data
%         [~, ~, gammaDist{pbe},~, ~] = fwdback(prior, transmatActiveRUN{ii}, B);
%     end

     
%     
%     [transmatActiveRUN{ii}, lambdaActiveRUN{ii}] = sortStates(transmatActiveRUN{ii}, lambdaActiveRUN{ii}, prior);
%     
    for sm = 1:2
        [HMMSeqScores(:, ii, sm), ~, dataLL(:, ii, sm)]  = modelCongruence(RUNbinnedPBEs.data, transmatActiveRUN{ii}, lambdaActiveRUN{ii}, nShuffle, directory, 'aRUN', 'RUN', 0, sm);
    end
end


transmatScores = HMMSeqScores(:,:,1);
BDscores = BDseqscore.RUN.data.wPBEtimeswap.replayScore.prctilescore;
figure;
for ii = 1:numel(timebins)
    
    currHMMscores = transmatScores(:, ii);
    
    subplot(7,2,ii)
    plot(BDscores, currHMMscores, '.', 'markersize', 3, 'color', 'k')
    
    xlim([0 100])
    ylim([0 100])
    
    xlabel('time swap radon transform scores', 'fontsize', 10)
    ylabel('HMM tmat shuffle scores', 'fontsize', 10)
    
    title(sprintf('time bin size = %d', timebins(ii)))
    axis square
    
    actualDiff = mean(abs(BDscores- currHMMscores));
    
    nIter = 1000;
    shuffDiff = zeros(nIter, 1);
    for jj = 1:nIter; shuffDiff(jj) = mean(abs(BDscores(randperm(nPBEs))-currHMMscores(randperm(nPBEs))));end
    
    text(10, 10, sprintf('scoreDiff = %.0f \n(%.0f percentile of shuffle)', actualDiff, length(find(shuffDiff > actualDiff))/nIter * 100), 'color', 'r')
    
end
    
    

%%
figure;
set(gcf, 'Units', 'centimeters', 'position', [0 0 25 21])


subplot(223); customCdfPlot(timebins, HMMSeqScores(:, : ,1)); xlabel('HMM sequence score', 'fontsize', 12); ylabel('Cumulative ratio', 'fontsize', 12); title('t-mat shuffle', 'fontsize', 12)
subplot(224); customCdfPlot(timebins, HMMSeqScores(:, : ,2)); xlabel('HMM sequence score', 'fontsize', 12); ylabel('Cumulative ratio', 'fontsize', 12); title('within-PBE time swap', 'fontsize', 12)

subplot(221); customCdfPlot(timebins, dataLL(:, : ,1));    xlabel('Likelihood', 'fontsize', 12); ylabel('Cumulative ratio', 'fontsize', 12); title('likelihood', 'fontsize', 12)

subDir = fullfile(directory, 'replayHMM');
mkdir(subDir);

filename = fullfile(subDir, 'activeRUNTimeBinSize_testedOnPOSTPBEs');

savepdf(gcf, filename, '-dsvg')
saveas(gcf, filename, 'epsc')

% save(fullfile(subDir, 'activeRUNTimeBinSize_testedOnPOSTPBEs.mat'), 'HMMSeqScores', 'dataLL')




%% transmat and lambda matrices corresponding to different time bin size of active RUN  

figure;
set(gcf, 'units', 'centimeters', 'position', [0 0 19 27])

for ii = 1:numel(timebins)/2
    
    subplot(numel(timebins)/2, 5, (ii-1)*5+1)
    imagesc(transmatActiveRUN{ii})
    
    set(gca, 'xticklabel', [], 'yticklabel', [])
    xlabel('state j', 'fontsize', 10)
    ylabel('state i', 'fontsize', 10)
    
    title(sprintf('bin size = %d', timebins(ii)), 'fontsize', 10)
%     axis square
    
    subplot(numel(timebins)/2, 5, (ii-1)*5+2)
    imagesc(lambdaActiveRUN{ii}-repmat(mean(lambdaActiveRUN{ii}, 2), [1, nStates]))
    set(gca, 'YDir', 'normal', 'xticklabel', [], 'yticklabel', [])
    xlabel('state', 'fontsize', 10)
    ylabel('unit', 'fontsize', 10)
%     axis square
    
    %%%
    subplot(numel(timebins)/2, 5, (ii-1)*5+4)
    imagesc(transmatActiveRUN{ii+numel(timebins)/2})
    
    set(gca, 'xticklabel', [], 'yticklabel', [])
    xlabel('state j', 'fontsize', 10)
    ylabel('state i', 'fontsize', 10)
    
    title(sprintf('bin size = %d', timebins(ii+numel(timebins)/2)), 'fontsize', 10)
%     axis square
    
    subplot(numel(timebins)/2, 5, (ii-1)*5+5)
    imagesc(lambdaActiveRUN{ii+numel(timebins)/2}-repmat(mean(lambdaActiveRUN{ii+numel(timebins)/2}, 2), [1, nStates]))
    set(gca, 'YDir', 'normal', 'xticklabel', [], 'yticklabel', [])
    xlabel('state', 'fontsize', 10)
    ylabel('unit', 'fontsize', 10)
%     axis square
    
end


colormap('jet')







%% functions

function customCdfPlot(timeBinSizes, data)


[~, ind] = max(median(data, 1));

nPBEs         = size(data, 1);
nTimeBinSizes = numel(timeBinSizes);

cm = colormap('jet');
set(gca, 'linewidth', 1.5)
hold on

for ii = 1: nTimeBinSizes
   
    x = data(:, ii);
    x = sort(x, 'ascend');
    
    notDup = ([diff(x(:)); 1] > 0);
    xCDF   = x(notDup);
    
    yCDF = (1:nPBEs)'/nPBEs;
    yCDF = [0; yCDF(notDup)];
    
    
    k = length(xCDF);
    n = reshape(repmat(1:k, 2, 1), 2*k, 1);
    xCDF    = [-Inf; xCDF(n); Inf];
    yCDF    = [0; 0; yCDF(1+n)];
    
    
    
    cl = cm(ceil(ii/nTimeBinSizes*64), :);
    
    if ii == ind
       disname = sprintf('%d ms (max median)', timeBinSizes(ii));
       lw = 4;
    else
       disname = sprintf('%d ms', timeBinSizes(ii));
       lw = 2;
    end
    
    plot(xCDF, yCDF, 'color', cl , 'linewidth', lw, 'DisplayName', disname)

    
end



xlim([min(data(:))-0.05*range(data(:)) max(data(:))])
ylim([-0.05 1])

legend show
legend('location', 'northwest', 'NumColumns', 1, 'fontsize', 8, 'box', 'off')

axis square


end


%% 




% 
% 
% 
% figure;
% % transmat shuffle scores
% 
% subplot(221); 
% 
% fill([timebins fliplr(timebins)]', [prctile(HMMSeqScores(:, :, 1), 25, 1) fliplr(prctile(HMMSeqScores(:, :, 1), 75, 1))]' , 'b')
% hold on; plot(timebins, median(HMMSeqScores(:, :, 1), 1), 'k', 'linewidth', 3)
% xlabel('Active RUN bin size (ms)', 'fontsize', 14); ylabel('Sequence Score', 'fontsize', 14); title('t-mat shuffle', 'fontsize', 14)
% 
% 
% subplot(223); fill([timebins fliplr(timebins)]', [prctile(dataLL(:, :, 1), 25, 1) fliplr(prctile(dataLL(:, :, 1), 75, 1))]' , 'g')
% hold on; plot(timebins, median(dataLL(:, :, 1), 1), 'k', 'linewidth', 3)
% xlabel('Active RUN bin size (ms)', 'fontsize', 14); ylabel('PBE likelihood', 'fontsize', 14); 
% 
% 
% subplot(222); fill([timebins fliplr(timebins)]', [prctile(HMMSeqScores(:, :, 2), 25, 1) fliplr(prctile(HMMSeqScores(:, :, 2), 75, 1))]' , 'b')
% hold on; plot(timebins, median(HMMSeqScores(:, :, 2), 1), 'k', 'linewidth', 3)
% xlabel('Active RUN bin size (ms)', 'fontsize', 14); ylabel('Sequence Score', 'fontsize', 14); title('within-PBE time swap', 'fontsize', 14)

