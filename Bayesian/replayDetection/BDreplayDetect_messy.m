function [BDprctile, PosteriorProbMatrix, begPosition, endPosition, slope, replaydirection, replayOrder] = BDreplayDetect(eventsBinnedfiring, tuningRL, tuningLR, activeUnits, binDur, FileBase, dataType)

% This function is intended to calculate the bayesian replay score for each
% event. The significance of the replay scores is evaluated through
% comparing it with a distribution of scores resulted from shuffle
% datasets, i.e. column cycle shuffle, time bin shuffle, etc. 



noEvents = size(eventsBinnedfiring, 1);
noPositionBins = size(tuningRL, 2);

noshufflingMethods = 2; %% number of methods for shuffling (time swap and column cycle shuffle here)
noShuffle = 200; %% number of shuffles in each of the methods 


if isempty(activeUnits)
   activeUnits = 1:size(tuningRL, 1); 
end

%%%%%%%%%%%% Initiation %%%%%%%%%%%%%


%%% posterior probability matrices
%%% the matrix comprised of posterior probability of the track positions (as rows) within each time bins of the event (as columns)


PosteriorProbMatrix = cell(noEvents,3); 

%%% 1st array -->  leftward(RL) trajectory only
%%% 2nd array --> rightward(LR) trajectory only
%%% 3rd array -->  summed over the two directions
%%% we can use either the sum or the two trajectories separately for replay
%%% detection (here we have used just the sum)



%%% best fitting line
slope = zeros(noEvents, 3); %% RL, LR, and both
rho = zeros(noEvents,3);  
begPosition = zeros(noEvents,2,3);   %% y intercept; the third dimention is comprised of RL,LR, and both, respectively
endPosition = zeros(noEvents,2,3);  

%%% replay score by Bayesian (goodness of fit of the best fitting line)
replayScore = zeros(noEvents,3);


%%% replay type
replaydirection = zeros(noEvents,1); %% if above zero, current event is replaying rightward (LR) trajectory template and if less than zero, leftward (RL) template
replayOrder = zeros(noEvents, 1); %% if above zero is a forward replay and if less than zero is a reverse replay

%%% the matrix for shuffled data 
% 
% replayScore_null = cell(1,3);
% 
% for ii = 1:3
%     replayScore_null{ii} = zeros(noEvents, noShuffle, noshufflingMethods);
% end
% 

%%% significance evaluation
BDprctile = zeros(noEvents, noshufflingMethods, 3); 
% pval = zeros(noEvents, noshufflingMethods, 3);


for evt =  1: noEvents
    
    currEvent = eventsBinnedfiring{evt, 2}(activeUnits, :); %%% using just the active units
    
    
    % exclude the time bins with no unit firing
    
    sumFiring = sum(currEvent, 1);
    idx = find(sumFiring == 0);
    
    %%% Bayesian Sequential Score

    %%% calculate the posterior probabilty matrix without any normalization 
    
    postprobRL = baysDecoder(currEvent, tuningRL(activeUnits, :), binDur); 
    postprobLR = baysDecoder(currEvent, tuningLR(activeUnits, :), binDur);

    
    postprobRL = postprobRL + eps*eps;
    postprobLR = postprobLR + eps*eps;
    
    
    %%% Summed probability distribution, every column is normalized to one
    %%% (this is the same method as used in Davidson et al. 2009)
    
%     PosteriorProbMatrix{evt, 3} = (postprobRL + postprobLR)./repmat(sum(postprobRL + postprobLR, 1), [noPositionBins, 1]);
%     PosteriorProbMatrix{evt, 3}(:, idx) = 0;

    %%% calculate the best fitting line and the replay score and replay
    %%% type for the current event
                       
%     [replayScore(evt, 3), slope(evt, 3), rho(evt, 3), begPosition(evt,:,3), endPosition(evt,:,3)] = replayevaluation(PosteriorProbMatrix{evt, 3});
    
    %%% calculating replay order(Forward or reverse) and replay direction
    %%% (LR or RL)
    
    
%     [replaydirection(evt), replayOrder(evt)] = replaytype(postprobRL, postprobLR, slope(evt, 3), rho(evt, 3)); 
    
    
    %%% normalize posterior probabilty matrices of each direction in case we are interested in
    %%% doing replay evaluation for each direction independently (as used by Grosmark and Buzsaki 2015 science paper)
        
    PosteriorProbMatrix{evt, 1} = postprobRL ./ repmat(sum(postprobRL, 1), [noPositionBins, 1]); %% normalized posterior probabilty matrix for each direction
    PosteriorProbMatrix{evt, 2} = postprobLR ./ repmat(sum(postprobLR, 1), [noPositionBins, 1]); 

  
    PosteriorProbMatrix{evt, 1}(:, idx) = eps*eps; 
    PosteriorProbMatrix{evt, 2}(:, idx) = eps*eps;
    
    
    [replayScore(evt, 2), slope(evt, 2), rho(evt, 2), begPosition(evt,:, 2), endPosition(evt,:, 2)] = replayevaluation(PosteriorProbMatrix{evt, 2}, 1);
    [replayScore(evt, 1), slope(evt, 1), rho(evt, 1), begPosition(evt,:, 1), endPosition(evt,:, 1)] = replayevaluation(PosteriorProbMatrix{evt, 1}, 1);
    
    
end


% -----------------



close all

figure;
set(gcf, 'Units', 'centimeters', 'position', [0 0 30 30])


t = 0;
for pbe = 1001:1100
    
    t = t+1;
    subplot(10,10, t)

    hold on

    currPPM = PosteriorProbMatrix{pbe, 1};

    maxProb = max(currPPM(:));
    minProb = min(currPPM(:));

    Clim_low = minProb;
    Clim_high = minProb + 0.75 *(maxProb - minProb);

    %%% plot the posterior probability matrix
    
    imagesc(currPPM, [Clim_low Clim_high])
    set(gca,'YDir','normal')
    colormap('jet')

    hold on

    %%%% plot the best fitting line
%         plot([begPosition_nullRL_cc(pbe, 1) endPosition_nullRL_cc(pbe, 1)], [begPosition_nullRL_cc(pbe, 2) endPosition_nullRL_cc(pbe, 2)],'color','w', 'linewidth', 2)

    
    plot([begPosition(pbe,1, 1) endPosition(pbe,1, 1)], [begPosition(pbe, 2, 1)+7 endPosition(pbe,2, 1)+7],'color','w', 'linestyle', '-.', 'linewidth', 0.25)
    plot([begPosition(pbe,1, 1) endPosition(pbe,1, 1)], [begPosition(pbe, 2, 1)-7 endPosition(pbe,2, 1)-7],'color','w', 'linestyle', '-.', 'linewidth', 0.25)

    title(sprintf('RS= %.2f', replayScore(pbe, 1)),'fontsize', 8, 'Units', 'normalized', 'Position', [0 1], 'HorizontalAlignment', 'left')
    set(gca, 'xtick',[],'ytick',[])

    xlim([0.5 size(currPPM, 2)+0.5])
    ylim([1 size(currPPM, 1)])
%         set(gca, 'visible', 'off')


end




%%% calculating the replay score for the shuffled data (with two
%%% shuffling methods)

for evt = 1100 %: noEvents
    
    
    currEvent = eventsBinnedfiring{evt, 2}(activeUnits, :); 

    notimeBins = size(currEvent, 2);
%     method = 1;
%     for method = 1 : noshufflingMethods
%         randShift = randi(noPositionBins-1, noShuffle, notimeBins);

 -

    replayScore_nullRL_cc = zeros(noShuffle, 1);
    replayScore_nullLR_cc = zeros(noShuffle, 1);

    %--------------------------------------------
    begPosition_nullRL_cc = zeros(noShuffle, 2);
    endPosition_nullRL_cc = zeros(noShuffle, 2);
    
    begPosition_nullRL_ts = zeros(noShuffle, 2);
    endPosition_nullRL_ts = zeros(noShuffle, 2);
    
    %--------------------------------------------
    
    
    replayScore_nullRL_ts = zeros(noShuffle, 1);
    replayScore_nullLR_ts = zeros(noShuffle, 1);

    PosteriorProbMatrixRL = PosteriorProbMatrix(evt, 1);
    PosteriorProbMatrixLR = PosteriorProbMatrix(evt, 2);
    
    %------------------------------------
    PosteriorProbMatrixRL_shuffle = cell(noShuffle, 1);
    PosteriorProbMatrixLR_shuffle = cell(noShuffle, 1);
    %------------------------------------
    
    
    for sn = 1 : noShuffle 


        % method one: column cycle shuffle
        randShift = randi(noPositionBins-1, 1, notimeBins);

%         PosteriorProbMatrixRL_shuffle{sn} = column_cycle_shuffle(PosteriorProbMatrixRL, randShift);
%         PosteriorProbMatrixLR_shuffle{sn} = column_cycle_shuffle(PosteriorProbMatrixLR, randShift);
        
%         [replayScore_nullRL_cc(sn), ~, ~, begPosition_nullRL_cc(sn, :), endPosition_nullRL_cc(sn, :)] = replayevaluation(PosteriorProbMatrixRL_shuffle{sn}, 1);
        
        
%         replayScore_nullRL_cc(sn) = replayevaluation(PosteriorProbMatrixRL_shuffle{sn}, 1);
%         replayScore_nullLR_cc(sn) = replayevaluation(PosteriorProbMatrixLR_shuffle{sn}, 1);

% 
%         % method 2: time swap
% 
%         PosteriorProbMatrixRL_shuffle = genTimeSwap(PosteriorProbMatrixRL);
%         PosteriorProbMatrixRL_shuffle = PosteriorProbMatrixRL_shuffle{1};
          

%-------------------------------

PosteriorProbMatrixRL_shuffle_temp = genTimeSwap(PosteriorProbMatrixRL);
PosteriorProbMatrixRL_shuffle{sn} = PosteriorProbMatrixRL_shuffle_temp{1};

%-------------------------------


%         PosteriorProbMatrixLR_shuffle = genTimeSwap(PosteriorProbMatrixLR);
%         PosteriorProbMatrixLR_shuffle = PosteriorProbMatrixLR_shuffle{1};
% 
%         replayScore_nullRL_ts(sn) = replayevaluation(PosteriorProbMatrixRL_shuffle, 1);
%         replayScore_nullLR_ts(sn) = replayevaluation(PosteriorProbMatrixLR_shuffle, 1);

[replayScore_nullRL_ts(sn), ~, ~, begPosition_nullRL_ts(sn, :), endPosition_nullRL_ts(sn, :)] = replayevaluation(PosteriorProbMatrixRL_shuffle{sn}, 1);


    end
    
 % --------------------------------------------------------- 
 
 
    figure;
    set(gcf, 'Units', 'centimeters', 'position', [0 0 30 30])
    
    for sn = 1:noShuffle
        
        subplot(10,10, sn)
        
        hold on
        
        currPPM = PosteriorProbMatrixRL_shuffle{sn};

        maxProb = max(currPPM(:));
        minProb = min(currPPM(:));

        Clim_low = minProb;
        Clim_high = minProb + 0.75 *(maxProb - minProb);

        %%% plot the posterior probability matrix
        imagesc(currPPM, [Clim_low Clim_high])
        set(gca,'YDir','normal')
        colormap('jet')

        hold on

        %%%% plot the best fitting line
%         plot([begPosition_nullRL_cc(sn, 1) endPosition_nullRL_cc(sn, 1)], [begPosition_nullRL_cc(sn, 2) endPosition_nullRL_cc(sn, 2)],'color','w', 'linewidth', 2)
        
        plot([begPosition_nullRL_ts(sn, 1) endPosition_nullRL_ts(sn, 1)], [begPosition_nullRL_ts(sn, 2)+7 endPosition_nullRL_ts(sn, 2)+7],'color','w', 'linestyle', '-.', 'linewidth', 0.25)
        plot([begPosition_nullRL_ts(sn, 1) endPosition_nullRL_ts(sn, 1)], [begPosition_nullRL_ts(sn, 2)-7 endPosition_nullRL_ts(sn, 2)-7],'color','w', 'linestyle', '-.', 'linewidth', 0.25)
        
        
        title(sprintf('RS= %.2f', replayScore_nullRL_ts(sn)),'fontsize', 8, 'Units', 'normalized', 'Position', [0 1], 'HorizontalAlignment', 'left')
        set(gca, 'xtick',[],'ytick',[])
        
        xlim([0.5 size(currPPM, 2)+0.5])
        ylim([0 size(currPPM, 1)])
%         set(gca, 'visible', 'off')
        
        
    end
    
    saveas(gcf, 'ExampleShufflePBEs', 'epsc')

    
    figure;
%     set(gca, 'Units', 'centimeters', 'position', [0 0 3 3])
    

    currPPM = PosteriorProbMatrixRL{:};

    maxProb = max(currPPM(:));
    minProb = min(currPPM(:));

    Clim_low = minProb;
    Clim_high = minProb + 0.75 *(maxProb - minProb);

    %%% plot the posterior probability matrix
    imagesc(currPPM, [Clim_low Clim_high])
    set(gca,'YDir','normal')
    colormap('jet')

    hold on

    %%%% plot the best fitting line
%         plot([begPosition_nullRL_cc(sn, 1) endPosition_nullRL_cc(sn, 1)], [begPosition_nullRL_cc(sn, 2) endPosition_nullRL_cc(sn, 2)],'color','w', 'linewidth', 2)

    plot([begPosition(evt,1, 1) endPosition(evt, 1, 1)], [begPosition(evt, 2, 1)+7 endPosition(evt, 2, 1)+7],'color','w', 'linestyle', '-.', 'linewidth', 0.25)
    plot([begPosition(evt,1, 1) endPosition(evt, 1, 1)], [begPosition(evt, 2, 1)-7 endPosition(evt, 2, 1)-7],'color','w', 'linestyle', '-.', 'linewidth', 0.25)


    title(sprintf('RS= %.2f', replayScore_nullRL_cc(sn)),'fontsize', 8, 'Units', 'normalized', 'Position', [0 1], 'HorizontalAlignment', 'left')
    set(gca, 'xtick',[],'ytick',[])

    xlim([0 size(currPPM, 2)])
    ylim([0 size(currPPM, 1)])
    
    
    
    

    replayScore_nullRL(:, 1) = replayScore_nullRL_cc;
    replayScore_nullRL(:, 2) = replayScore_nullRL_ts;

    replayScore_nullLR(:, 1) = replayScore_nullLR_cc;
    replayScore_nullLR(:, 2) = replayScore_nullLR_ts;

    %%% replayScore(z) or zrScore is defined as how distanct the
    %%% replay score of raw data is from median of the shuffle
    %%% distribution normalized to the distance of 95 percentile of the
    %%% distribution to the median. We do the same for each direction independently as well

    for method = 1:2

        BDprctile(evt, method, 1) = length(find(replayScore_nullRL(:,1) < replayScore(evt, 1)))/noShuffle;
        BDprctile(evt, method, 2) = length(find(replayScore_nullLR(:, 1) < replayScore(evt, 2)))/noShuffle;
    end

%         for ii = 1:2 %% direction
%             
%             BDprctile(evt, method, ii) = length(find(replayScore_null{ii}(evt,:,method) < replayScore(evt, ii)))/noShuffle;
%             
% %             [BDprctile(evt, method, ii), pval(evt, method, ii)] = rawScoreCompare2null(replayScore(evt, ii), replayScore_null{ii}(evt,:,method)); 
%         end
%         



    if mod(evt, 50) == 0
        fprintf('replay scores of event #%d are %.2f (R2L) and %.2f (L2R)\n', evt, BDprctile(evt,1,1), BDprctile(evt, 1, 2))
    end
    
end


subfolder = fullfile(FileBase, dataType);
mkdir(subfolder)

save(fullfile(subfolder, 'BDresults'), 'BDprctile', 'PosteriorProbMatrix', 'replayScore', 'begPosition', 'endPosition', 'slope') % , 'pval' , 'replayScore_null'

end

