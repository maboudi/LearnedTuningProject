function [HMMprctile, gammaDist, binResolvedScores] = HMMcongruence_crossvalid(data, testEvts, transmat, lambda, nShuffle, dataType, fileBase, shuffleMethod)
% [HMMScoreZ, binsCongSignificance, maxCongruenceLength, totalCongruenceLength] = measureHMMcongruence(data, testEvts, transmat, lambda, nShuffle, dataType, fileinfo, whichPart)


nEvents = size(data, 1);
nFolds = length(testEvts);
nStates = size(transmat, 1);



dataLL = zeros(1, nEvents); %% log likelihood of the raw data
nullLL = zeros(nEvents, nShuffle); %% log likelihood of shuffle data

% HMMScoreZ = zeros(nEvents, 1); %% 2 if we add again time swap
HMMprctile = zeros(nEvents, 1);

gammaDist = cell(nEvents, 1);
binResolvedScores = cell(nEvents, 1);

% partialLikeli = cell(1, nEvents); 
% bins_nullDistPercentiles = cell(1, nEvents);
% binsCongSignificance = cell(1, nEvents);
% maxCongruenceLength = zeros(1, nEvents);
% totalCongruenceLength = zeros(1, nEvents);

for fold = 1 : nFolds
    
    curr_transmat = transmat(:, :, fold);
    curr_lambda = lambda(:, :, fold);
    
    %
%     arrivalSparsity = ginicoeff(curr_transmat, 1);
%     departureSparsity = ginicoeff(curr_transmat, 2);
%     
%     diffSparsity = abs(departureSparsity' - arrivalSparsity);
%     
%     setaside = find(diffSparsity > prctile(diffSparsity, 75)); % keep aside the states with diffsparisties in the 4th quartile
%     
%     setaside = randi(nStates, 1, length(setaside)); % if we are keeping aside some random states
    
    setaside = [];
    shuffleTransmats = genshuffles(curr_transmat, nShuffle, setaside);

    %
    
    
    concatEvts = cell2mat(data(testEvts{fold}, 2)');
    concatB = poisson_prob(concatEvts, curr_lambda,1);
    
    totNbins = size(concatEvts, 2);
    for evt = testEvts{fold}

        
        currEvent = data{evt, 2};
        nTimeBins = size(currEvent, 2);
%         nTimebins = size(currEvent, 2);
        
        B = poisson_prob(currEvent, curr_lambda,1);
        prior = 1/nStates * ones(nStates,1); %%% a uniform prior probability distribution over the states
        
        %%% observation likelihood for the raw data
        [~, ~, gammaDist{evt},  dataLL(evt), ~] = fwdback(prior, curr_transmat, B);
        
        
       shuffleB = zeros(size(B,1), size(B,2), nShuffle);
       
       for ii = 1:nShuffle
          
           shuffleB(:, :, ii) = B(:, randperm(nTimeBins)); 
           
       end
        
        %%%% null distribution of observation likelihhod
        %%%% shuffling the transmat matrix
        
%         shuffled_transmat = zeros(size(curr_transmat));
        v = zeros(1, nShuffle);
        
        if shuffleMethod == 1
            parfor sn = 1:nShuffle

    %             for s = 1 : nStates
    %                 
    %                 %%% redistribute the transition probabilities within each 
    %                 %%% row with an exception that the self-transitions (matrix diagonal elements) are not changed
    %                 %
    %                 exclude = [setaside s];
    %                 initialIdx = setdiff(1:nStates, exclude);
    %                 
    %                 shuffleSet = initialIdx(randperm(length(initialIdx)));
    %                 
    %                 shuffleInd = zeros(nStates, 1);
    %                 
    %                 shuffleInd(exclude) = exclude;
    %                 shuffleInd(initialIdx) = shuffleSet;
    %                 
    %                 
    % %                 randgen = randperm(nStates)';
    % %                 randgen(randgen == s) = [];
    % % 
    % %                 shuffleInd = [randgen(1:s-1); s; randgen(s : end)]; 
    % 
    %                 shuffled_transmat(s,:) = curr_transmat(s, shuffleInd); % shuffling within the rows
    %                 
    % %                 shuffled_transmat(:,s) = curr_transmat(shuffleInd, s); % shuffling within the columns
    %             end             

    %             [~, ~, ~, nullLL(evt,sn), ~] = fwdback(prior, shuffled_transmat, B); %%% B and Prior are the same as in the case of raw model
                [~, ~, ~, v(sn), ~] = fwdback(prior, shuffleTransmats(:, :, sn), B);
            end
        else
            parfor sn = 1:nShuffle
                
                [~, ~, ~, v(sn), ~] = fwdback(prior, curr_transmat, shuffleB(:, :, sn));
            end
            
        end
        
        nullLL(evt, :) = v;
        
%         
%         %%%%%%% shuffling the events (time swap)
% HMMprctile
%         for sn = 1 : nShuffle
%             
%             %%% generating the time swapped data 
%             
%             rndgen = randperm(size(currEvent,2));
%             currEvent_sh = currEvent(:, rndgen); 
%     
% 
%             
%             B = poisson_prob(currEvent_sh, curr_lambda,1);
% 
%             [~, ~, ~, nullLL(evt,sn,2), ~] = fwdback(prior, curr_transmat, B);
% 
%         end
        
        
        for method = 1%:2

%             HMMScoreZ(evt, method) = rawScoreCompare2null(dataLL(evt), nullLL(evt,:,method)); 
            HMMprctile(evt, method) = length(find(nullLL(evt,:,method) < dataLL(evt)))/nShuffle * 100;

        end
        
        
        
%         % time-resolved analysis
%         
%         binResolvedScores{evt} = zeros(1, nTimeBins);
%         
%         for bin = 1: nTimeBins
%             
%             B_singleBinShuffle = B;
%             
%             tempShuffLL = zeros(nShuffle, 1);
%             for jj = 1:nShuffle
%               
% %                 B_singleBinShuffle(:, bin) = concatB(:, randi(totNbins));
%                 B_singleBinShuffle(:, bin) = B(:, randi(nTimeBins));
% 
%                 [~, ~, ~, tempShuffLL(jj), ~] = fwdback(prior, curr_transmat, B_singleBinShuffle);
%                 
%             end
%             
%             binResolvedScores{evt}(bin) = length(find(tempShuffLL < dataLL(evt)))/nShuffle;
%             
%         end
        
        
        
        if mod(evt, 100) == 0
            fprintf(1, 'cross-validation_%s_event %d, HMMPercentile = %f\n', dataType, evt, HMMprctile(evt, 1));
        end
%         
%         percentile = 5;
%         
%         [partialLikeli{evt}, bins_nullDistPercentiles{evt}] = likelihoodProfile(currEvent, curr_transmat, curr_lambda, percentile, nShuffle);
%         
%         binsCongSignificance{evt} = (partialLikeli{evt}' - bins_nullDistPercentiles{evt}(:,2))./(bins_nullDistPercentiles{evt}(:,3) - bins_nullDistPercentiles{evt}(:,2)); 
%         
%         [maxCongruenceLength(evt), totalCongruenceLength(evt)] = measureCongLen(partialLikeli{evt}, bins_nullDistPercentiles{evt});



    end

end

% save(fullfile(fileBase, [dataType '.mat']), 'HMMprctile', 'dataLL', 'nullLL')

% save([FileBase '/' dataType '.mat'], 'HMMScoreZ', 'dataLL', 'nullLL', 'partialLikeli', 'bins_nullDistPercentiles', 'binsCongSignificance', 'maxCongruenceLength', 'totalCongruenceLength')


end

function shuffleTransmats = genshuffles(transmat, nShuffle, setaside)


nStates = size(transmat, 1);
shuffleTransmats = zeros(nStates, nStates, nShuffle);

for sn = 1 : nShuffle
    
    
    shuffle_transmat = zeros(nStates, nStates);
    for s = 1 : nStates

        %%% redistribute the transition probabilities within each 
        %%% row with an exception that the self-transitions (matrix diagonal elements) are not changed
        %
        
        exclude = s;
        
%         exclude = [setaside s];
        initialIdx = setdiff(1:nStates, exclude);

        shuffleSet = initialIdx(randperm(length(initialIdx)));

        shuffleInd = zeros(nStates, 1);

        shuffleInd(exclude) = exclude; % keep the same indices for states within the exclude set
        shuffleInd(initialIdx) = shuffleSet; % use the randomized indices for the remained


%                 randgen = randperm(nStates)';
%                 randgen(randgen == s) = [];
% 
%                 shuffleInd = [randgen(1:s-1); s; randgen(s : end)]; 

%         shuffleTransmats(s,:, sn) = transmat(s, shuffleInd); % shuffling within the rows

        shuffle_transmat(:, s) = transmat(shuffleInd, s); % shuffling within the columns
    end
    
    % if doing column-wise shuffle
    
    shuffleTransmats(:, :, sn) = shuffle_transmat./repmat(sum(shuffle_transmat, 2), [1, nStates]);
%       shuffleTransmats(:, :, sn) = shuffle_transmat;
    
end

end