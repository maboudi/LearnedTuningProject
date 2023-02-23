function HMMprctile = HMMcongruence_crossvalid_pts(data, testEvts, transmat, lambda, nShuffle, dataType, fileBase, shuffleMethod)


noEvents    = size(data, 1);
numofFolds  = length(testEvts);
numofStates = size(transmat, 1);



dataLL = zeros(1, noEvents); %% log likelihood of the raw data
nullLL = zeros(noEvents, nShuffle); %% log likelihood of shuffle data

% HMMScoreZ = zeros(noEvents, 1); %% 2 if we add again time swap
HMMprctile = zeros(noEvents, 1);

% partialLikeli = cell(1, noEvents); 
% bins_nullDistPercentiles = cell(1, noEvents);
% binsCongSignificance = cell(1, noEvents);
% maxCongruenceLength = zeros(1, noEvents);
% totalCongruenceLength = zeros(1, noEvents);

for fold = 1 : numofFolds
    
    curr_transmat = transmat(:, :, fold);
    curr_lambda   = lambda(:, :, fold);
    
    %
%     arrivalSparsity = ginicoeff(curr_transmat, 1);
%     departureSparsity = ginicoeff(curr_transmat, 2);
%     
%     diffSparsity = abs(departureSparsity' - arrivalSparsity);
%     
%     setaside = find(diffSparsity > prctile(diffSparsity, 75)); % keep aside the states with diffsparisties in the 4th quartile
%     
%     setaside = randi(numofStates, 1, length(setaside)); % if we are keeping aside some random states
    
    setaside = [];
    shuffleTransmats = genshuffles(curr_transmat, nShuffle, setaside);

    %
    
    % testing the model on shuffled test PBEs (pooled timeswap surrogate dataset)
    testPBEs = data(testEvts{fold}, :);
    
%     pts_testPBEs = cell(size(testPBEs));
%     pts_testPBEs(:,1) = []; % 1 ms bins, we don't need to compute it here
    
    cnctPBEs = cell2mat(testPBEs(:,2)');
    
    % resample each surrogate PBE (or better to say the bins) from the concatenated PBEs
    
    rndgen = randperm(size(cnctPBEs, 2));
    shuffled_cnctPBEs = cnctPBEs(:, rndgen); 
    
    for evt = testEvts{fold}

        
        currEvent = data{evt, 2};
        nTimeBins = size(currEvent, 2);
        
        surrogateEvent = shuffled_cnctPBEs(:, 1:nTimeBins);
        shuffled_cnctPBEs(:, 1:nTimeBins) = []; % the array should exhaust when we are done with all of the events

        B = poisson_prob(surrogateEvent, curr_lambda,1);
        prior = 1/numofStates * ones(numofStates,1); %%% a uniform prior probability distribution over the states
        
        %%% observation likelihood for the raw data
        [~, ~, ~,  dataLL(evt), ~] = fwdback(prior, curr_transmat, B);
        
        
       shuffleB = zeros(size(B,1), size(B,2), nShuffle);
       
       for ii = 1:nShuffle
          
           shuffleB(:, :, ii) = B(:, randperm(nTimeBins)); 
           
       end
        
        
        %%%% null distribution of observation likelihhod
        %%%% shuffling the transmat matrix
        
%         shuffled_transmat = zeros(size(curr_transmat));
        v = zeros(nShuffle, 1);
        
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
        
        
        if mod(evt, 100) == 0
            fprintf(1, 'cross-validation_%s_event %d, HMMPercentile = %f\n', dataType, evt, HMMprctile(evt, 1));
        end
        
        %%%% a time-resolved analysis
%         
%         percentile = 5;
%         
%         [partialLikeli{evt}, bins_nullDistPercentiles{evt}] = likelihoodProfile(currEvent, curr_transmat, curr_lambda, percentile, nShuffle);
%         
%         binsCongSignificance{evt} = (partialLikeli{evt}' - bins_nullDistPercentiles{evt}(:,2))./(bins_nullDistPercentiles{evt}(:,3) - bins_nullDistPercentiles{evt}(:,2)); 
%         
%         [maxCongruenceLength(evt), totalCongruenceLength(evt)] = measureCongLen(partialLikeli{evt}, bins_nullDistPercentiles{evt});
%         
    end

end

% save(fullfile(fileBase, [dataType '.mat']), 'HMMprctile', 'dataLL', 'nullLL')

% save([FileBase '/' dataType '.mat'], 'HMMScoreZ', 'dataLL', 'nullLL', 'partialLikeli', 'bins_nullDistPercentiles', 'binsCongSignificance', 'maxCongruenceLength', 'totalCongruenceLength')


end
