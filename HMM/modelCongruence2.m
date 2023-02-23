
function [HMMprctile, gammaDist, dataLL] = modelCongruence2(testData, transmat, lambda, nShuffle, fileBase, trainPeriod, testPeriod, ifShuffle)


% For generating pooled time swap data

cnctPBEs = cell2mat(testData(:,2)');

rndgen = randperm(size(cnctPBEs, 2)); 
shuffled_cnctPBEs = cnctPBEs(:, rndgen); % randomizing the order of time bins within the concatenated data  


nEvents = size(testData, 1);
nStates = size(transmat, 1);

dataLL = zeros(1, nEvents); %% log likelihood of the raw data
nullLL = zeros(nEvents, nShuffle); %% log likelihood of shuffle data

HMMprctile = zeros(nEvents, 1);


setaside = [];
shuffleTransmats = genshuffles(transmat, nShuffle, setaside);

gammaDist = cell(nEvents, 1);        


for evt = 1 : nEvents
    
    nTimebins = size(testData{evt, 2}, 2);
    
    if ifShuffle
        currEvent = shuffled_cnctPBEs(:, 1:nTimebins);
        shuffled_cnctPBEs(:, 1:nTimebins) = [];
    else
        currEvent = testData{evt, 2};
    end
    
      
    B = poisson_prob(currEvent, lambda,1);
    prior = 1/nStates * ones(nStates,1); %%% a uniform prior probability distribution over the states

    %%% observation likelihood for the actual data
    [~, ~, gammaDist{evt},  dataLL(evt), ~] = fwdback(prior, transmat, B);
    
    %%%% null distribution of observation likelihhod
    %%%% shuffling the transmat matrix

%     shuffled_transmat = zeros(size(transmat));
    v = zeros(1, nShuffle);
    parfor sn = 1 : nShuffle

%         for s = 1 : nStates
% 
%             %%% redistribute the transition probabilities within each 
%             %%% row with an exception that the self-transitions (matrix diagonal elements) are not changed
% 
%             randgen = randperm(nStates)';
%             randgen(randgen == s) = [];
% 
%             shuffleInd = [randgen(1:s-1); s; randgen(s : end)]; 
% 
% %             shuffled_transmat(s,:) = transmat(s, shuffleInd); 
%             
%             shuffled_transmat(:,s) = transmat(shuffleInd, s);
% 
%         end             

%         [~, ~, ~, nullLL(evt,sn), ~] = fwdback(prior, shuffled_transmat, B); %%% B and Prior are the same as in the case of raw model

%         [~, ~, ~, v(sn), ~] = fwdback(prior, shuffleTransmats(:,:,sn), B); %%% B and Prior are the same as in the case of raw model

          [~, ~, ~, v(sn), ~] = fwdback(prior, transmat, B(:, randperm(size(B, 2))));
        
    end
    nullLL(evt,:) = v;



    for method = 1%:2

        HMMprctile(evt, method) = length(find(nullLL(evt,:,method) < dataLL(evt)))/nShuffle *100;
    end
    
    
    if ifShuffle
        dataType = 'pooledTimeSwap';
    else
        dataType = 'actual';
    end
    
    
    if mod(evt, 100) == 0
        fprintf(1, [trainPeriod '-' testPeriod '-' dataType '__event %d, HMMPercentile = %f\n'], evt, HMMprctile(evt, 1));
    end
    
    
end


% 
% currDir = pwd;
% FileBase = [currDir '/' fileinfo.name  '/HMM/sequenceDetection/bw_periods_congruence/'];
% mkdir(FileBase)
% 
% mkdir(fullfile(fileBase, testPeriod))
% 
% save(fullfile(fileBase, testPeriod,  ['/Train_' trainPeriod '___Test_' testPeriod '__' dataType '.mat']), 'HMMprctile', 'dataLL', 'nullLL')



end