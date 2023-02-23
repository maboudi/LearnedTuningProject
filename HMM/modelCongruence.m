
function [HMMprctile, gammaDist, dataLL] = modelCongruence(testData, transmat, lambda, nShuffle, fileBase, trainPeriod, testPeriod, ifShuffle, shuffleMethod)


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
    
    nTimeBins = size(testData{evt, 2}, 2);
    
    if ifShuffle
        currEvent = shuffled_cnctPBEs(:, 1:nTimeBins);
        shuffled_cnctPBEs(:, 1:nTimeBins) = [];
    else
        currEvent = testData{evt, 2};
    end
    
      
    B = poisson_prob(currEvent, lambda,1);
    prior = 1/nStates * ones(nStates,1); %%% a uniform prior probability distribution over the states

    %%% observation likelihood for the actual data
    [~, ~, gammaDist{evt},  dataLL(evt), ~] = fwdback(prior, transmat, B);
    
    
    
    
    %%%% null distribution of observation likelihhod
    %%%% shuffling the transition matrix

%     shuffled_transmat = zeros(size(transmat));
    v = zeros(1, nShuffle);
    
    if shuffleMethod == 1 % shuffling the transition matrix
        
        parfor sn = 1:nShuffle
            [~, ~, ~, v(sn), ~] = fwdback(prior, shuffleTransmats(:,:,sn), B); %%% B and Prior are the same as in the case of raw model
        end
        
    else % within-PBE time swap
        
        shuffleB = zeros(size(B,1), size(B,2), nShuffle);
        for ii = 1:nShuffle
            shuffleB(:, :, ii) = B(:, randperm(nTimeBins)); 
        end
        
        parfor sn = 1:nShuffle
            [~, ~, ~, v(sn), ~] = fwdback(prior, transmat, shuffleB(:, :, sn));
        end
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
        fprintf(1, ['\n' trainPeriod '-' testPeriod '-' dataType '__event %d, HMMPercentile = %f'], evt, HMMprctile(evt, 1));
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