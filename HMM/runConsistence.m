function [runLikeli, runZScore, runLikeli_s] = runConsistence(runBinnedfiring, transmat, lambda, numTransmatShuffle)

numofThetaBouts = length(runBinnedfiring);
subsize = numofThetaBouts;

runLikeli = zeros(subsize, 1); %%% we just score subsize bouts per run session
runLikeli_s = zeros(subsize, numTransmatShuffle); %%% scores from transmat shuffling
runZScore = zeros(subsize, 1);

numofStates = size(transmat, 1);
t = 0;

for bout = 1:subsize
    
    t = t+1
    
    currLap = runBinnedfiring{bout};

    B = poisson_prob(currLap, lambda,1);
    prior = 1/numofStates * ones(numofStates,1); %%% a uniform prior probability distribution over the states
    
    
    %%% observation likelihood for the real data
    [~, ~, ~,  runLikeli(t), ~] = fwdback(prior, transmat, B);
    
    
    
    shuffled_transmat = zeros(size(transmat));
    
    for shuffle = 1 : numTransmatShuffle
        
        for s = 1 : numofStates
            
            %%% redistribute the transition probabilities within each
            %%% column with an exception that the self-transitions (matrix diagonal elements) are not changed
            
            randgen = randperm(numofStates)';
            randgen(randgen == s) = [];
            shuffleInd = [randgen(1:s-1); s; randgen(s : end)];
            shuffled_transmat(:,s) = transmat(shuffleInd, s);
        end
        [~, ~, ~,  runLikeli_s(t, shuffle), ~] = fwdback(prior, shuffled_transmat, B);
        
    end
    
    runZScore(t) = (runLikeli(t) - mean(runLikeli_s(t, :)))/std(runLikeli_s(t,:));
end
