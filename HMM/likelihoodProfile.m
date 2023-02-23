function [partialLikeli, bins_nullDistPercentiles] = likelihoodProfile(eventData, transmat, lambda, percentile, noshuffle)

%%% likelihoodProfile calculates the degree to which the observation in each
%%% time bin is significant given the model (or one may say how much the transition
%%% to the most probable state in that time bin is significant). To do this the data is
%%% truncated at a certain time bin. The likelihood of observing the
%%% truncated data given the model is calulated using forward recursion 
%%% algorithm (the forward parameter is calculated in each time bin using 
%%% the one for the previous time bin and the observation likelihood of the 
%%% data is calculated by marginalizing the parameter over the states). 
%%% Then, for measuring the significance of observation (or transition to 
%%% the most probable corresponding state) of the current time bin, we 
%%% preserved the model for estimating the forward parameter till the bin 
%%% prior to the certain time bin and shuffled just the transition matrix 
%%% in calcualting the forwardparameter corresponding to the last time bin. 


noTimeBins = size(eventData, 2);

numofStates = size(transmat, 1);



B = poisson_prob(eventData, lambda); 

prior = 1/numofStates * ones(numofStates,1); %% a uniform prior probability of the states


partialLikeli = zeros(1, noTimeBins);
partialLikeli_null = zeros(noshuffle, noTimeBins);

for bin = 1 : noTimeBins
    
    [~, ~, ~,  partialLikeli(bin), ~] = fwdback(prior, transmat, B(:, 1:bin), 'fwd_only', 1); %% we need just the forward part to calculate the observation likelihood
    
    
    for sn = 1 : noshuffle
        
        %%%% generating the shuffled transtion matrix 

        shuffled_transmat = zeros(size(transmat));

        for s = 1 : numofStates

           randgen = randperm(numofStates);
           randgen(randgen == s) = [];

           shuffleInd= [randgen(1:s-1)'; s; randgen(s : end)']; %% redistributing the transition probabilities in each column except the self-transitions which is preserved
           shuffled_transmat(:,s) = transmat(shuffleInd, s);

        end 
        
        %%% here the forwrad-backward algorithm is modified in a way to to consider the shuffle transition matrix for the last time bin in the recursion algorithm
        [~, ~, ~,  partialLikeli_null(sn, bin), ~] = fwdback2(prior, transmat, shuffled_transmat, B(:, 1:bin), 'fwd_only', 1); 
        
    end
    
end

bins_nullDistPercentiles = [prctile(partialLikeli_null, percentile)' median(partialLikeli_null)' prctile(partialLikeli_null, 100 - percentile)' prctile(partialLikeli_null, 75)']; %% a row corresponding to each bin
    
end
