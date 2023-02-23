function  [transmat, lambda, prior] = trainHMM(eventsBinnedfiring, activeUnits, numofStates)


if ~isempty(activeUnits) 
    eventsBinnedfiring = justActiveUnits(eventsBinnedfiring, activeUnits);
end


% initializing parameters of the model

prior0 = normalise(rand(numofStates,1));
transmat0 = mk_stochastic(rand(numofStates,numofStates));
lambda0 = rand(length(activeUnits), numofStates);

% training the model
[~, prior, transmat, lambda] = phmm_em(eventsBinnedfiring, prior0, transmat0, lambda0, 'max_iter', 200);


end

