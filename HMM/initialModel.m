function [transmat0, lambda0, prior0] = initialModel(numofStates, noActiveUnits)

prior0 = normalise(rand(numofStates,1));
transmat0 = mk_stochastic(rand(numofStates,numofStates));
lambda0 = rand(noActiveUnits, numofStates);

end