function decodedPositions = lsPFPositionDecode(obsProb, transmat, prior, train_lsPFs, xposcenters,yposcenters, shuffleMode)



noTestBins = size(obsProb, 2);
% numofStates = size(transmat, 1);
shuffledOrder = randperm(noTestBins, noTestBins);

    
if strcmp(shuffleMode, 'shuffledOrder')
    [~, ~, gamma(:, shuffledOrder)] = fwdback(prior, transmat, obsProb(:, shuffledOrder));
else
    [~, ~, gamma] = fwdback(prior, transmat, obsProb);
end


decodedPositions = zeros(noTestBins, 2);
for jj = 1: noTestBins
    
        
%     if strcmp(shuffleMode, 'shuffledGamma')
%         stateProbability = gamma(randperm(numofStates, numofStates), jj);
%     else
%         stateProbability = gamma(:, jj);
%     end
    
    stateProbability = gamma(:, jj);
    
    probabilityOverPosition = sum(repmat(permute(stateProbability, [3 2 1]), length(yposcenters), length(xposcenters)) .* train_lsPFs, 3); % averaging the states' lsPFs weighted by the probability of the states given in the current time bin

    decodedPositions(jj, 1) = yposcenters(ceil(sum((1:length(yposcenters))' .* sum(probabilityOverPosition, 2))));
    decodedPositions(jj, 2) = xposcenters(ceil(sum((1:length(xposcenters)) .* sum(probabilityOverPosition, 1))));
end


end
