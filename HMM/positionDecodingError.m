function [PosDecError, PosDecError_shuffledlsPFs, actualPositions, decodedPositions_actual, probabilityOverPosition] = positionDecodingError(runBinnedfiring, transmat, lambda, posbinIdx, xposcenters, yposcenters, numofStates)
% [PosDecError, PosDecError_shuffledOrder, PosDecError_shuffledlsPFs] = positionDecodingError(runBinnedfiring, transmat, lambda, posbinIdx, xposcenters, yposcenters, numofStates)


excludeBins = find(posbinIdx(:,1).* posbinIdx(:,2) == 0);

posbinIdx(excludeBins, :) = [];
runBinnedfiring(:, excludeBins) = [];


noTimeBins = length(posbinIdx);

% 5-fold cross-validation

noFolds = 5;
foldSize = floor(noTimeBins/ noFolds);

decodedPositions_actual = zeros(length(posbinIdx), 2);
% decodedPositions_shuffledOrder = zeros(length(posbinIdx), 2);
decodedPositions_shuffledlsPFs = zeros(length(posbinIdx), 2);

% probabilityOverPosition = zeros(length(yposcenters), length(xposcenters), length(posbinIdx));


for k = 1:noFolds
    
    
    if k == noFolds
        testSet = (k-1)*foldSize+1 : noTimeBins;
    else
        testSet = (k-1)*foldSize+1 : k*foldSize;
    end
    
    trainSet = setdiff(1:noTimeBins, testSet);
    
    
    %% calculate lsPFs for the train set
    
    
    train_timeBins = runBinnedfiring(:,trainSet);
    
    train_posbinIdx = posbinIdx(trainSet, :); % x and y coordiantes as columns
    
    train_lsPFs = HMM_StatesPlaceField_2D(train_timeBins, transmat, lambda, train_posbinIdx, numofStates, xposcenters, yposcenters);
    
    
%     train_shuffle_posbinIdx = train_posbinIdx(randperm(size(trainSet, 1), size(trainSet, 1)), :);
    train_shuffle_posbinIdx = train_posbinIdx(randperm(length(trainSet)), :); % the occupancy will stay the same
    
    shuffle_train_lsPFs = HMM_StatesPlaceField_2D(train_timeBins, transmat, lambda, train_shuffle_posbinIdx, numofStates, xposcenters, yposcenters);
    
    
    % Decode the states in RUN data 
    
    test_timeBins = runBinnedfiring(:,testSet);
    
    obsProb = poisson_prob(test_timeBins, lambda,1);
      
    prior = 1/numofStates * ones(numofStates,1); % a uniform prior
    
    
    %(1) Decoding track positions based on actual lsPFs
    
    decodedPositions_actual(testSet, :) = lsPFPositionDecode(obsProb, transmat, prior, train_lsPFs, xposcenters,yposcenters, 'actual'); 
    
%     probabilityOverPosition(:,:,testSet) = rr;
    %--- shuffled time bin orders
    
%     decodedPositions_shuffledOrder(testSet, :) = lsPFPositionDecode(obsProb, transmat, prior, train_lsPFs, xposcenters,yposcenters, 'shuffledOrder'); 
    
    
    %--- shuffled gamma: it supposed to take away all the order and
    %coactivity information regarding the states
    
%     decodedPositions_shuffledlsPFs(testSet, :) = lsPFPositionDecode(obsProb, transmat, prior, train_lsPFs, xposcenters,yposcenters, 'shuffledlsPFs'); 
    

    %(2) Decoding track positions based on shuffled lsPFs (lsPFs calculated after shuffling the position indices)
    
    decodedPositions_shuffledlsPFs(testSet, :) = lsPFPositionDecode(obsProb, transmat, prior, shuffle_train_lsPFs, xposcenters,yposcenters, 'actual'); 

end


posbinIdx(find(posbinIdx(:,2) > length(yposcenters)), 2) = length(yposcenters);
posbinIdx(find(posbinIdx(:,1) > length(xposcenters)), 1) = length(xposcenters);


actualPositions = [yposcenters(posbinIdx(:,2))' xposcenters(posbinIdx(:,1))'];


PosDecError = sqrt((decodedPositions_actual(:, 1) - actualPositions(:, 1)).^2  + (decodedPositions_actual(:, 2) - actualPositions(:, 2)).^2);
% PosDecError(isnan(PosDecError)) = [];


% PosDecError_shuffledOrder = sqrt((decodedPositions_shuffledOrder(:, 1) - yposcenters(posbinIdx(:,2))').^2  + (decodedPositions_shuffledOrder(:, 2) - xposcenters(posbinIdx(:,1))').^2);

PosDecError_shuffledlsPFs = sqrt((decodedPositions_shuffledlsPFs(:, 1) - yposcenters(posbinIdx(:,2))').^2  + (decodedPositions_shuffledlsPFs(:, 2) - xposcenters(posbinIdx(:,1))').^2);    
% PosDecError_shuffledlsPFs(isnan(PosDecError_shuffledlsPFs)) = [];



end
