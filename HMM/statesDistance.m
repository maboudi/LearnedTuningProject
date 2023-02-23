function [DistanceWeighMean, distances_probThresholded] = statesDistance(transmat, corrPositions, varargin)
%%% This function calculate the distance between positions corresponding to
%%% each pair of states.  
%%% Applying a threshold on the transition probabilities, it finds a
%%% distribution of distances between the involving states in those
%%% transitions 

%%% It calculates mean of the distances wighthed by the transition
%%% probabilities

if nargin > 2
   transProbThresh = varargin{1};
else
   transProbThresh = prctile(transmat(:), 95); 
end


%% defining the states distance matrix

numofStates = size(transmat, 1);

distmat = zeros(size(transmat));

for ii = 1 : numofStates
    for jj = 1 : numofStates
        distmat(ii, jj) = abs(corrPositions(ii) - corrPositions(jj));
    end
end


%% calculating the mean wieghted distance

transmat = transmat .* (ones(numofStates) - eye(numofStates)); %% removing the slef-transitions

weighted_distmat = distmat .* transmat;
DistanceWeighMean = sum(weighted_distmat(:)) / sum(transmat(:));


%% thresholding

highProbTransitions = find(transmat > transProbThresh);
distances_probThresholded = distmat(highProbTransitions); %%% ./ transmat(highProbTransitions)

end

