
function [sortedTransmat, sortedLambda] = sortStates(transmat, lambda, prior)

% This function sort the transition probability matrix calcualted using HMM
% Starting from an initial state (with max prior probability), it
% calculates the most probable state the has not yet been visited

[~, startIdx] = max(prior);

[~, nStates] = size(lambda);
visited = startIdx;

available = setdiff(1:nStates, visited);
allVisited = 0;
nextState = startIdx;
while ~allVisited
   
   [~, idx] = max(transmat(nextState, available)); 
   nextState = available(idx);
   
   visited = [visited nextState];
   
   available = setdiff(1:nStates, visited);
   if isempty(available)
      allVisited = 1;
   end
end

sortInd = visited;

sortedTransmat = transmat(sortInd, sortInd);
sortedLambda = lambda(:, sortInd);


end