function [transmat_departure_sparsity, lambda_state_sparsity, lambda_unit_sparsity] = modelSparsity(events, activeUnits, numofStates, FileBase, numIter)


% This function is intended to calculate the sparisty for the model
% (transition matrix and Lambda matrix).
% It also measures the degree of connectivity within a graph corresponding
% to the transition matrix, with the states and their transitions as the 
% graph's nodes and edges. In other words, the function quantifies the 
% degree to which there is a speicific pattern within the graph,
% which could correspond to dominant sequential transitions between states.
% The network connectivity is measured by: mixing time, bottleneck ratio,
% and also path analysis (specifically, the longest sequence of states 
% connected together in a thresholded graph in which weaker connections are
% omitted). 

numofUnits = length(activeUnits);

transmats = zeros(numofStates, numofStates, numIter);
lambdas = zeros(numofUnits, numofStates, numIter);
priors = zeros(numofStates, numIter);

transmat_departure_sparsity = zeros(numIter, numofStates);
transmat_arrival_sparsity = zeros(numIter, numofStates);

lambda_state_sparsity = zeros(numIter, numofStates);
lambda_unit_sparsity = zeros(numIter, length(activeUnits));

% mixingTime = zeros(numIter, 1);
% bottleneckRatio = zeros(numIter, 1);
% 
% % path analysis
% thresh_range = 0.2:0.05:0.8;
% maxPathLen = zeros(numIter, length(thresh_range));

% Results for differnt initializations

for ii = 1:numIter
    
    
    if mod(ii, 10) == 0
        fprintf(1, 'Iteration %d has started\n', ii);
    end
    
        
    [curr_transmat, curr_lambda, curr_prior] = trainHMM(events, activeUnits, numofStates);
    
    [~, maxInd] = max(curr_prior);
    
    sortInd = sortStates(curr_transmat, maxInd);
    curr_transmat = curr_transmat(sortInd, sortInd);
    
    
    transmats(:, :, ii) = curr_transmat;
    lambdas(:, :, ii) = curr_lambda(:, sortInd);
    priors(:, ii) = curr_prior(sortInd);
    
    
    
    % Sparsity
    
    % (1)Transition Matrix
    transmat_departure_sparsity(ii, :) = ginicoeff(curr_transmat, 2);
    transmat_arrival_sparsity(ii, :) = ginicoeff(curr_transmat, 1);
    
    % (2)Lambda Matrix
    lambda_state_sparsity(ii, :) = ginicoeff(curr_lambda, 1); % the degree to which cofirings corresponding to individual states are sparse  
    lambda_unit_sparsity(ii,:) = ginicoeff(curr_lambda, 2); % the degree to which single units participate sparsly across the states
    
    
    % Connectivity
    
    % (1)Path Analysis
    
    % for making an adjacency matrix we need to consider only the
    % transition probabiities above a certain threshold. Since the results
    % may depend on the defined threholds, we use a range of thresholds.
    
%     
%     maxPathLen(ii, :) = pathLen_thresh(curr_transmat, thresh_range)
% 
%     % (2) Mixing time and Bottleneck ratio
%     
%     [mixingTime(ii), bottleneckRatio(ii)] = connectivity(curr_transmat);
%     
end

% save([FileBase '.mat'], 'transmats', 'lambdas', 'priors', 'transmat_departure_sparsity', ...
%     'transmat_arrival_sparsity', 'lambda_state_sparsity', 'lambda_unit_sparsity', 'thresh_range', 'maxPathLen', 'mixingTime', 'bottleneckRatio')

save([FileBase '.mat'], 'transmats', 'lambdas', 'priors', 'transmat_departure_sparsity', ...
    'transmat_arrival_sparsity', 'lambda_state_sparsity', 'lambda_unit_sparsity')


end
