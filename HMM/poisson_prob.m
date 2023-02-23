function B = poisson_prob( data, lambda, factor )
% posisson_prob evaluate the pdf of a conditional multiple possion
% distribution

% lambda is the matrix of rate given the neuron and state
% data(:,t) represent the observation at each time
% factor determine whether the factorial term of Poisson pdf is considered in the
% calculation of the probabilities (by default is set to zero meaning the
% factorial term is ommitted)

if nargin < 3 
    factor = 0; 
end

[~, T] = size(data); %% N number of neurons and T number of time bins
[~, Q] = size(lambda);

ff = factorial(data);

logprob = repmat(-sum(lambda',2), [1, T]) + log(lambda)'* data - repmat(sum(log(ff), 1), [Q, 1]);
B = exp(logprob);

end


