function [ICs, ICweights, PCeigenVals, thresh, mu, sigma] = pca_ica(data, varargin)


[nUnits, nTimeBins] = size(data);

% data = (data - repmat(mean(data, 2), [1 nTimeBins]))./repmat(std(data, 0, 2), [1 nTimeBins]);

[data, mu, sigma] = zscore(data, [], 2); % in order to remove the bias due to differences in the average firing rates
data(isnan(data)) = 0; % in case if we have any non-firing unit

[PCweights, PCs, PCeigenVals] = pca(data'); % the principal components are sorted based on their eigenvalues
% columns of input correspond to mixed variables/signals
% columns of PCs correspond to principal components
% Each column of PCweights contains coefficients for one principal
% component (eigen Vectors)


   
if nargin > 1
    nPCs = varargin{1};
else
    thresh  = (1+sqrt(nUnits/nTimeBins))^2; % theoretical upperbound of eigen value for iid variables 
    nPCs = find(PCeigenVals > thresh, 1, 'last');
%       nPCs = find(cumsum(PCeigenVals) > 0.9*sum(PCeigenVals), 1, 'first');
end


% prejecting back again to the observation (mixed signal) space and doing ICA
% datahat = PCs(:, 1:nPCs) * PCweights(:, 1:nPCs)';
% [ICs, ~, ICweights] = fastica(datahat', 'verbose', 'off'); 

signPCs = PCs(:, 1:nPCs);


[ICs, ~, unmixingMatrix] = fastica(signPCs', 'verbose', 'off'); % I need to eliminate the extra whitening in this code

ICweights = PCweights(:, 1:nPCs)*unmixingMatrix;


ICweights = ICweights';

end
