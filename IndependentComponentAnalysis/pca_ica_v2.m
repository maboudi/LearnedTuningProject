function [ICs, ICweights, PCeigenVals, thresh] = pca_ica_v2(data, varargin)



[nUnits, nTimeBins] = size(data);

% parameters

if nargin == 4 
   cumEVthresh = varargin{1}; % threshold on cumulative PC eigenvalues (as a ratio of sum)
   firingMean  = varargin{2};
   firingStd   = varargin{3};
   
elseif nargin == 3 
   firingMean  = varargin{1};
   firingStd   = varargin{2};

elseif nargin == 2    
   cumEVthresh = varargin{1};
end


% z-scoring in order to remove the bias due to differences in the average firing rates

if exist('firingMean', 'var')
    data = (data - repmat(firingMean, [1 nTimeBins]))./repmat(firingStd, [1 nTimeBins]);
else  
    data = zscore(data, [], 2); 
end
data(isnan(data)) = 0; % in case if we have any non-firing unit

 
% [PCweights, PCs, PCeigenVals] = pca(data'); % the principal components are sorted based on their eigenvalues (highest to the lowest)
% [PCweights, PCs, PCeigenVals2] = pca(data', 'Algorithm', 'eig', 'Rows', 'pairwise');

% columns of input correspond to mixed variables/signals
% columns of PCs correspond to principal components
% Each column of PCweights contains coefficients for one principal
% component (eigen Vectors)

 
% if exist('cumEVthresh', 'var') 
%     nPCs   = find(cumsum(PCeigenVals) <= cumEVthresh*sum(PCeigenVals), '1', 'last');
%     thresh = PCeigenVals(nPCs);
% else
%     thresh = (1+sqrt(nUnits/nTimeBins))^2; % theoretical upperbound of eigen value for iid variables 
%     nPCs   = find(PCeigenVals > thresh, 1, 'last');
% %       nPCs = find(cumsum(PCeigenVals) > 0.9*sum(PCeigenVals), 1, 'first');
% end



thresh = (1+sqrt(nUnits/nTimeBins))^2; % theoretical upperbound of eigen value for iid variables 
% signPCs = PCs(:, 1:nPCs);

% [ICs, ~, unmixingMatrix] = fastica(signPCs', 'verbose', 'off'); % I need to eliminate the extra whitening in this code
[ICs, ~, ICweights, PCeigenVals] = fastica(data, 'verbose', 'off', 'eigvaluethresh', thresh); % I need to eliminate the extra whitening in this code
% , 'only', 'white'

% ICweights = PCweights(:, 1:nPCs)*unmixingMatrix;

% ICweights = ICweights';

end
