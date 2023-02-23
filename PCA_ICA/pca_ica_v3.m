function [ICs, ICweights, PCeigenVals, eigValueThresh] = pca_ica_v3(data, varargin)

% parameters

analysis = 'all';

for ii = 1:2:length(varargin)-1 

    switch lower (varargin{ii})
        case 'ncomp'
            nComp = varargin{ii+1};
        case 'firingmean'
            firingMean = varargin{ii+1};
        case 'firingstd'
            firingStd  = varargin{ii+1};
        case 'onlypca'
            
            if varargin{ii+1} == 1
                analysis = 'white'; % only PCA
            else
                analysis = 'all'; % PCA and ICA
            end
    end
end


[nUnits, nTimeBins] = size(data);

% z-scoring in order to remove the bias due to differences in the average firing rates

if exist('firingMean', 'var')
    data = (data - repmat(firingMean, [1 nTimeBins]))./repmat(firingStd, [1 nTimeBins]);
else  
    data = zscore(data, [], 2); 
end
data(isnan(data)) = 0; % in case if we have any non-firing unit



eigValueThresh = (1+sqrt(nUnits/nTimeBins))^2; % theoretical upperbound of eigen value for iid variables 

if exist('nComp', 'var')
    [ICs, ~, ICweights, PCeigenVals] = fastica(data, 'verbose', 'off', 'lastEig', nComp, 'only', analysis); % lastEig specifies the index of the smallest eigenvalue or principal component
else
    [ICs, ~, ICweights, PCeigenVals] = fastica(data, 'verbose', 'off', 'eigvaluethresh', eigValueThresh, 'only', analysis); % I need to eliminate the extra whitening in this code
end

end
