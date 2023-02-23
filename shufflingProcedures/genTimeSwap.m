function [outPBEs, shuffleOrder] = genTimeSwap(inPBEs, varargin)

noPBEs = size(inPBEs, 1);


if iscell(inPBEs) % for pooled time swap
    cnctPBEs = cell2mat(inPBEs');
else
    cnctPBEs = inPBEs; % just one PBE, corresponding to within-PBE time swap
end

    

if nargin > 1
    shuffleOrder = varargin{1};
else
    shuffleOrder = randperm(size(cnctPBEs, 2));
end

Tswap_cnctPBEs = cnctPBEs(:, shuffleOrder);


if iscell(inPBEs)

    outPBEs = cell(size(inPBEs));
    for pbe = 1: noPBEs

        nTimeBins = size(inPBEs{pbe}, 2);

        outPBEs{pbe} = Tswap_cnctPBEs(:, 1:nTimeBins);
        Tswap_cnctPBEs(:, 1:nTimeBins) = [];

    end
else
    outPBEs = Tswap_cnctPBEs;
end

end
