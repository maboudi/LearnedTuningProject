function [maxCongruenceLength, totalCongruenceLength]= measureCongLen(partialLikelihood, lastbinShuffledLikeliDist)


isbinSignificant = zeros(1,length(partialLikelihood));

isbinSignificant(find(partialLikelihood > lastbinShuffledLikeliDist(:,4)')) = 1; %%% the forth element is 75 percentile
isbinSignificant(find(partialLikelihood < lastbinShuffledLikeliDist(:,4)' | isnan(partialLikelihood))) = 0;
    

%%% find the significance epochs

significantBins = find(isbinSignificant);

if isempty(significantBins)
    
    maxCongruenceLength = 0;
    totalCongruenceLength = 0;
    
else
    sigEpochsBeg = significantBins(find([0 diff(significantBins)] > 1));
    sigEpochsBeg = [significantBins(1) sigEpochsBeg];

    sigEpochsEnd = significantBins(find([0 diff(significantBins)] > 1) - 1);
    sigEpochsEnd = [sigEpochsEnd significantBins(end)];

    sigEpochsDuration = sigEpochsEnd - sigEpochsBeg + 1;

    %%% remove one bin-length epochs
    epochs2remove = find(sigEpochsDuration == 1);

    sigEpochsBeg(epochs2remove) = [];
    sigEpochsEnd(epochs2remove) = [];
    sigEpochsDuration(epochs2remove) = [];


    %%% merge the epochs (tolerating between-epochs single non-significane bins)
    
    if isempty(sigEpochsDuration)
        
        maxCongruenceLength = 0;
        totalCongruenceLength = 0;
        
    else
        
        for epoch = 1 : length(sigEpochsBeg)-1

            if  sigEpochsBeg(epoch + 1) - sigEpochsEnd(epoch) - 1 == 1 %%% if the distance between the two epochs is just one bin merge them
                sigEpochsEnd(epoch) = 0;
                sigEpochsBeg(epoch + 1) = 0;
            end

        end
        
        sigEpochsDuration = sigEpochsEnd - sigEpochsBeg + 1;
        
        maxCongruenceLength = max(sigEpochsDuration);
        totalCongruenceLength = sum(sigEpochsDuration);
        
    end

end




