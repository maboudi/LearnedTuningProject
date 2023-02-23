
function shuffleTransmats = genshuffles(transmat, noShuffle, setaside)


numofStates = size(transmat, 1);
shuffleTransmats = zeros(numofStates, numofStates, noShuffle);

for sn = 1 : noShuffle
    
    
    shuffle_transmat = zeros(numofStates, numofStates);
    for s = 1 : numofStates

        %%% redistribute the transition probabilities within each 
        %%% row with an exception that the self-transitions (matrix diagonal elements) are not changed
        %
        
        exclude = s;
        
%         exclude = [setaside s];
        initialIdx = setdiff(1:numofStates, exclude);

        shuffleSet = initialIdx(randperm(length(initialIdx)));

        shuffleInd = zeros(numofStates, 1);

        shuffleInd(exclude) = exclude; % keep the same indices for states within the exclude set
        shuffleInd(initialIdx) = shuffleSet; % use the randomized indices for the remained


%                 randgen = randperm(numofStates)';
%                 randgen(randgen == s) = [];
% 
%                 shuffleInd = [randgen(1:s-1); s; randgen(s : end)]; 

%         shuffleTransmats(s,:, sn) = transmat(s, shuffleInd); % shuffling within the rows

%         shuffle_transmat(:, s) = transmat(shuffleInd, s); % shuffling within the columns
        
        shuffle_transmat(:, s) = transmat(shuffleInd, s);
    end
    
    % if doing column-wise shuffle
    
    shuffleTransmats(:, :, sn) = shuffle_transmat./repmat(sum(shuffle_transmat, 2), [1, numofStates]);
    
end

end
