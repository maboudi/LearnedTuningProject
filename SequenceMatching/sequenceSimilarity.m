function [seqMatchingIdx, matchProb] = sequenceSimilarity(PBEsequence, template)

% inputs:
% PBEsequence: PBE sequence constructed by convoling the spikes trains with
% a gaussian kernel (sigma = 180 ms for hippocampal cells based on JiWilson,NatNeurosci2007)
% template: template sequence during run


% outputs
% seqMatchingIdx: sequence matching index
% matchProb: matching probability

[m, n] = size(PBEsequence); 
if m < n
   PBEsequence = transpose(PBEsequence);
end
PBEsequence(~ismember(PBEsequence, template)) = []; %PBE sequence elements not existing in the tempalte sequence are omitted. 
                                                    %including the non-template elements doesn't make a difference regarding the final results, but increases the processing time


PBEseqLen = numel(PBEsequence);

if PBEseqLen <= 4
    seqMatchingIdx = nan;
    matchProb      = nan;
else

[m, n] = size(template); 
if m < n
   template = transpose(template);
end


% data

seqMatchingIdx = measureSeqMatchingIdx(PBEsequence, template);


if nargout > 1

% shuffles

if PBEseqLen < 8
    
   allperms = perms(PBEsequence);
    
elseif PBEseqLen >= 8
   
    % generate 10000 random permutations
    
    allperms = zeros(10000, PBEseqLen);
    
    for ii = 1:10000
        allperms(ii, :) = PBEsequence(randperm(PBEseqLen));
    end
  
end


seqMatchingIdx_null = zeros(size(allperms, 1), 1);

parfor ii = 1:size(allperms, 1)    
    currSeq = allperms(ii, :)';
    
    if isequal(currSeq, PBEsequence) || isequal(currSeq, PBEsequence(end:-1:1))
       seqMatchingIdx_null(ii) = nan;
    else
       seqMatchingIdx_null(ii) = measureSeqMatchingIdx(currSeq, template);  
    end
end

matchProb = length(find(abs(seqMatchingIdx_null) >= abs(seqMatchingIdx)))/size(allperms, 1);

end

end
end


function seqMatchingIdx = measureSeqMatchingIdx(PBEsequence, template)


% PBEseqLen          = numel(PBEsequence);
PBEsequence_sorted = sortPBEseq(PBEsequence, template); 

% totalPairNum = PBEseqLen*(PBEseqLen-1)/2;


seqMatrix = repmat(transpose(PBEsequence_sorted), [length(PBEsequence_sorted) 1]);
diffMatrix = seqMatrix - transpose(seqMatrix);
diffMatrix = tril(diffMatrix, -1);

rankDifferences = diffMatrix(:);
rankDifferences(rankDifferences == 0) = []; 

consNumofPairs   = numel(find(rankDifferences < 0)); % number of pairs with order between them consistent with their order in the template sequence
inconsNumofpairs = numel(find(rankDifferences > 0)); % the same as above but for inconsistent pairs

seqMatchingIdx = (consNumofPairs - inconsNumofpairs)/(consNumofPairs + inconsNumofpairs);

end

function sequence_sorted = sortPBEseq(sequence, template)

% sequence(~ismember(sequence, template)) = []; %PBE sequence elements not existing in the tempalte sequence are omitted

sequence_sorted = zeros(size(sequence));

for ii = 1:numel(sequence)
    
    sequence_sorted(ii) = find(template == sequence(ii));
end
    
end

        