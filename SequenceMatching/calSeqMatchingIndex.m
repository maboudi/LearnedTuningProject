function [seqMatchingIdx, PBEsequence] = calSeqMatchingIndex(eventsBinnedfiring, template_LR, template_RL, activeUnits)

%%% inputs 

% PREbinnedPBEs: population burst events binned with 20 ms or 1 ms time
% bins

% template_LR: right to left run sequence template
% template_RL: left to right run sequence template

% activeUnits: in general, this enables selecting a subset of units instead of all

%%% outputs

% seqMatchingIdx: sequence similarity index based on the ratio of unit
% pairs firing in the same order as in the run template

% seqMatchingProb: matching probability characterizing the extent to which the PBE
% sequence is unlikely. This is calulated by comparing the sequence
% similarity index with a null distribution resulted from random
% permutation of the PBE sequence. This is like doing unit ID shuffle. 

if isempty(activeUnits)
   activeUnits = 1:size(eventsBinnedfiring{1}, 1); 
end


nEvents = size(eventsBinnedfiring, 1);
PBEsequence = cell(nEvents, 1);


if ~isempty(template_RL) % linear tracks

    seqMatchingIdx = struct('LR', [], 'RL', []);
    [seqMatchingIdx.LR, seqMatchingIdx.RL] = deal(struct('index', zeros(nEvents, 1), 'probability', zeros(nEvents, 1)));

    for pbe = 417%1:nEvents

        currEvent   = eventsBinnedfiring{pbe, 1}(activeUnits, :); % one milisecond time bins
        PBEsequence{pbe} = extractPBESequence(currEvent);

        % for the actual PBE find the matching index between the PBE sequence
        % (frame sequence as in JiWillsonNatNeurosci2007) and each of the
        % (run) templates; leftright (template_LR) nad rightleft (template_RL)

        [seqMatchingIdx.LR.index(pbe), seqMatchingIdx.LR.probability(pbe)] = sequenceSimilarity(PBEsequence{pbe}, template_LR);
        [seqMatchingIdx.RL.index(pbe), seqMatchingIdx.RL.probability(pbe)] = sequenceSimilarity(PBEsequence{pbe}, template_RL);


        if mod(pbe, 50) == 0 
            fprintf('\nPBEs %d to %d were processed: mean seqMatchIdx and seqMatchProb were %.2f and %.2f', pbe-49, pbe, ...
                nanmean([seqMatchingIdx.LR.index(pbe-49:pbe); seqMatchingIdx.RL.index(pbe-49:pbe)]), ...
                nanmean([seqMatchingIdx.LR.probability(pbe-49:pbe); seqMatchingIdx.RL.probability(pbe-49:pbe)]))
        end

    end
    
else % circular tracks
    
    
    seqMatchingIdx = deal(struct('index', zeros(nEvents, 1), 'probability', zeros(nEvents, 1)));

    for pbe = 1:nEvents

        currEvent   = eventsBinnedfiring{pbe, 1}(activeUnits, :); % one milisecond time bins
        PBEsequence{pbe} = extractPBESequence(currEvent);

        % for the actual PBE find the matching index between the PBE sequence
        % (frame sequence as in JiWillsonNatNeurosci2007) and each of the
        % (run) templates; leftright (template_LR) nad rightleft (template_RL)

        [seqMatchingIdx.index(pbe), seqMatchingIdx.probability(pbe)] = sequenceSimilarity(PBEsequence{pbe}, template_LR);
        
        if mod(pbe, 50) == 0 
            fprintf('\nPBEs %d to %d were processed: mean seqMatchIdx and seqMatchProb were %.2f and %.2f', pbe-49, pbe, ...
                nanmean(seqMatchingIdx.index(pbe-49:pbe)), ...
                nanmean(seqMatchingIdx.probability(pbe-49:pbe)))
        end

    end
end

end