function [A_B_corresp_pval,  sortAstates, gammagivenA, gammagivenB, obsProbA, obsProbB] = stateCorrespondence(testData, lastPBEoffirstPeriod, lambdaA, transmatA, lambdaB, transmatB, uniformTranprobs, seqContentCorrespondence)


if ~isempty(lastPBEoffirstPeriod) % for generating the pooled time swap data, I generate the shuffled data separately for each period and then consider them together
    
    testData1 = testData(1:lastPBEoffirstPeriod, :);
    testData2 = testData(lastPBEoffirstPeriod+1 : end, :);
end


noShuffles = 100;

% Decoding the states in the test period given model A and B (trained on the periods A and B, resp.)

% noTimebins = size(testData, 2);

% Decoded states given model A
numofStatesA = size(lambdaA, 2);

if ~uniformTranprobs
    currTransmatA = transmatA;
else
    currTransmatA = 1/numofStatesA * ones(numofStatesA);% to make the transition probabilities irrelevant 
end


priorA = 1/numofStatesA * ones(numofStatesA,1);


nPBEs = size(testData, 1);

obsProbA = cell(size(testData));
gammagivenA = cell(size(testData));
for pbe = 1: nPBEs
    obsProbA{pbe} =  poisson_prob(testData{pbe}, lambdaA, 1);
    [~, ~, gammagivenA{pbe}] = fwdback(priorA, currTransmatA, obsProbA{pbe});
end

cnctgammagivenA = cell2mat(gammagivenA');
cnctobsProbgivenA = cell2mat(obsProbA');

noTimebins = size(cnctgammagivenA, 2);


% Decoded states given model B
numofStatesB = size(lambdaB, 2);

if ~uniformTranprobs
    currTransmatB = transmatB;
else
    currTransmatB = 1/numofStatesB * ones(numofStatesB);
end


priorB = 1/numofStatesB * ones(numofStatesB,1);

obsProbB = cell(size(testData));
gammagivenB = cell(size(testData));
for pbe = 1: nPBEs
    obsProbB{pbe} =  poisson_prob(testData{pbe}, lambdaB, 1);
    [~, ~, gammagivenB{pbe}] = fwdback(priorB, currTransmatB, obsProbB{pbe});
end

cnctgammagivenB = cell2mat(gammagivenB');
cnctobsProbgivenB = cell2mat(obsProbB');


AstatesProbDis_chance = zeros(numofStatesA, noTimebins, noShuffles);
BstatesProbDis_chance = zeros(numofStatesB, noTimebins, noShuffles);

switch seqContentCorrespondence
    
    case 0
        
        AstatesProbDis = cnctobsProbgivenA ./ repmat(sum(cnctobsProbgivenA, 1), [numofStatesA 1]);
        BstatesProbDis = cnctobsProbgivenB ./ repmat(sum(cnctobsProbgivenB, 1), [numofStatesB 1]);
        
        for ii = 1: noShuffles
            AstatesProbDis_chance(:,:, ii) = AstatesProbDis(:, randperm(noTimebins));
            BstatesProbDis_chance(:,:, ii) = BstatesProbDis(:, randperm(noTimebins));
        end
    
    case 1 
        
        AstatesProbDis = cnctgammagivenA;
        BstatesProbDis = cnctgammagivenB;
        
        for ii = 1: noShuffles
            AstatesProbDis_chance(:,:, ii) = cnctgammagivenA(:, randperm(noTimebins));
            BstatesProbDis_chance(:,:, ii) = cnctgammagivenB(:, randperm(noTimebins));
        end
    
    case 2
        
        AstatesProbDis = cnctgammagivenA;
        BstatesProbDis = cnctgammagivenB;
        
        % obs probs given A
        obsProb_A1 = cell(size(testData1));
        for pbe = 1: size(testData1, 1) % we calculate the observation probabilities only once and for the shuffle data we shuffle the obs probs instead of the time-binned firing rates
            obsProb_A1{pbe} =  poisson_prob(testData1{pbe}, lambdaA, 1);
        end

        obsProb_A2 = cell(size(testData2));
        for pbe = 1:size(testData2, 1)
            obsProb_A2{pbe} = poisson_prob(testData2{pbe}, lambdaA, 1);
        end

        % obs probs given B
        obsProb_B1 = cell(size(testData1));
        for pbe = 1: size(testData1, 1) 
            obsProb_B1{pbe} =  poisson_prob(testData1{pbe}, lambdaB, 1);
        end

        obsProb_B2 = cell(size(testData2));
        for pbe = 1:size(testData2, 1)
            obsProb_B2{pbe} = poisson_prob(testData2{pbe}, lambdaB, 1);
        end



        for ii = 1: noShuffles  % it's like pooled time swap
            ii
            [obsProb_A1_ts, shuffleOrder1] = genTimeSwap(obsProb_A1);
            [obsProb_A2_ts, shuffleOrder2] = genTimeSwap(obsProb_A2);
            obsProb_A_ts = [obsProb_A1_ts; obsProb_A2_ts];

            obsProb_B1_ts = genTimeSwap(obsProb_B1, shuffleOrder1);
            obsProb_B2_ts = genTimeSwap(obsProb_B2, shuffleOrder2);
            obsProb_B_ts = [obsProb_B1_ts; obsProb_B2_ts];


            % given model A
            nPBEs = size(obsProb_A_ts, 1);

            gammagivenA_chance = cell(size(obsProb_A_ts));
            for pbe = 1: nPBEs
                [~, ~, gammagivenA_chance{pbe}] = fwdback(priorA, currTransmatA, obsProb_A_ts{pbe});
            end
            AstatesProbDis_chance(:,:, ii) = cell2mat(gammagivenA_chance');


            % given model B
            nPBEs = size(obsProb_B_ts, 1);

            gammagivenB_chance = cell(size(obsProb_B_ts));
            for pbe = 1: nPBEs
                [~, ~, gammagivenB_chance{pbe}] = fwdback(priorB, currTransmatB, obsProb_B_ts{pbe});
            end
            BstatesProbDis_chance(:,:, ii) = cell2mat(gammagivenB_chance');
        end
        
    case 3
        
        AstatesProbDis = cnctgammagivenA;
        BstatesProbDis = cnctgammagivenB;
        
        shuffledTmatA = genshuffles(currTransmatA, noShuffles, []);
        shuffledTmatB = genshuffles(currTransmatB, noShuffles, []);
        for ii = 1: noShuffles  % transition matrix shuffle
            ii
            % gamma for shuffled data given model A
            gammagivenA_chance = cell(size(testData));
            for pbe = 1: nPBEs
                obsProb =  poisson_prob(testData{pbe}, lambdaA, 1);
%                 [~, ~, gammagivenA_chance{pbe}] = fwdback(priorA, shuffledTmatA(:,:,randi(noShuffles)), obsProb);
                [~, ~, gammagivenA_chance{pbe}] = fwdback(priorA, shuffledTmatA, obsProb);
            end
            AstatesProbDis_chance(:,:, ii) = cell2mat(gammagivenA_chance');

            % gamma for shuffled data given model B
            gammagivenB_chance = cell(size(testData));
            for pbe = 1: nPBEs
                obsProb =  poisson_prob(testData{pbe}, lambdaB, 1);
%                 [~, ~, gammagivenB_chance{pbe}] = fwdback(priorB, shuffledTmatB(:,:,randi(noShuffles)), obsProb);
                [~, ~, gammagivenB_chance{pbe}] = fwdback(priorB, currTransmatB, obsProb);
            end
            BstatesProbDis_chance(:,:, ii) = cell2mat(gammagivenB_chance');

        end
end
    

A_B_corresp = (AstatesProbDis * BstatesProbDis')/noTimebins; % A states in rows and B states in columns

A_B_corresp_chance = zeros(numofStatesA, numofStatesB, noShuffles);

for ii = 1:noShuffles
    A_B_corresp_chance(:,:, ii) = (AstatesProbDis_chance(:,:, ii) * BstatesProbDis_chance(:,:, ii)')/noTimebins;
end


A_B_corresp_pval = zeros(numofStatesA, numofStatesB);
for istate = 1: numofStatesA
    for jstate = 1:numofStatesB
        
        [~, pval] = ttest2(A_B_corresp(istate, jstate), squeeze(A_B_corresp_chance(istate, jstate, :)), 'Tail', 'right', 'alpha', 0.01);
        A_B_corresp_pval(istate, jstate) = -log10(pval);
        
    end
end


% sorting the states in one period based on the maximum corresponding
% states from anoter period

[~, maxBstates] = max(A_B_corresp_pval, [], 2);

[~, sortAstates] = sort(maxBstates);

end

