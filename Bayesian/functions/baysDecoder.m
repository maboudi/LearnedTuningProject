% function postPr = baysDecoder( binnedFiring, tuning, binDur)
% 
% % exclude the times bins without any cell firing
% 
% % sumBinnedFiring = sum(binnedFiring, 1);
% % idx = find(sumBinnedFiring);
% 
% 
% noPosBins = size(tuning, 2); 
% nTimeBins = size(binnedFiring, 2);
% % noFiringTimeBins = length(idx);
% 
% postPr = zeros(noPosBins, nTimeBins);
% 
% 
% for pos = 1 : noPosBins
%     
%     
%     tuneTerm   = binDur*tuning(:,pos);
%     expTerm    = exp(-binDur * sum(tuning(:,pos)));
%     
% 
%     aa = zeros(nTimeBins, 1);
%     for time = 1 : nTimeBins  
%         
%         firingCounts = binnedFiring(:, time);
%         try
%             aa(time) = prod(tuneTerm.^ firingCounts) * expTerm ; %/ prod(factorial(firingCounts(firingCounts > 0))) : we can think about a look up table if we are interested to include the factorial terms
%         catch
%             aa(time) = prod(tuneTerm.^ firingCounts) * expTerm;
%         end
%     
%     end
%     postPr(pos,:) = aa;
% end


%%%%%%%%%%%%%%%%%%%
% a little better %
%%%%%%%%%%%%%%%%%%%

% function postPr = baysDecoder( binnedFiring, tuning, binDur)
% 
% 
% % The difference with the version above is that we calculate the log
% % probabilities (logProbs) and at the end with revert it back to praobility values by
% % calcualting exponentials og the logProb values
% 
% % logProb = -\tau*(\lambda1 + \lambda2 + ... + \lambdaN) + n1*log(\tau * \lambda1) + n2*log(\tau * \lambda2) + ... + nN*log(\tau * \lambdaN)
% % logProb = firstTerm + secondterm;
% % prob    = exp(logProb);
% 
% 
% noPosBins = size(tuning, 2); 
% nTimeBins = size(binnedFiring, 2);
% 
% 
% logTunings      = log(binDur * tuning);
% tuningSummation = - binDur * summation(tuning, 1);
% 
% 
% postPr = zeros(noPosBins, nTimeBins);
% 
% 
% for pos = 1 : noPosBins
%     
%     firstTerm  = tuningSummation(:, pos);
%     logLambdas = logTunings(:, pos);
%     
%     
%     
%     if nTimeBins > 1000 % in cases other than individual PBEs, for example continuous decoding during a sleep period
%         
%         secondTerm = zeros(nTimeBins, 1);
%         for time = 1 : nTimeBins  
% 
%             firingCounts = binnedFiring(:, time);
% 
%             secondTerm(time) =  firingCounts' * logLambdas;
% 
%         end
%         postPr(pos,:) = exp(firstTerm + secondTerm);
% 
%     else % for individual PBEs
%         
%         secondTerm    = sum(firingCounts .* repmat(logLambdas, [1 nTimeBins], 1);
%         postPr(pos,:) = exp(firstTerm + secondTerm);
%     end
%     
% end
% 
% end



function postPr = baysDecoder(binnedFiring, tuning, binDur)


% The difference with the version above is that we calculate the log
% probabilities (logProbs) and at the end with revert it back to praobility values by
% calcualting exponentials og the logProb values

% logProb = -\tau*(\lambda1 + \lambda2 + ... + \lambdaN) + n1*log(\tau * \lambda1) + n2*log(\tau * \lambda2) + ... + nN*log(\tau * \lambdaN)
% logProb = firstTerm + secondterm;
% prob    = exp(logProb);


noPosBins = size(tuning, 2); 
nTimeBins = size(binnedFiring, 2);


logTunings      = log(binDur * tuning);
tuningSummation = - binDur * sum(tuning, 1);



if nTimeBins < 1000 % this is suitbale if we are decoding bins within individual PBEs

    logPr = sum(repmat(permute(logTunings, [1 3 2]), [1 nTimeBins 1]) .* repmat(binnedFiring, [1 1 noPosBins]), 1) + ...
         repmat(permute(tuningSummation, [1 3 2]), [1 nTimeBins 1]);

    logPr  = permute(logPr, [3 2 1]);
    postPr = exp(logPr);
    
else
    
    
    maxBins   = 1000;
    nSegments = ceil(nTimeBins/maxBins);
    
    postPr = nan(noPosBins, nTimeBins);
    
    for ii = 1:nSegments
        
        currTimeBins  = (ii- 1)*maxBins+1 : min(ii* maxBins, nTimeBins); 
        nCurrTimeBins = numel(currTimeBins);
        
        logPr = sum(repmat(permute(logTunings, [1 3 2]), [1 nCurrTimeBins 1]) .* repmat(binnedFiring(:, currTimeBins), [1 1 noPosBins]), 1) + ...
         repmat(permute(tuningSummation, [1 3 2]), [1 nCurrTimeBins 1]);

        logPr  = permute(logPr, [3 2 1]);
        postPr(:, currTimeBins) = exp(logPr);
        
    end
    
end


end

