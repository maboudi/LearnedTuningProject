



if ~isempty(activeUnits) 
    eventsBinnedfiring = justActiveUnits(eventsBinnedfiring, activeUnits);
    tuningRL = tuningRL(activeUnits, :);
    tuningLR = tuningLR(activeUnits, :);
end


trainDataAmount = floor(0.2*length(eventsBinnedfiring));

trainData = eventsBinnedfiring(1:trainDataAmount, 2);
testData = eventsBinnedfiring(trainDataAmount+1 : length(eventsBinnedfiring), 2);



numofStates = 60;  %% predetermined number of states

prior0 = normalise(rand(numofStates,1));
transmat0 = mk_stochastic(rand(numofStates,numofStates));
lambda0 = rand(length(activeUnits), numofStates);

[~, ~, transmat, lambda] = phmm_em(trainData, prior0, transmat0, lambda0, 'max_iter', 30);



loglike = zeros(length(testData), 1);
noBins = zeros(length(testData), 1);
conScore = zeros(length(testData), 1);

for evt = 1:length(testData)
    currEvt = testData{evt};
    noBins(evt) = size(currEvt, 2);
    
    B = poisson_prob(currEvt, lambda,1);
    
    conScore(evt) = log(prod(mean(B,1)));
    prior = 1/numofStates * ones(numofStates,1); %%% a uniform prior probability distribution over the states

    %%% observation likelihood for the raw data
    [~, ~, ~, loglike(evt) , ~] = fwdback(prior, transmat, B);
end



NumChunks = 5;
NumEvtPerChunk = floor(length(testData)/NumChunks);
LL = zeros(NumChunks, NumEvtPerChunk);
NBins = zeros(NumChunks, NumEvtPerChunk);

distLL = zeros(NumChunks, 3);
distNBins = zeros(NumChunks, 3);

for i = 1 : NumChunks
    
    LL(i,:) = loglike((i-1)*NumEvtPerChunk+1 : i*NumEvtPerChunk);
    NBins(i,:) = noBins((i-1)*NumEvtPerChunk+1 : i*NumEvtPerChunk);
    conSS(i, :) = conScore((i-1)*NumEvtPerChunk+1 : i*NumEvtPerChunk);
    
    distLL(i, :) = [median(LL(i,:)) prctile(LL(i,:), 25) prctile(LL(i,:), 75)];
    distNBins(i, :) = [median(NBins(i,:)) prctile(NBins(i,:), 25) prctile(NBins(i,:), 75)];
    distConSS(i, :) = [median(conSS(i,:)) prctile(conSS(i,:), 25) prctile(conSS(i,:), 75)]; 
    
    
end



figure;

subplot(1,3,1)
plot(distLL(:,1), 'r', 'linewidth', 3)
hold on
plot(distLL(:,2), '--r')
hold on
plot(distLL(:,3), '--r')

set(gca, 'fontsize', 16)
xlabel('Division of Data', 'fontsize', 16)
ylabel('Log Likelihood', 'fontsize', 16)
ylim([min(distLL(:)) max(distLL(:))])

subplot(1,3,2)
plot(distNBins(:,1), 'k', 'linewidth', 3)
hold on
plot(distNBins(:,2), '--k')
hold on
plot(distNBins(:,3), '--k')

set(gca, 'fontsize', 16)
ylabel('Number of Bins', 'fontsize', 16)
xlabel('Division of Data', 'fontsize', 16)

subplot(1,3,3)
plot(distConSS(:,1), 'k', 'linewidth', 3)
hold on
plot(distConSS(:,2), '--k')
hold on
plot(distConSS(:,3), '--k')

set(gca, 'fontsize', 16)
ylabel('Context Score', 'fontsize', 16)
xlabel('Division of Data', 'fontsize', 16)

    
    