function   HMMScaling(eventsBinnedfiring, activeUnits, fileinfo, whichPart)


%%% folder for results

currDir = pwd;
FileBase = [currDir '/' fileinfo.name '/data part' num2str(whichPart) '/HMM/Scaling' ];
mkdir(FileBase)

noEvents = length(eventsBinnedfiring);

%%% preprocess the candidate events and tunings: including just the active units

if ~isempty(activeUnits) 
    eventsBinnedfiring = justActiveUnits(eventsBinnedfiring, activeUnits);
end


wholeTrainSize = floor(0.7 * noEvents); %% assigning the events randomely to train and test sets
randgen = randperm(noEvents);

WholeTrainData = eventsBinnedfiring(randgen(1:wholeTrainSize), 2);
testData = eventsBinnedfiring(randgen(wholeTrainSize+1 : noEvents), 2);



trainDataSizes = 25:25:wholeTrainSize;

randgen = randperm(wholeTrainSize); %% generating the random for once

fprintf('total number of training sets is %d\n', length(trainDataSizes))


%%% initialization of HMM
numofStates = 60; 

prior0 = normalise(rand(numofStates,1));
transmat0 = mk_stochastic(rand(numofStates,numofStates));
lambda0 = rand(length(activeUnits), numofStates);

transmat = cell(length(trainDataSizes), 1);
lambda = cell(length(trainDataSizes), 1);
distLL = zeros(3, length(trainDataSizes));
loglike = zeros(length(testData), length(trainDataSizes));

for ii = 1 : length(trainDataSizes)
    
    fprintf('\n\n\ntraining a model by %d randomely selected events ...\n', trainDataSizes(ii))
    
    currTrainSet = WholeTrainData(randgen(1:trainDataSizes(ii)));
    
    [~, ~, transmat{ii}, lambda{ii}] = phmm_em(currTrainSet, prior0, transmat0, lambda0, 'max_iter', 30);
    
    
    fprintf('\nevaluating the trained model on the validation set ...\n')
    

    for evt = 1:length(testData)
        
        currEvt = testData{evt};

        B = poisson_prob(currEvt, lambda{ii},1);
        prior = 1/numofStates * ones(numofStates,1); %%% a uniform prior probability distribution over the states

        %%% observation likelihood for the raw data
        [~, ~, ~, loglike(evt,ii) , ~] = fwdback(prior, transmat{ii}, B);
    end
    
    distLL(:,ii) = [median(loglike(:,ii)); prctile(loglike(:,ii), 25); prctile(loglike(:,ii), 75)];
    
    fprintf('the median likelihood of validation set is %.2f\n', median(loglike(:,ii)))
end

figure('visible', 'off');

plot(trainDataSizes, distLL(1, :), 'r', 'linewidth', 3)
hold on
plot(trainDataSizes, distLL(2, :), '--r')
hold on
plot(trainDataSizes, distLL(3, :), '--r')

set(gca, 'fontsize', 16)
xlabel('Training size', 'fontsize', 20)
ylabel('Log Likelihood', 'fontsize', 20)
ylim([min(distLL(:)) max(distLL(:))])

sigDiffp = zeros(length(trainDataSizes), 1);
for ii = 2 : length(trainDataSizes)
    [~, sigDiffp(ii)] = ttest(loglike(:,ii), loglike(:, ii-1));
end

saveas(gcf, [FileBase '/LLvstrainSize.fig'])
save([FileBase '/trainedModels.mat'], 'loglike', 'transmat', 'lambda', 'sigDiffp')


end
    
    


