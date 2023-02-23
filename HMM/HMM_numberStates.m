function  HMM_numberStates(eventsBinnedfiring, activeUnits, Kfold, maxNumStates, stepSize, FileBase)

% This function is intended to investigate the effect of pre-determined
% number of states in an HMM model on the model's performance. The number 
% of states in increased gradually in certain steps and performance of the 
% model is evaluated through K-fold cross-validation. 


noEvents = length(eventsBinnedfiring);

if ~isempty(activeUnits) 
    eventsBinnedfiring = justActiveUnits(eventsBinnedfiring, activeUnits);
end


randEvtIdx = randperm(noEvents); % randomize the event indices
segmentSize = floor(noEvents/ Kfold); % the size of the folds

numberofDims = floor(maxNumStates/stepSize);

trainingLogLike = zeros(numberofDims, Kfold);
trainingLogLike_n = zeros(numberofDims, Kfold);

validationLike = zeros(noEvents, numberofDims);
validationLike_n = zeros(noEvents, numberofDims);

inum = 0;
for numofStates = stepSize:stepSize:maxNumStates
    
    inum = inum+1;
    
    for k = 1 : Kfold
        
        validationSet = randEvtIdx((k-1)*segmentSize+1 : k*segmentSize);
        trainingSet = setdiff(randEvtIdx, validationSet);


        validationData = eventsBinnedfiring(validationSet, 2);
        trainingData = eventsBinnedfiring(trainingSet, 2);
        
        trainingBins = size(cell2mat(trainingData'), 2); % number of training bins, we need this for normalization
        
        
        % training an HMM on training PBEs
        
        prior0 = normalise(rand(numofStates,1));
        transmat0 = mk_stochastic(rand(numofStates,numofStates));
        lambda0 = rand(length(activeUnits), numofStates);
        
        [trainingLogLike(inum, k), ~, transmat, lambda] = phmm_em(trainingData, prior0, transmat0, lambda0, 'max_iter', 30, 'verbose', 0);
        
        trainingLogLike_n(inum, k) = trainingLogLike(inum, k)/trainingBins;
        
        for ievt = 1:length(validationSet)
            
            evt = validationSet(ievt);

            currEvt = validationData{ievt};

            B = poisson_prob(currEvt, lambda,1);
            prior = 1/numofStates * ones(numofStates,1); %%% a uniform prior probability distribution over the states

            %%% likelihood of the validation events
            [~, ~, ~, validationLike(evt, inum) , ~] = fwdback(prior, transmat, B);
            validationLike_n(evt, inum) = validationLike(evt, inum) / size(currEvt, 2);
            
        end
        
    end
    
end

save([FileBase '/stateSensitivity.mat'], 'trainingLogLike' , 'trainingLogLike_n', 'validationLike', 'validationLike_n')
    

% plot the results


validationLike(find(validationLike(:,1) == 0),:) = []; % remove the event which were not in any validation sets
validationLike_n(find(validationLike_n(:,1) == 0),:) = [];


figure('visible', 'off');
hold on

sem = std(validationLike_n)/sqrt(length(validationLike_n)); % standard error

upperBound = mean(validationLike_n) + sem;
lowerBound = mean(validationLike_n) - sem;

xAxisSet = stepSize:stepSize:maxNumStates;

fill([xAxisSet'; flipud(xAxisSet')], [lowerBound'; flipud(upperBound')], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none')


plot(xAxisSet, mean(trainingLogLike_n, 2), 'color', [0.7 0.7 0.7], 'linewidth', 2)


xlabel('number of states', 'fontsize', 20)
ylabel('LogLikelihood (normalized)', 'fontsize', 20)
title('Effect of number of states', 'fontsize', 20)
set(gca,'fontsize', 16)

saveas(gcf, [FileBase '/sensitivity2noStates.fig'])

end

% 
% set(gcf,'Units','Inches');
% pos = get(gcf,'Position');
% set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(gcf,'Sensitivity_noStates_pin01','-dpdf','-r0')


    


