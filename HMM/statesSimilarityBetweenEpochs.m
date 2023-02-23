

nStatesPRE = size(lambdaPRE, 2);
nStatesRUN = size(lambdaRUN, 2);
nStatesPOST = size(lambdaPOST, 2);


% generate simulated data

nSimulations = 100;
nUnits = size(lambdaPRE, 1);


%% PRE
simulationsPRE = zeros(nUnits, nStatesPRE, nSimulations);

for ii = 1: nStatesPRE
    simulationsPRE(:, ii, :) = permute(simulatePoisson(lambdaPRE(:,ii), nSimulations), [1 3 2]);
end


%% RUN
simulationsRUN = zeros(nUnits, nStatesRUN, nSimulations);

for ii = 1: nStatesRUN
    simulationsRUN(:, ii, :) = permute(simulatePoisson(lambdaRUN(:,ii), nSimulations), [1 3 2]);
end

%% POST
simulationsPOST = zeros(nUnits, nStatesPOST, nSimulations);

for ii = 1: nStatesPOST
    simulationsPOST(:, ii, :) = permute(simulatePoisson(lambdaPOST(:,ii), nSimulations), [1 3 2]);
end



%%

RUNgivenPRE_instances = zeros(nStatesPRE, nStatesRUN);
for ii = 1: nSimulations
    RUNgivenPRE_instances(:,:,ii) = poisson_prob(simulationsRUN(:,:,ii), lambdaPRE, 1);
end

RUNgivenPRE = mean(RUNgivenPRE_instances, 3);



POSTgivenRUN_instances = zeros(nStatesRUN, nStatesPOST);
for ii = 1: nSimulations
    POSTgivenRUN_instances(:,:,ii) = poisson_prob(simulationsPOST(:,:,ii), lambdaRUN, 1);
end

POSTgivenRUN = mean(POSTgivenRUN_instances, 3);


POSTgivenPRE_instances = zeros(nStatesPRE, nStatesPOST);
for ii = 1: nSimulations
    POSTgivenPRE_instances(:,:,ii) = poisson_prob(simulationsPOST(:,:,ii), lambdaPRE, 1);
end

POSTgivenPRE = mean(POSTgivenPRE_instances, 3);



%% finding the similartity between the states through inner products



sim_PRE_RUN = innersimilarity(lambdaPRE, lambdaRUN);

sim_PRE_POST = innersimilarity(lambdaPRE, lambdaPOST);

sim_RUN_POST = innersimilarity(lambdaRUN, lambdaPOST);




%% Functions

function simulations = simulatePoisson(meanFirings, nSimulations)
        
nUnits = length(meanFirings);
simulations = zeros(nUnits, nSimulations);

for ii = 1: nUnits
    
    simulations(ii, :) = poissrnd(meanFirings(ii), 1, nSimulations);
    
end

end


function simMat = innersimilarity(lambdaA, lambdaB)


simMat = zeros(size(lambdaA, 2), size(lambdaB, 2));

for ii = 1: size(lambdaA, 2)
    for jj = 1: size(lambdaB, 2)
        simMat(ii, jj) = (lambdaA(:, ii))' * lambdaB(:, jj)/ norm(lambdaA(:, ii))/ norm(lambdaB(:, jj));      
    end
end

end