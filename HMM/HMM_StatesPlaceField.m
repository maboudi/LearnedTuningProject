function HMM_StatesPlaceField(eventsBinnedfiring, runBinnedfiring, posbinIdx, activeUnits, numofStates, posbincenters, FileBase)
%     function [SpatialInfo_uniOccProb, SpatialInfo_actOccProb] = HMM_StatesPlaceField(eventsBinnedfiring, runBinnedfiring, posbinIdx, posbin_Occ, activeUnits, numofStates, posbincenters, FileBase) 
% This function uses the HMM trained on PBEs to decode the states in time bins of run
% data. 

% posbin_Occ = posbin_Occ./sum(posbin_Occ);

if isempty(activeUnits)
    activeUnits = 1:size(eventsBinnedfiring{1,2}, 1);
end

% Train an HMM on PBEs
[transmat, lambda] = trainHMM(eventsBinnedfiring(:,2), activeUnits, numofStates);


% Decode the states in RUN data
runData = justActiveUnits(runBinnedfiring, activeUnits);

B = poisson_prob(runData{2}, lambda,1);
prior = 1/numofStates * ones(numofStates,1); % a uniform prior 

% [~, ~, gamma,  runData_likelihood] = fwdback(prior, transmat, B); 
[~, ~, gamma] = fwdback(prior, transmat, B); 
% Note that the gamma icludes the probability of each state within each time bin  

noPosBins = length(posbincenters);
gamma_avg = zeros(numofStates, noPosBins); % doing a summation over the probabilty distributions corresponding to each position bin


for ii  = 1 : noPosBins
    idx = find(posbinIdx == ii);
    gamma_avg(:, ii) = sum(gamma(:,idx), 2)/length(idx);
end


% smooth the gamma
sigma = 2; % the same smoothing parameters used for place fields calculation of units 
halfwidth = 5 * sigma;
smoothwin = gausswindow(sigma, halfwidth);

peak = zeros(numofStates, 1); % calculating the peak for each states, which then will be used for sorting the states

for jj = 1 : numofStates
    gamma_avg(jj, :) = conv(gamma_avg(jj, :), smoothwin, 'same');
    [~, peak(jj)] = max(gamma_avg(jj, :));
end

[~, sortIdx] = sort(peak, 'ascend');

gamma_sorted = gamma_avg(sortIdx, :); % sort the states
gamma_sorted = gamma_sorted ./ repmat(sum(gamma_sorted, 2), [1 noPosBins]); % Normalize ecah row (belonging to a state) to one


%%% note that from now on we are dealing with the sorted states and the
%%% indices have changed

% calcualte the spatial Information, with and without considering
% occupancy probability
% 
% SpatialInfo_uniOccProb = zeros(numofStates, 1); % spatial entropy for each state (w/o considering occupancy probabilities)
% SpatialInfo_actOccProb = zeros(numofStates, 1);
% 
% for jj = 1 : numofStates
% 
%     SpatialInfo_uniOccProb(jj) = sum((1/noPosBins) * (gamma_sorted(jj, :)/mean(gamma_sorted(jj, :))) .* log2(gamma_sorted(jj, :)/mean(gamma_sorted(jj, :))));
% 
%     SpatialInfo_actOccProb(jj) = sum(posbin_Occ' .* (gamma_sorted(jj, :)/mean(gamma_sorted(jj, :))) .* log2(gamma_sorted(jj, :)/mean(gamma_sorted(jj, :))));
% 
%     
% end

% save([FileBase '/stateSpaceField_numStates' num2str(numofStates) '.mat'], 'transmat', 'lambda', ...
%     'gamma', 'gamma_avg', 'gamma_sorted', 'runData_likelihood', 'SpatialInfo_uniOccProb', 'SpatialInfo_actOccProb')

save([FileBase '/stateSpaceField_numStates' num2str(numofStates) '.mat'], 'transmat', 'lambda', ...
    'gamma', 'gamma_avg', 'gamma_sorted')

 
% plot the states' place fields

figure('visible', 'off');

x0=0;
y0=0;
width=400;
height=400* numofStates/30;
set(gcf,'units','points','position',[x0,y0,width,height])


for jj = 1 : numofStates

    fill([posbincenters fliplr(posbincenters)], [0.05*jj+gamma_sorted(jj, :) fliplr(0.05*jj*ones(size(gamma_sorted(jj, :))))], 'k','LineStyle','none')
    hold on

    plot(posbincenters, 0.05*jj+gamma_sorted(jj, :),'color', 'k','linewidth', 0.1);

    alpha(0.25)

    set(gca, 'YTick', [], 'YTickLabel', [], 'color', 'none', 'YColor', 'none', 'box', 'off')

end


difpos = posbincenters(2)- posbincenters(1);
xlim([posbincenters(1)-difpos/2 posbincenters(end)+difpos/2])

xlabel('Position(cm)', 'fontsize', 10)
%     text(posbincenters(floor(noPosBins/2)), 0, 'Position(cm)', 'fontsize', 10, 'HorizontalAlignment', 'center')

h = text(posbincenters(1)-3*difpos, 0.05*numofStates/2, 'State', 'fontsize', 10, 'HorizontalAlignment', 'center');
set(h, 'rotation', 90)

pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[pos(3), pos(4)])
print(gcf,[FileBase '/stateSpaceField_numStates' num2str(numofStates)],'-dpdf','-r0')


end
