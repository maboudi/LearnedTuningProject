function [gamma_avg, gamma] = HMM_StatesPlaceField_2D(runBinnedfiring, transmat, lambda, posbinIdx, numofStates, xposBins, yposBins)
% [gamma_avg, gamma] = HMM_StatesPlaceField_2D(runBinnedfiring, transmat, lambda, posbinIdx, activeUnits, numofStates, xposBins, yposBins, normalise_lsPFs)

% This function uses the HMM trained on PBEs to decode the states in time bins of run
% data. 



% Decode the states in RUN data
% runData = justActiveUnits(runBinnedfiring, activeUnits);

% B = poisson_prob(runData{2}, lambda,1);

B = poisson_prob(runBinnedfiring, lambda,1);

prior = 1/numofStates * ones(numofStates,1); % a uniform prior 

% [~, ~, gamma,  runData_likelihood] = fwdback(prior, transmat, B); 
[~, ~, gamma] = fwdback(prior, transmat, B); 
% Note that the gamma icludes the probability of each state within each time bin  

noxPosBins = length(xposBins);
noyPosBins = length(yposBins);


gamma_avg = zeros(noyPosBins, noxPosBins, numofStates); % doing a summation over the state probabilty distributions corresponding to each position bin


for ii  = 1 : noyPosBins
    for jj = 1: noxPosBins
        
        idx = find(posbinIdx(:,2) == ii & posbinIdx(:,1) == jj);
        
        if ~isempty(idx)
            gamma_avg(ii, jj, :) = permute(sum(gamma(:,idx), 2)/length(idx), [3,2,1]); % the first dimension is the gamma matrix needs to be the 3rd in the avg_gamma(as a function of position matrix), hence [3,2,1]
        end
    
    end
end


% smooth the gamma
sigma = 3; % the same smoothing parameters used for place fields calculation of units 
halfwidth = 3 * sigma;
smoothwin = gausswindow(sigma, halfwidth);
smoothwin2D = smoothwin' * smoothwin;


for ii = 1 : numofStates
    gamma_avg(:, :, ii) = conv2(gamma_avg(:, :, ii), smoothwin2D, 'same');
end


% if  normalise_lsPFs
    gamma_avg = gamma_avg ./ repmat(sum(sum(gamma_avg, 1),2), [size(gamma_avg, 1) size(gamma_avg, 2)]); % normalize for each state
% end

% 
% % plot the states' place fields
% 
% figure('visible', 'off');
% 
% x0=0;
% y0=0;
% width=400;
% height=400* numofStates/30;
% set(gcf,'units','points','position',[x0,y0,width,height])
% 
% 
% for jj = 1 : numofStates
% 
%     fill([posbincenters fliplr(posbincenters)], [0.05*jj+gamma_sorted(jj, :) fliplr(0.05*jj*ones(size(gamma_sorted(jj, :))))], 'k','LineStyle','none')
%     hold on
% 
%     plot(posbincenters, 0.05*jj+gamma_sorted(jj, :),'color', 'k','linewidth', 0.1);
% 
%     alpha(0.25)
% 
%     set(gca, 'YTick', [], 'YTickLabel', [], 'color', 'none', 'YColor', 'none', 'box', 'off')
% 
% end
% 
% 
% difpos = posbincenters(2)- posbincenters(1);
% xlim([posbincenters(1)-difpos/2 posbincenters(end)+difpos/2])
% 
% xlabel('Position(cm)', 'fontsize', 10)
% %     text(posbincenters(floor(noPosBins/2)), 0, 'Position(cm)', 'fontsize', 10, 'HorizontalAlignment', 'center')
% 
% h = text(posbincenters(1)-3*difpos, 0.05*numofStates/2, 'State', 'fontsize', 10, 'HorizontalAlignment', 'center');
% set(h, 'rotation', 90)
% 
% pos = get(gcf,'Position');
% set(gcf,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[pos(3), pos(4)])
% print(gcf,[FileBase '/stateSpaceField_numStates' num2str(numofStates)],'-dpdf','-r0')


end
