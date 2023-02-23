function varargout = HMMParametersVisualiztion(lambda, transmat, tuningRL, tuningLR, binDur, positionBins, FileBase)


noPositionBins = size(tuningRL, 2);
noActiveUnits = size(tuningRL, 1);
numofStates = size(transmat, 1);
%%%% the Proabability distribution of positions for each direction
%%%% separately; in this step we don't normalize the distributions to one
%%%% as we will be normalizing the distribution across the directions as
%%%% well

PositionsLikelihoodRL = baysDecoder(lambda, tuningRL, binDur);
PositionsLikelihoodLR = baysDecoder(lambda, tuningLR, binDur); 

DirectionalLikratio = sum(PositionsLikelihoodRL, 1) ./ sum(PositionsLikelihoodLR, 1); 
%%%  ratio > 1 if the state is most likley representing leftward place fields (sum of the likelihood over all positions, distribution were not normalized)
%%%  ratio < 1 if rightward place fields are being most likely represented by the state




PositionsLikelihoodMargin = (PositionsLikelihoodRL + PositionsLikelihoodLR)./repmat(sum(PositionsLikelihoodRL + PositionsLikelihoodLR, 1), [noPositionBins, 1]);
%%% the likelihhod distributions marginalized over the directions (summed over the directions and then normalized to sum of one)


%%% sorting the states based on the maximum likely posoitions corresponding
%%% to the states

[~, mlPositionMa] = max(PositionsLikelihoodMargin);

[~, statesSortIndMa] = sort(mlPositionMa, 'ascend'); %%% states sorting index based on the marginalized distribution

PositionsLikelihoodMargin = PositionsLikelihoodMargin(:, statesSortIndMa);

lambdaMa = lambda(:, statesSortIndMa); 
transmatMa = transmat(statesSortIndMa, statesSortIndMa);



%%%%% sorting based on directional distributions


%%% sorting based on RL

PositionsLikelihoodRL = PositionsLikelihoodRL ./ repmat(sum(PositionsLikelihoodRL, 1), [noPositionBins, 1]); %%% normalizing the distributions

[~, mlPositionRL] = max(PositionsLikelihoodRL);
[~, statesSortIndRL] = sort(mlPositionRL, 'ascend');


lambdaRL = lambda(:, statesSortIndRL);
transmatRL = transmat(statesSortIndRL, statesSortIndRL);


%%% sorting based on LR

PositionsLikelihoodLR = PositionsLikelihoodLR ./ repmat(sum(PositionsLikelihoodLR, 1), [noPositionBins, 1]); %%% normalizing the distributions

[~,mlPositionLR] = max(PositionsLikelihoodLR);  
[~,statesSortIndLR] = sort(mlPositionLR, 'ascend');

lambdaLR = lambda(:, statesSortIndLR);
transmatLR = transmat(statesSortIndLR, statesSortIndLR);


%%%% plotting the sorted transition matrix and lambda (three sorting based on:
%%%% Marginalized(over directions), RL, and LR sorting indices)

figure('Visible','off', 'units','normalized','position',[.1 .1 1.2 1.2]);

subplot(3,4,1)
imagesc(lambdaMa)
ylabel({'Marginalized over directions';[]; 'Units'},'fontsize', 20)
xlabel('States','fontsize', 20)
title({'Lambda Matrix'; '(sorted states)'} ,'fontsize', 20)
set(gca,'fontsize', 16);

subplot(3,4,2)
imagesc(1:numofStates, positionBins, PositionsLikelihoodMargin)
ylabel('Position(cm)','fontsize', 20)
xlabel('States','fontsize', 20)
set(gca,'YDir','normal')
title({'positions Likelihood'; '(sorted states)'} ,'fontsize', 20)
set(gca,'fontsize', 16);


subplot(3,4,3)
imagesc(transmatMa)
ylabel('States','fontsize', 20)
xlabel('States','fontsize', 20)
title({'Transition Matrix'; '(sorted states)'} ,'fontsize', 20)
set(gca,'fontsize', 16);



%%% RL
subplot(3,4,5)
imagesc(lambdaRL)
title({'Lambda Matrix'; '(sorted states)'} ,'fontsize', 20)
ylabel({'RL Place Fields';[]; 'Units'},'fontsize', 20)
set(gca,'fontsize', 16);

subplot(3,4,6)
y = PositionsLikelihoodRL(:, statesSortIndRL);
imagesc(1:numofStates, positionBins, y)
set(gca,'YDir','normal')
title({'positions Likelihood'; '(sorted states)'} ,'fontsize', 20)
ylabel('Position(cm)','fontsize', 20)
set(gca,'fontsize', 16);

subplot(3,4,7)
imagesc(transmatRL)
title({'Transition Matrix'; '(sorted states)'} ,'fontsize', 20)
ylabel('States','fontsize', 20)
set(gca,'fontsize', 16);

subplot(3,4,8)
imagesc(positionBins, 1:noActiveUnits, tuningRL)
title('Units tunings', 'fontsize', 16)
ylabel('Units','fontsize', 16)
set(gca,'fontsize', 16);



%%% LR
subplot(3,4,9)
imagesc(lambdaLR)
title({'Lambda Matrix'; '(sorted states)'} ,'fontsize', 20)
ylabel({'LR Place Fields';[]; 'Units'},'fontsize', 20)
set(gca,'fontsize', 16);

subplot(3,4,10)

y = PositionsLikelihoodLR(:, statesSortIndLR);
imagesc(1:numofStates, positionBins, y)
set(gca,'YDir','normal')
title({'positions Likelihood'; '(sorted states)'} ,'fontsize', 20)
ylabel('Position(cm)','fontsize', 20)
set(gca,'fontsize', 16);

subplot(3,4,11)
imagesc(transmatLR)
title({'Transition Matrix'; '(sorted states)'} ,'fontsize', 20)
ylabel('States','fontsize', 20)
set(gca,'fontsize', 16);

subplot(3,4,12)
imagesc(positionBins, 1:noActiveUnits, tuningLR)
title('Units tunings', 'fontsize', 16)
ylabel('Units','fontsize', 16)
set(gca,'fontsize', 16);


colormap(flipud(colormap('bone')))

saveas(gcf, [FileBase '/TrainedHMM.fig'])

% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% 
% print(h,[FileBase '/TrainedHMM'],'-dpdf','-r0')



%%%%% overlapped plots of the position likelihood distribution over
%%%%% leftward and rightward place fields; the sorting of the states was
%%%%% done based on the summation of likelihood over the two directions 


figure('Visible','off');

%%% plot colors
colormap1 = [ones(64,1) linspace(0,1,64)' linspace(0,1,64)']; %% graded shades of red
colormap2 = [linspace(0,1,64)' linspace(0,1,64)' ones(64,1)]; %% graded shades of blue


% colormax = prctile([PositionsLikelihoodLR(:); PositionsLikelihoodRL(:)], 95); 
colormax = min(max(PositionsLikelihoodLR(:)), max(PositionsLikelihoodRL(:))); %% color limits of imagesc (we want both the patterns be clear to the same degree)

ax1 = axes;
y = PositionsLikelihoodLR(:,statesSortIndMa);
p2 = imagesc(1:numofStates, positionBins, y, [0 colormax]);
p2.AlphaData = 1;

ax2 = axes;
x = PositionsLikelihoodRL(:,statesSortIndMa);
p1 = imagesc(1:numofStates, positionBins, x,[0 colormax]);
p1.AlphaData = 0.5;
title('blue: LR, red: RL','fontsize', 16)


linkaxes([ax1,ax2])
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];

colormap(ax1, flipud(colormap(colormap2)))
colormap(ax2, flipud(colormap(colormap1)))

saveas(gcf, [FileBase '/position_decoding_of_states(overlapped).fig'])


%%% Ratio of the likelihoods (RL/LR)

DirectionalLikratio = DirectionalLikratio(statesSortIndMa); %%% sorting the states based on the most likely position
mostLikelyDirection = zeros(size(DirectionalLikratio)); %%% the direction that a state is most probably representing

mostLikelyDirection(find(DirectionalLikratio >= 1)) = 1; %%% if the state is most likely represnting RL place fields
mostLikelyDirection(find(DirectionalLikratio < 1)) = -1; %%% same for LR place fields


figure('Visible','off');
subplot(6,1,1)
imagesc(1:numofStates, 1, mostLikelyDirection,[-1.7 1.7])
colormap jet
ylabel('Max Dir')
set(gca, 'XTick', [], 'YTick', [])

subplot(6,1, 2:6)
stairs(0:numofStates-1, log10(DirectionalLikratio), 'color', 'k')
hold on
line([1 numofStates],[1 1],'color', 'r')
hold on
line([1 numofStates],[-1 -1],'color', 'b')
grid on
ylabel('logLR(left / right)')
xlabel('States')

saveas(gcf, [FileBase '/states_most_Likely_corresponding_direction.fig'])



if nargout >= 1
    varargout{1} = mlPositionMa;
end

if nargout >= 2
    varargout{2} = statesSortIndMa;
end

if nargout >= 3
    varargout{3} = statesSortIndRL;
end

if nargout >= 4
    varargout{4} = statesSortIndLR;
end


%% added recently
if nargout >= 5
    varargout{5} = PositionsLikelihoodRL;
end

if nargout >= 6
    varargout{6} = PositionsLikelihoodLR;
end


end

