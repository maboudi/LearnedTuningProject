for ii = 1:10

figure; 
x0=0;
y0=0;
width=1200;
height=900;
set(gcf,'units','points','position',[x0,y0,width,height])

colormap(flipud(colormap('gray')))
% ploting the transmat and lambda matrices

clims_t = [0 1];
clims_l = [0 max([lambda(:);lambda_c(:);lambda_ts(:);lambda_p(:);lambda_i(:)])];

curr_transmat = transmat(:,:,ii); 
curr_transmat_c = transmat_c(:,:,ii);
curr_transmat_ts = transmat_ts(:,:,ii);
curr_transmat_p = transmat_p(:,:,ii);
curr_transmat_i = transmat_i(:,:,ii);

curr_lambda = lambda(:,:,ii); 
curr_lambda_c = lambda_c(:,:,ii);
curr_lambda_ts = lambda_ts(:,:,ii);
curr_lambda_p = lambda_p(:,:,ii);
curr_lambda_i = lambda_i(:,:,ii);


%%% sorting the states

[~, startIdx] = max(prior_r(:,:,ii)); 
sortInd_r = sortStates(curr_transmat, startIdx);
curr_transmat = curr_transmat(sortInd_r, sortInd_r);
curr_lambda = curr_lambda(:, sortInd_r);


[~, startIdx] = max(prior_r(:,:,ii)); 
sortInd_c = sortStates(curr_transmat_c, startIdx);
curr_transmat_c = curr_transmat_c(sortInd_r, sortInd_r);
curr_lambda_c = curr_lambda_c(:, sortInd_r);


[~, startIdx] = max(prior_r(:,:,ii)); 
sortInd_ts = sortStates(curr_transmat_ts, startIdx);
curr_transmat_ts = curr_transmat_ts(sortInd_r, sortInd_r);
curr_lambda_ts = curr_lambda_ts(:, sortInd_r);


[~, startIdx] = max(prior_r(:,:,ii)); 
sortInd_p = sortStates(curr_transmat_p, startIdx);
curr_transmat_p = curr_transmat_p(sortInd_r, sortInd_r);
curr_lambda_p = curr_lambda_p(:, sortInd_r);


[~, startIdx] = max(prior_r(:,:,ii)); 
sortInd_i = sortStates(curr_transmat_i, startIdx);
curr_transmat_i = curr_transmat_i(sortInd_r, sortInd_r);
curr_lambda_i = curr_lambda_i(:, sortInd_r);




numStates = size(curr_transmat, 1);
numUnits = size(curr_lambda, 1);

subplot(4,12,[1 2])
imagesc(curr_transmat, clims_t)
set(gca, 'fontsize', 10, 'TickDir','out')
set(gca,'XTick',0:20:numStates, 'YTick',0:20:numStates)

xlabel('State', 'fontsize', 10)
ylabel('State', 'fontsize', 10)
title('Real data', 'fontsize', 10)

subplot(4,12,[3 4])
imagesc(curr_transmat_c, clims_t)
set(gca, 'fontsize', 10, 'TickDir','out')
set(gca,'XTick',0:20:numStates, 'YTick',0:20:numStates, 'YTickLabel', [])


xlabel('State', 'fontsize', 10)
% ylabel('State', 'fontsize', 10)
title('Incoherent time shuffle', 'fontsize', 10)

subplot(4,12,[5 6])
imagesc(curr_transmat_ts, clims_t)
set(gca, 'fontsize', 10, 'TickDir','out')
set(gca,'XTick',0:20:numStates, 'YTick',0:20:numStates, 'YTickLabel', [])

xlabel('State', 'fontsize', 10)
% ylabel('State', 'fontsize', 10)
title('Coherent time shuffle', 'fontsize', 10)

subplot(4,12,[7 8])
imagesc(curr_transmat_p, clims_t)
set(gca, 'fontsize', 10, 'TickDir','out')
set(gca,'XTick',0:20:numStates, 'YTick',0:20:numStates, 'YTickLabel', [])

xlabel('State', 'fontsize', 10)
% ylabel('State', 'fontsize', 10)
title('Poisson simulated', 'fontsize', 10)

subplot(4,12,[9 10])
imagesc(curr_transmat_i, clims_t)
set(gca, 'fontsize', 10, 'TickDir','out')
set(gca,'XTick',0:20:numStates, 'YTick',0:20:numStates, 'YTickLabel', [])

xlabel('State', 'fontsize', 10)
% ylabel('State', 'fontsize', 10)
title('Unit identity shuffle', 'fontsize', 10)

subplot(4,12,[13 14])
imagesc(curr_lambda, clims_l)
set(gca, 'fontsize', 10, 'TickDir','out')
set(gca,'XTick',0:20:numStates, 'YTick',0:20:numUnits)

xlabel('State', 'fontsize', 10)
ylabel('Unit', 'fontsize', 10)

subplot(4,12,[15 16])
imagesc(curr_lambda_c, clims_l)
set(gca, 'fontsize', 10, 'TickDir','out')
set(gca,'XTick',0:20:numStates, 'YTick',0:20:numUnits, 'YTickLabel', [])

xlabel('State', 'fontsize', 10)
% ylabel('Unit', 'fontsize', 10)

subplot(4,12,[17 18])
imagesc(curr_lambda_ts, clims_l)
set(gca, 'fontsize', 10, 'TickDir','out')
set(gca,'XTick',0:20:numStates, 'YTick',0:20:numUnits, 'YTickLabel', [])

xlabel('State', 'fontsize', 10)
% ylabel('Unit', 'fontsize', 10)

subplot(4,12,[19 20])
imagesc(curr_lambda_p, clims_l)
set(gca, 'fontsize', 10, 'TickDir','out')
set(gca,'XTick',0:20:numStates, 'YTick',0:20:numUnits, 'YTickLabel', [])

xlabel('State', 'fontsize', 10)
% ylabel('Unit', 'fontsize', 10)

subplot(4,12,[21 22])
imagesc(curr_lambda_i, clims_l)
set(gca, 'fontsize', 10, 'TickDir','out')
set(gca,'XTick',0:20:numStates, 'YTick',0:20:numUnits, 'YTickLabel', [])

xlabel('State', 'fontsize', 10)
% ylabel('Unit', 'fontsize', 10)


subplot(4,12,[25 26 27])
transmatEntropy(transmat, transmat_c, transmat_ts, transmat_p, transmat_i, ii)

subplot(4,12,[37 38 39])
transmatEntropy2(transmat, transmat_c, transmat_ts, transmat_p, transmat_i, ii)


subplot(4,12,[28 29 30])
transmatGini(transmat, transmat_c, transmat_ts, transmat_p, transmat_i, ii)


subplot(4,12,[40 41 42])
transmatGini2(transmat, transmat_c, transmat_ts, transmat_p, transmat_i, ii)


subplot(4,12,[31 32 33])
lambdaGini_state(lambda, lambda_c, lambda_ts, lambda_p, lambda_i, ii)

subplot(4,12,[43 44 45])
lambdaGini_unit(lambda, lambda_c, lambda_ts, lambda_p, lambda_i, ii)


subplot(4,12,[34 35 36])
lambdaSparsity(lambda, lambda_c, lambda_ts, lambda_p, lambda_i, ii)

% subplot(4,12,[46 47 48])
% pathLength(transmat, transmat_c, transmat_ts, transmat_p, transmat_i, ii)

subplot(4,12,[46 47 48])
pathLen_thresh(transmat, transmat_c, transmat_ts, transmat_p, transmat_i, ii)

pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[pos(3), pos(4)])
print(gcf,['Urethane_' num2str(ii)],'-dpdf','-r0')


end


% 
% 
% for ii = 10
%     transmat2gephi(transmat(:,:,ii), ['real_' num2str(ii)])
%     transmat2gephi(curr_transmat_c, ['incoherent_' num2str(ii)])
%     transmat2gephi(curr_transmat_ts, ['coherent_' num2str(ii)])
%     transmat2gephi(curr_transmat_p, ['poisson_' num2str(ii)])
%     transmat2gephi(curr_transmat_i, ['unitIdentity_' num2str(ii)])
% end


