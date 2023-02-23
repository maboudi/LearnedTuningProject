function plotExampleHMM(transmatPRE, lambdaPRE, transmatRUN, lambdaRUN, transmatPOST, lambdaPOST, transmatActiveRUN, lambdaActiveRUN, mainDir)


% set the color limits

tmat_temp = [transmatPRE(:); transmatRUN(:); transmatPOST(:); transmatActiveRUN(:)];
tmat_clim = [min(tmat_temp) max(tmat_temp)];


omat_temp = [lambdaPRE(:); lambdaRUN(:); lambdaPOST(:); lambdaActiveRUN(:)];
omat_clim = [min(omat_temp) max(omat_temp)];


figure;
set(gcf, 'units', 'centimeters', 'position', [0 0 9 15])


% PRE
subplot(4,2,1)
imagesc(transmatPRE) % , tmat_clim
% set(gca, 'YDir', 'normal')
xlabel('State j', 'fontsize', 8)
ylabel({'PRE';'';'State i'}, 'fontsize', 8)
title('Transition matrix', 'fontsize', 8)

subplot(4,2,2)
imagesc(lambdaPRE) % , omat_clim
set(gca, 'YDir', 'normal')
xlabel('State', 'fontsize', 8)
ylabel('Unit', 'fontsize', 8)
title('Observation matrix', 'fontsize', 8)


% RUN
subplot(4,2,3)
imagesc(transmatRUN)
% set(gca, 'YDir', 'normal')
xlabel('State j', 'fontsize', 8)
ylabel({'RUN';'';'State i'}, 'fontsize', 8)


subplot(4,2,4)
imagesc(lambdaRUN)
set(gca, 'YDir', 'normal')
xlabel('State', 'fontsize', 8)
ylabel('Unit', 'fontsize', 8)


% POST
subplot(4,2,5)
imagesc(transmatPOST)
% set(gca, 'YDir', 'normal')
xlabel('State j', 'fontsize', 8)
ylabel({'POST';'';'State i'}, 'fontsize', 8)


subplot(4,2,6)
imagesc(lambdaPOST)
set(gca, 'YDir', 'normal')
xlabel('State', 'fontsize', 8)
ylabel('Unit', 'fontsize', 8)



% Active RUN
subplot(4,2,7)
imagesc(transmatActiveRUN)
% set(gca, 'YDir', 'normal')
xlabel('State j', 'fontsize', 8)
ylabel({'Active RUN';'';'State i'}, 'fontsize', 8)


subplot(4,2,8)
imagesc(lambdaActiveRUN)
set(gca, 'YDir', 'normal')
xlabel('State', 'fontsize', 8)
ylabel('Unit', 'fontsize', 8)



colormap('jet')


mkdir(fullfile(mainDir, 'exampleHMMs'))

filename = fullfile(mainDir, 'exampleHMMs', 'exampleHMMS');

savepdf(gcf, filename, '-dsvg')
saveas(gcf, filename, 'epsc')

end