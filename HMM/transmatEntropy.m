function transmatEntropy(transmat, transmat_c, transmat_ts, transmat_p, transmat_i, ii)

numStates = size(transmat, 1);

real = transmat(:,:,ii);
unitCycle = transmat_c(:,:,ii);
timeSwap = transmat_ts(:,:,ii);
poiss = transmat_p(:,:,ii);
identity = transmat_i(:,:,ii);

% remove the self-transitions
% real = real - real .* eye(size(real));
% real = real./repmat(sum(real, 2), [1 numStates]);
% real = real+1e-4;
% 
% unitCycle = unitCycle - unitCycle .* eye(size(unitCycle));
% unitCycle = unitCycle./repmat(sum(unitCycle, 2), [1 numStates]);
% unitCycle = unitCycle + 1e-4;
% 
% timeSwap = timeSwap - timeSwap .* eye(size(timeSwap));
% timeSwap = timeSwap./repmat(sum(timeSwap, 2), [1 numStates]);
% timeSwap = timeSwap + 1e-4;
% 
% 
% poiss = poiss - poiss .* eye(size(poiss));
% poiss = poiss./repmat(sum(poiss, 2), [1 numStates]);
% poiss = poiss + 1e-4;
% 
% identity = identity - identity .* eye(size(identity));
% identity = identity./repmat(sum(identity, 2), [1 numStates]);
% identity = identity + 1e-4;

% calculating entropy for each state

realEntropy = sum(-real .* log2(real), 2); % probability of transitions from each state
unitCycleEntropy = sum(-unitCycle .* log2(unitCycle), 2);
timeSwapEntropy = sum(-timeSwap .* log2(timeSwap), 2);
poissEntropy = sum(-poiss .* log2(poiss), 2);
identityEntropy = sum(-identity .* log2(identity), 2);



% histograms

temp = [realEntropy unitCycleEntropy timeSwapEntropy poissEntropy identityEntropy];
bins = linspace(min(temp(:)), max(temp(:)), 50);

realH = hist(realEntropy, bins);
UnitCycleH = hist(unitCycleEntropy, bins);
timeSwapH = hist(timeSwapEntropy, bins);
poissH = hist(poissEntropy, bins);
identityH = hist(identityEntropy, bins);


% statistical tests for comparing the distributions

realVScycleP = ranksum(realEntropy, unitCycleEntropy);
realVStimeSwapP = ranksum(realEntropy, timeSwapEntropy);
realVSpoissP = ranksum(realEntropy, poissEntropy);
realVSidentityP = ranksum(realEntropy, identityEntropy);


% smoothing and plot the distributions
sigma = 2;
halfwidth = 2 * sigma;
smoothwin = gausswindow(sigma, halfwidth);

realH = conv(realH, smoothwin, 'same');
UnitCycleH = conv(UnitCycleH, smoothwin, 'same');
timeSwapH = conv(timeSwapH, smoothwin, 'same');
poissH = conv(poissH, smoothwin, 'same');
identityH = conv(identityH, smoothwin, 'same');

% realH = realH ./ sum(realH);
% UnitCycleH = UnitCycleH ./ sum(UnitCycleH);
% timeSwapH = timeSwapH ./ sum(timeSwapH);
% poissH = poissH ./ sum(poissH);
% identityH = identityH ./ sum(identityH);

% figure
hold on

plot(bins, realH, 'k', 'linewidth', 1.5)
plot(bins, UnitCycleH, 'c','linewidth', 1)
plot(bins, timeSwapH, 'b','linewidth', 1)
plot(bins, poissH, 'g','linewidth', 1)
plot(bins, identityH, 'm','linewidth', 1)


% significance bars
[realPy, realPx]= max(realH);
[UnitCyclePy, UnitCyclePx]= max(UnitCycleH);
[timeSwapPy, timeSwapPx]= max(timeSwapH);
[poissPy, poissPx]= max(poissH);
[identityPy, identityPx]= max(identityH);

firstmax = max([realPy, UnitCyclePy, timeSwapPy, poissPy, identityPy])*1.05;
textsize = 8;
% Wilcoxon ranksum test
if realVScycleP < 0.05
    y = firstmax;
    plot([bins(realPx) bins(UnitCyclePx)],[y y],'k','linewidth', 1.5)
    
    if realVScycleP < 0.01
        psign = '**';
    else
        psign = '*';
    end
        
    text(mean([bins(realPx), bins(UnitCyclePx)]), y+0.05, psign, 'fontsize', textsize, 'HorizontalAlignment', 'center')

end

if realVStimeSwapP < 0.05
    y = firstmax+0.5;
    plot([bins(realPx) bins(timeSwapPx)],[y y],'k','linewidth', 1.5)
    
    if realVStimeSwapP < 0.01
        psign = '**';
    else
        psign = '*';
    end
        
    text(mean([bins(realPx), bins(timeSwapPx)]), y+0.05, psign, 'fontsize', textsize, 'HorizontalAlignment', 'center')

end

if realVSpoissP < 0.05
    y = firstmax+1;
    plot([bins(realPx) bins(poissPx)],[y y],'k','linewidth', 1.5)
    
    if realVSpoissP < 0.01
        psign = '**';
    else
        psign = '*';
    end
        
    text(mean([bins(realPx), bins(poissPx)]), y+0.05, psign, 'fontsize', textsize, 'HorizontalAlignment', 'center')

    
end

if realVSidentityP < 0.05
    y = firstmax+1.5;
    plot([bins(realPx) bins(identityPx)],[y y],'k','linewidth', 1.5)
    
    if realVSidentityP < 0.01
        psign = '**';
    else
        psign = '*';
    end
        
    text(mean([bins(realPx), bins(identityPx)]), y+0.05, psign, 'fontsize', textsize, 'HorizontalAlignment', 'center')

end

hold off

set(gca, 'fontsize', 10)
xlabel('Entropy', 'fontsize', 10)
ylabel('Counts', 'fontsize', 10)
title('States transition sparsity', 'fontsize', 10)