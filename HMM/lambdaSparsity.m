function lambdaSparsity(lambda, lambda_c, lambda_ts, lambda_p, lambda_i, ii)


real = lambda(:,:,ii);
unitCycle = lambda_c(:,:,ii);
timeSwap = lambda_ts(:,:,ii);
poiss = lambda_p(:,:,ii);
identity = lambda_i(:,:,ii);

% Calculate the probability of each unit firing at least once given each
% state (:: 1-exp(-lambda))

real = 1 - exp(-real);
unitCycle = 1 - exp(-unitCycle);
timeSwap = 1 - exp(-timeSwap);
poiss = 1 - exp(-poiss);
identity = 1 - exp(-identity);

% find the joint probability of units 
% Pr(u1 > 1, u2 > 1, ... , un > 1) = Pr(u1 > 1). Pr(u2 > 1). ... Pr(un > 1)

realjoint = log10(prod(real));
unitCyclejoint = log10(prod(unitCycle));
timeSwapjoint = log10(prod(timeSwap));
poissjoint = log10(prod(poiss));
identityjoint = log10(prod(identity));


% histograms

temp = [realjoint unitCyclejoint timeSwapjoint poissjoint identityjoint];
bins = linspace(min(temp(:)), max(temp(:)), 50);

realH = hist(realjoint, bins);
UnitCycleH = hist(unitCyclejoint, bins);
timeSwapH = hist(timeSwapjoint, bins);
poissH = hist(poissjoint, bins);
identityH = hist(identityjoint, bins);


% statistical tests for comparing the distributions

realVScycleP = ranksum(realjoint, unitCyclejoint);
realVStimeSwapP = ranksum(realjoint, timeSwapjoint);
realVSpoissP = ranksum(realjoint, poissjoint);
realVSidentityP = ranksum(realjoint, identityjoint);


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

firstmax = max([realPy, UnitCyclePy, timeSwapPy, poissPy, identityPy])*1.02;
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
xlabel('Probability (log)', 'fontsize', 10)
ylabel('Counts', 'fontsize', 10)
title('Lambda sparsity', 'fontsize', 10)





