

seqMatchingIdx = struct('PRE', [], 'RUN', 'POST');
seqMatchingIdx = structfun(@(x)struct('data', [], 'poisson', []), seqMatchingIdx,'UniformOutput',false);

PBEsequence = struct('PRE', [], 'RUN', [], 'POST', []);

% PRE
[seqMatchingIdx.PRE.data, PBEsequence.PRE] = calSeqMatchingIndex(PREbinnedPBEs.data, runTemplate_LR, runTemplate_RL, []);
seqMatchingIdx.PRE.poisson = calSeqMatchingIndex(PREbinnedPBEs.p, runTemplate_LR, runTemplate_RL, []);

% RUN
[seqMatchingIdx.RUN.data, PBEsequence.RUN] = calSeqMatchingIndex(RUNbinnedPBEs.data, runTemplate_LR, runTemplate_RL, []);
seqMatchingIdx.RUN.poisson = calSeqMatchingIndex(RUNbinnedPBEs.p, runTemplate_LR, runTemplate_RL, []);

% POST
[seqMatchingIdx.POST.data, PBEsequence.POST] = calSeqMatchingIndex(POSTbinnedPBEs.data, runTemplate_LR, runTemplate_RL, []);
seqMatchingIdx.POST.poisson = calSeqMatchingIndex(POSTbinnedPBEs.p, runTemplate_LR, runTemplate_RL, []);

% 
% 
% % PRE
% [seqMatchingIdx.PRE.data, PBEsequence.PRE] = calSeqMatchingIndex(PREbinnedPBEs.data, runTemplate, [], []);
% seqMatchingIdx.PRE.poisson = calSeqMatchingIndex(PREbinnedPBEs.p, runTemplate, [], []);
% 
% % RUN
% [seqMatchingIdx.RUN.data, PBEsequence.RUN] = calSeqMatchingIndex(RUNbinnedPBEs.data, runTemplate, [], []);
% seqMatchingIdx.RUN.poisson = calSeqMatchingIndex(RUNbinnedPBEs.p, runTemplate, [], []);
% 
% % POST
% [seqMatchingIdx.POST.data, PBEsequence.POST] = calSeqMatchingIndex(POSTbinnedPBEs.data, runTemplate, [], []);
% seqMatchingIdx.POST.poisson = calSeqMatchingIndex(POSTbinnedPBEs.p, runTemplate, [], []);


bins = linspace(0, 100, 100);

periods = {'PRE';'RUN';'POST'};



figure;
set(gcf, 'Units', 'centimeters', 'position', [0 0 17.6 9])

for ii=1:3
    subplot(3,3,[3+ii 6+ii])

    data = (1-[seqMatchingIdx.(periods{ii}).data.LR.probability;seqMatchingIdx.(periods{ii}).data.RL.probability])*100;
    poisson = (1-[seqMatchingIdx.(periods{ii}).poisson.LR.probability;seqMatchingIdx.(periods{ii}).poisson.RL.probability])*100;
    
%     data = (1-seqMatchingIdx.(periods{ii}).data.probability)*100;
%     poisson = (1-seqMatchingIdx.(periods{ii}).poisson.probability)*100;
   
    significanceMarker_poisson = calSigMarker(data, poisson);
    
    countsDATA  = hist(data , bins); countsDATA = countsDATA./sum(countsDATA); cumCountsDATA = cumsum(countsDATA);
    countsPOISS = hist(poisson, bins); countsPOISS = countsPOISS./sum(countsPOISS); cumCountsPOISS = cumsum(countsPOISS);

    plot(bins, cumCountsPOISS, 'color', [.8 .8 .8], 'linewidth', 3)
    hold on
    plot(bins, cumCountsDATA, 'color', [220 85 85]/255, 'linewidth', 3)
    

    xlabel('match probability (%)', 'fontsize', 10)
    title({periods{ii}; sprintf('(n=%d)', numel(data))}, 'fontsize', 10)
    
    if ii==1; ylabel('cumulative ratio of PBEs', 'fontsize', 10); end
    
    
    ylim([-0.05 1])
    xlim([-5 100])

    axis square
    set(gca, 'fontsize', 10, 'linewidth', 2, 'TickDir', 'out', 'TickLength',[0.01, 0.01])
    set(gca, 'XTick', bins([1 25 50 75 100]), 'YTick', linspace(0, 1, 5), ...
    'XTickLabels', {'0', '', '', '', '100'}, 'YTickLabels', {'0', '', '', '', '1'}, 'box', 'off')

    grid on
    text(2, 0.9, {sprintf('vs poisson  %s', significanceMarker_poisson)}, 'fontsize', 6) 
    

    dim = [0.005 0.85 0.5 .1];

    annotation('textbox',dim,'String',{fileinfo.name; 'Sequence matching index as a percentile of all or 10000 unit permutations'; ''}, 'FitBoxToText','on', 'Interpreter', 'none');

end


function significanceMarker = calSigMarker(seqScore1, seqScore2)

pval = ranksum(seqScore1, seqScore2, 'tail', 'right');

if pval < 0.001
    significanceMarker = '< 0.001';
elseif pval < 0.01
    significanceMarker = '< 0.01';
elseif pval < 0.05
    significanceMarker = '< 0.05';
else
    significanceMarker = sprintf('%.2f', pval);
end

end

 