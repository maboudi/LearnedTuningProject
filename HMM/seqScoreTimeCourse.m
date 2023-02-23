function seqScoreTimeCourse(seqscore, PBEs, wholePeriod, currTitle, fileinfo)

hold on 

nChuncks = 5;
edges    = linspace(wholePeriod(1), wholePeriod(2), nChuncks+1);


[~, idx] = histc(PBEs(:, 3), edges); 

med = zeros(5,1);
% lq  = zeros(5,1);
% hq  = zeros(5,1);

for ii = 1:nChuncks
   
    chunkScores = seqscore(idx == ii, :);
    chunkScores = chunkScores(:);
    
    med(ii) = median(chunkScores);
%     lq(ii)  = prctile(chunkScores, 25);
%     hq(ii)  = prctile(chunkScores, 75);    
    
end

centers  = edges(1:nChuncks) + (edges(2) - edges(1))/2 - wholePeriod(1);
centers  = centers/fileinfo.Fs/60/60; % in hours

plot(centers, med, '-square', 'color', [0 100 0]/255,  'linewidth', 2);
% plot(centers, lq , 'color', 'k', 'linestyle', ':', 'linewidth', 1);
% plot(centers, hq , 'color', 'k', 'linestyle', ':', 'linewidth', 1);


xlabel('time (hr)')
title(currTitle)

ylim([min(med)-10 min(max(med)+10, 100)])

set(gca, 'fontsize', 10, 'linewidth', 2, 'TickDir', 'out', 'TickLength',[0.01, 0.01])

axis square

% 
% for ii = 1:nChuncks
%     
%    scores1 = seqscore(idx == ii, :);
%    for jj = 1 :nChuncks
%        
%        scores2 = seqscore(idx == jj, :);
%        [p(ii,jj), h(ii, jj)] = ranksum(scores1(:), scores2(:));
%        
%    end 
% end

end