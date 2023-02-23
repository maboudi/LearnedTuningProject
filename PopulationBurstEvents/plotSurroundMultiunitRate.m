function plotSurroundMultiunitRate(sdat, PBEs , fileinfo, subfolder)


noPBEs2plot = min([100 length(PBEs)]);
randIdx = randperm(length(PBEs));
randIdx = randIdx(1:noPBEs2plot);

PBEs2 = floor((PBEs - fileinfo.tbegin)*1000/fileinfo.Fs);

nRos = ceil(sqrt(noPBEs2plot));
nCols = nRos;
cornerPanel = (nRos-1)*nCols+1;

addedflankingPeriod = 200; % in ms 


h = figure; 
set(gcf, 'Units', 'centimeters', 'position', [0 0 18.4 11])

for ii = 1: noPBEs2plot
   subplot(nRos, nCols, ii)
   hold on
   
   str = PBEs2(randIdx(ii), 1);
   strPlus = str - addedflankingPeriod;
   
   endd = PBEs2(randIdx(ii), 2);
   enddPlus = endd + addedflankingPeriod;
   
   try
   plot(0:(enddPlus-strPlus), sdat(strPlus:enddPlus), 'k', 'linewidth', 1)
   catch
    plot(0:(enddPlus-strPlus), sdat(strPlus+1:enddPlus), 'k', 'linewidth', 1)
      
   end
   
   
   yAxisRange = ylim;
   line([addedflankingPeriod addedflankingPeriod], [yAxisRange(1) yAxisRange(2)], 'linestyle', '-.', 'color', 'r', 'linewidth', 1)
   
   line([endd-str+addedflankingPeriod endd-str+addedflankingPeriod], [yAxisRange(1) yAxisRange(2)], 'linestyle', '-.', 'color', 'r', 'linewidth', 1)
   
   hold off
   
   xticks([addedflankingPeriod endd-str+addedflankingPeriod])
   xticklabels({'0', sprintf('   %.0f', endd - str)})
   
   yticks([yAxisRange(1) yAxisRange(2)])
   yticklabels({sprintf('%.0f', yAxisRange(1)), sprintf('%.0f', yAxisRange(2))})
   
   set(gca, 'fontsize', 5)
   
   if ii == cornerPanel    
       xlabel('time (ms)', 'fontsize', 5)
       ylabel('pop. rate (a.u.)', 'fontsize', 5)
   end
   
end


filename = fullfile(subfolder, 'PBEs_MultiunitRateProfiles');
savepdf(gcf, filename, '-dpdf')



end