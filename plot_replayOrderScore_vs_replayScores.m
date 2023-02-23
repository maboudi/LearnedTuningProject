close all;

epochNames = {'pre'; 'run'; 'post'};
methods = {'rt_ts'; 'rt_ui'; 'rt_pf'; 'rt_ds'; 'wc_ts'; 'wc_ui'; 'wc_pf'; 'wc_ds'};
quartiles = {'q1'; 'q2'; 'q3'; 'q4'};


figure;
for ie = 1:numel(epochNames)
   
    currEpoch = epochNames{ie};
    PBEs.(currEpoch) = PBEInfo_replayScores(strcmp({PBEInfo_replayScores.epoch}, currEpoch));
    
    
    replayOrderScores_all = abs([PBEs.(currEpoch).replayOrderScore_prctile]');
    
    
    for im = 1: numel(methods)
        
        ax = subplot(3,8, (ie-1)*numel(methods)+im);
        currMethod = methods{im};
        
        replayScores_all = [PBEs.(currEpoch).(currMethod)]';
        for iq = 1:numel(quartiles)
            
            idx = find(replayScores_all >= (iq-1)*25.01 & replayScores_all < iq*25.01);
            replayOrderScore.(currEpoch).(methods{im}).(quartiles{iq}) = replayOrderScores_all(idx);
        end
        
        h = violinplot(replayOrderScore.(currEpoch).(methods{im}));
        for ii = 1:4
            h(ii).ViolinColor = 'k';
            h(ii).EdgeColor   = 'none';
            h(ii).ShowData    = 0;
            
            fn = fieldnames(h(ii));
            
            for ff = 1:numel(fn)
                if isprop(h(ii).(fn{ff}), 'XData')
                    try
                        h(ii).(fn{ff}).XData = (h(ii).(fn{ff}).XData - ii)*30 + ii*25 - 12.5;
                    catch
                        continue
                    end
                end
            end
            
            h(ii).BoxWidth = h(ii).BoxWidth * 5;
            h(ii).BoxPlot.FaceColor = [0.4 0.4 0.4];
            h(ii).BoxPlot.EdgeColor = 'none';
            
        end
        
        xticks([12.5 37.5 67.5 87.5])
        xticklabels({'0-25%'; '25-50%'; '50-75%'; '75-100%'})
        xtickangle(45)
        
        if im > 1
            ax.YAxis.Visible = 'off';
        end
        
        % correlation
        [rho, pval] = corr(replayScores_all, replayOrderScores_all, 'type', 'spearman');
        
        yl = ylim;
        
        text(10, 0.95*yl(2), sprintf('rho = %.2f', rho), 'color', 'b')
        if pval < 0.001
            text(10, 0.90*yl(2), 'pval < 0.001', 'color', 'b')
        else
            text(10, 0.90*yl(2), sprintf('pval = %.3f', pval), 'color', 'b')
        end
        
        % regression line
        X = [ones(length(replayScores_all), 1) replayScores_all];
        b = X\replayOrderScores_all;
        
        yCalc = X*b;
        plot(replayScores_all, yCalc, 'r')
        
        
        if ie == 1
            tl  = currMethod;
            tl(tl == '_') = '-';
            title(tl)
        end
        
        if im == 1
            ylabel({currEpoch; 'directional replay index (percentile)'})
        end
        
        grid on
    end
    
end

% idx = [PRE_PBEs.rt_ts] < 25;
% median(abs([PRE_PBEs(idx).replayOrderScore]))







