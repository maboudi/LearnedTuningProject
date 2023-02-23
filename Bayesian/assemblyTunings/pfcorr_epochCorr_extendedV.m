

% to examine if the learning tuning of individual units is correlated
% between different epochs


sessionNumbers = [1:5];


parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/concat_GreatLakes_datasets_temp/';


%% preconfigureing the plots


% all sessions

plotheight = 25; %% in cm
plotwidth = 20;
fontsize = 4;

nor = 3;
noc = 3; %% number of events per column

leftmargin = 2;  rightmargin = 2;    topmargin = 2;    bottommargin = 2;

gapc = 2;
gapr = 3;

sub_pos = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, gapr, gapc);


f=figure;
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
    





% for each session

plotwidth_sess = (plotwidth - (rightmargin+leftmargin) - (noc-1) * gapc)/noc;
plotheight_sess = (plotheight - (topmargin+bottommargin) - (nor-1) * gapr)/nor;

nor_sess = 6;
noc_sess = 4;

leftmargin_sess = 0;  rightmargin_sess = 0;     bottommargin_sess = 0;    topmargin_sess = 0;
gapr_sess = 0;    gapc_sess = 0;

sub_pos_sess = subplot_pos(plotwidth_sess, plotheight_sess, leftmargin_sess, rightmargin_sess, bottommargin_sess, topmargin_sess, noc_sess, nor_sess, gapr_sess, gapc_sess); %% this gives us the positions for the subplots, each corresponded to an event







% for Bayesian Learned tunings

plotwidth_blt  = sub_pos_sess{1,4}(1) + sub_pos_sess{1,4}(3) - sub_pos_sess{1,1}(1);
plotheight_blt = sub_pos_sess{1,1}(2) + sub_pos_sess{1,1}(4) - sub_pos_sess{3,1}(2);

nor_blt = 1;
noc_blt = 4;

leftmargin_blt = 0;  rightmargin_blt = 0;     bottommargin_blt = 0;    topmargin_blt = 0;
gapr_blt = 0.015;    gapc_blt = 0.015;

sub_pos_blt = subplot_pos(plotwidth_blt, plotheight_blt, leftmargin_blt, rightmargin_blt, bottommargin_blt, topmargin_blt, noc_blt, nor_blt, gapr_blt, gapc_blt); %% this gives us the positions for the subplots, each corresponded to an event



for ii = 1 : nor_blt
    for jj = 1 : noc_blt

        sub_pos_blt{ii,jj}(1) = (sub_pos_blt{ii,jj}(1) * plotwidth_blt + sub_pos_sess{1}(1) * plotwidth_sess + sub_pos{1}(1) * plotwidth) / plotwidth;
        sub_pos_blt{ii,jj}(2) = (sub_pos_blt{ii,jj}(2) * plotheight_blt + sub_pos_sess{1}(2) * plotheight_sess + sub_pos{1}(2) * plotheight) / plotheight;
        sub_pos_blt{ii,jj}(3) = sub_pos_blt{ii,jj}(3) * plotwidth_blt / plotwidth;
        sub_pos_blt{ii,jj}(4) = sub_pos_blt{ii,jj}(4) * plotheight_blt / plotheight;
    end
end



axes('position',sub_pos_blt{1,1},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','off','Layer','top');
axes('position',sub_pos_blt{1,2},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','off','Layer','top');
axes('position',sub_pos_blt{1,3},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','off','Layer','top');
axes('position',sub_pos_blt{1,4},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','off','Layer','top');



% for assemblyTuning histogram

plotwidth_at  = plotwidth_sess;
plotheight_at = plotheight_blt;

sub_pos_at(1) = (sub_pos_sess{4,1}(1) * plotwidth_sess + sub_pos{1}(1) * plotwidth)/plotwidth;
sub_pos_at(2) = (sub_pos_sess{4,1}(2) * plotheight_sess + sub_pos{1}(2) * plotheight)/plotheight;
sub_pos_at(3) = (sub_pos_sess{4,3}(1) + sub_pos_sess{4,3}(3) - sub_pos_sess{4,1}(1))*plotwidth_sess/plotwidth;
sub_pos_at(4) = sub_pos_sess{4,4}(2)*plotheight_sess/plotheight;

axes('position',sub_pos_at,'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','off','Layer','top');






replayExist = 1;

%% assembly tuning correlation across different epochs

for isession = 1:numel(sessionNumbers)
    

sessionNumber = sessionNumbers(isession);

rr = dir(parentDir);
rr = rr([rr.isdir]);

sessionName = rr(sessionNumber+2).name;

basePath = fullfile(parentDir, sessionName);


load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTunings_allPBEs.mat']))

try
    load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTuning_vs_replayScores.mat']))
catch
    replayExist = 0;
end

epochNames = {'pre'; 'run'; 'post'};
nEpochs    = numel(epochNames);



assemblyTuning_epochCorr       = nan(nEpochs, nEpochs);
pval                           = nan(nEpochs, nEpochs);
assemblyTuning_epochCorr_prctl = nan(nEpochs, nEpochs);


assemblyTuning_epochCorr_rep       = nan(nEpochs, nEpochs);
pval_rep                           = nan(nEpochs, nEpochs);
assemblyTuning_epochCorr_rep_prctl = nan(nEpochs, nEpochs);


sub_pos2 = subplot_pos(plotwidth2, plotheight2, leftmargin2, rightmargin2, bottommargin2, topmargin2, noc2, nor2, gapr2, gapc2); 
for ii = 1 : nor2
    for jj = 1 : noc2

        sub_pos2{ii,jj}(1) = (sub_pos2{ii,jj}(1) * plotwidth2 + sub_pos{isession}(1) * plotwidth) / plotwidth;
        sub_pos2{ii,jj}(2) = (sub_pos2{ii,jj}(2) * plotheight2 + sub_pos{isession}(2) * plotheight) / plotheight;
        sub_pos2{ii,jj}(3) = sub_pos2{ii,jj}(3) * plotwidth2 / plotwidth;
        sub_pos2{ii,jj}(4) = sub_pos2{ii,jj}(4) * plotheight2 / plotheight;
    end
end


tt = 0;

for ie1 = 1:nEpochs-1
    for ie2 = ie1+1:nEpochs
        
        tt = tt+1;
        %%% correlatin of assembly tunings across different epochs 
        %%% calculated based on all PBEs
        
        
        % real data
        [assemblyTuning_epochCorr(ie1, ie2), pval(ie1, ie2)] = calcorr(assemblyTuningPFcorr.(epochNames{ie1}).data, assemblyTuningPFcorr.(epochNames{ie2}).data); 
        
        
        % plot the scatter for the real data
        axes('position',sub_pos2{tt,1},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','off','Layer','top');
        
        
        
        var1 = assemblyTuningPFcorr.(epochNames{ie1}).data;
        var2 = assemblyTuningPFcorr.(epochNames{ie2}).data;
        
        scatter_regline(var1, var2)
        xlabel((epochNames{ie1}), 'fontsize', fontsize)
        ylabel((epochNames{ie2}), 'fontsize', fontsize)
        text(0.7, 0.7, sprintf('rho=%.2f\n(p=%.3f)', assemblyTuning_epochCorr(ie1, ie2), pval(ie1, ie2)), 'fontsize', 6, 'color','r')
        
        if ie1 == 1 && ie2 == 2
            
            sessionName(sessionName == '_') = '-';
            title({sessionName; '';'all PBEs'}, 'fontsize', 8, 'fontweight', 'normal')
            
        end
        
        % --------------------
        
        
        % ui surrogate
        assemblyTuning_epochCorr_ui = zeros(100, 1);
        for is = 1:100
            assemblyTuning_epochCorr_ui(is) = calcorr(assemblyTuningPFcorr.(epochNames{ie1}).ui(:,is), assemblyTuningPFcorr.(epochNames{ie2}).ui(:, is));
        end
        
        assemblyTuning_epochCorr_prctl(ie1, ie2) = numel(find(assemblyTuning_epochCorr_ui < assemblyTuning_epochCorr(ie1, ie2)))/100;
        
        
        % plot the histogram for ui surrogates and real data as a percntile
        % of the ui data
        
        insetPos = sub_pos2{tt, 1};
        
        
        insetPos(1) = insetPos(1)+0.005;
        insetPos(2) = insetPos(2)+0.005;
        insetPos(3) = insetPos(3)*0.5;
        insetPos(4) = insetPos(4)*0.25;
        
        axes('position', insetPos,'XGrid','off','XMinorGrid','off','FontSize', 6,'Box','off','Layer','top');
        
        plotbars(assemblyTuning_epochCorr_ui, assemblyTuning_epochCorr(ie1, ie2), assemblyTuning_epochCorr_prctl(ie1, ie2))
        
        % ------------------
        
        if replayExist
        
            %%% correlatin of assembly tunings across different epochs 
            %%% calculated based on only PBEs with highest replay scores (>75%)


            % real data
            [assemblyTuning_epochCorr_rep(ie1, ie2), pval_rep(ie1, ie2)] = calcorr(asTuningPFcorr_sub.(epochNames{ie1}).data.rt_ui(:,4), asTuningPFcorr_sub.(epochNames{ie2}).data.rt_ui(:,4)); 


            %  plot the scatter for real data

            axes('position',sub_pos2{tt,2},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','off','Layer','top');



            var1 = asTuningPFcorr_sub.(epochNames{ie1}).data.rt_ui(:,4);
            var2 = asTuningPFcorr_sub.(epochNames{ie2}).data.rt_ui(:,4);

            scatter_regline(var1, var2)
            xlabel((epochNames{ie1}), 'fontsize', fontsize)
    %         ylabel((epochNames{ie2}), 'fontsize', fontsize)
            text(0.7, 0.7, sprintf('rho=%.2f\n(p=%.3f)', assemblyTuning_epochCorr_rep(ie1, ie2), pval_rep(ie1, ie2)), 'fontsize', 6, 'color', 'r')

            if ie1 == 1 && ie2 == 2
                title('replay PBEs', 'fontsize', 8, 'fontweight', 'normal')
            end

            % -------------------------------


            % ui surrogates
            assemblyTuning_epochCorr_rep_ui = zeros(100, 1);

            for is = 1:100
                assemblyTuning_epochCorr_rep_ui(is) = calcorr(asTuningPFcorr_sub.(epochNames{ie1}).ui.rt_ui(:, 4, is), asTuningPFcorr_sub.(epochNames{ie2}).ui.rt_ui(:, 4, is));
            end

            assemblyTuning_epochCorr_rep_prctl(ie1, ie2) = numel(find(assemblyTuning_epochCorr_rep_ui < assemblyTuning_epochCorr_rep(ie1, ie2)))/100;


            % plot the histogram for ui surrogates and real data as a
            % percentile of the uni data

            insetPos = sub_pos2{tt, 2};


            insetPos(1) = insetPos(1)+0.005;
            insetPos(2) = insetPos(2)+0.005;
            insetPos(3) = insetPos(3)*0.5;
            insetPos(4) = insetPos(4)*0.25;

            axes('position', insetPos,'XGrid','off','XMinorGrid','off','FontSize', 6,'Box','off','Layer','top');

            plotbars(assemblyTuning_epochCorr_rep_ui, assemblyTuning_epochCorr_rep(ie1, ie2), assemblyTuning_epochCorr_rep_prctl(ie1, ie2))


            % --------------------------------
        end

    end
    
end


end




%% functions

function [corrCoeff, pval] = calcorr(var1, var2)

idx = ~isnan(var1) & ~isnan(var2);

[corrCoeff, pval] = corr(var1(idx), var2(idx), 'type', 'spearman');

end

function scatter_regline(var1, var2)


idx = ~isnan(var1) & ~isnan(var2);

var1 = var1(idx);
var2 = var2(idx);

hold on
scatter(var1, var2, 5, 'filled', 'k', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none')

% b1 = var1\var2;
% var2Calc = var1*b1;

VAR1 = [ones(size(var1)) var1];
b = VAR1\var2;

var2Calc = VAR1*b;
plot(var1, var2Calc, '--r')

xlim([-1, 1])
ylim([-1, 1])

end

function plotbars(shuffles, real, percentile)

ax1 = gca;


bins = linspace(min(shuffles), max(shuffles), 10);

counts = hist(shuffles, bins);

counts = counts/max(counts);


bar(bins, counts, 'EdgeColor', 'k', 'FaceColor', 'none', 'FaceAlpha', 0.5)

hold on

line([real real], [0 1], 'color', 'r', 'linewidth', 1)

ax1.FontSize = 4;
ax1.YAxis.Visible = 'off';
ax1.Box = 'off';


text(real, 0.7, sprintf('%%%.1f', percentile*100), 'fontsize', 6, 'color', 'r')
text(min(shuffles), 1, 'ui', 'color', 'k', 'fontsize', 6)



end

