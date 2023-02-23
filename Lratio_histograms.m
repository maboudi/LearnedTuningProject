% distribution of L-ratio by pooling across sessions


clear; clc; close all

parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/assemblyTuning_finalResults/';

sub = dir(parentDir);


currSessions = [1:5 7:9 10 13:15 6];
nSessions = numel(currSessions);

Lratio_allSessions = cell(nSessions);

for iSession = 1:nSessions

    currSession = currSessions(iSession);

    sessionName = sub(currSession+2).name;

    load(fullfile(parentDir, sessionName, 'spikes', [sessionName '.clusterQuality.mat']))

    Lratios = [];

    for sh = 1:numel(clusterQuality)
        Lratios = [Lratios; clusterQuality(sh).Lratio(:)];
    end
    Lratio_allSessions{iSession} = Lratios(~isnan(Lratios));

end


%% for RatN calculate the overall distribution of L-ratio

plotwidth  = 240;
plotheight = 100;
fontsize   = 6;

close all
f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

hold on


Lratio_exampleSession = Lratio_allSessions{6};
Lratio_exampleSession_median = nanmedian(Lratio_exampleSession); 

Lratio_pooled = cell2mat(Lratio_allSessions);
Lratio_pooled_median = nanmedian(Lratio_pooled);




[cdf_example, x_example] = ecdf(Lratio_allSessions{6});
[x_example, idx] = unique(x_example);
cdf_example = cdf_example(idx);


[cdf_pooled, x_pooled] = ecdf(Lratio_pooled);
[x_pooled, idx] = unique(x_pooled);
cdf_pooled = cdf_pooled(idx);



x = logspace(-20, 10);
cdf_example = interp1(x_example, cdf_example, x);
cdf_pooled  = interp1(x_pooled, cdf_pooled, x);


h1 = semilogx(x, cdf_example, 'color', [0.3 0.3 1], 'linewidth', 1);
h2 = semilogx(x, cdf_pooled, 'color', [1 0.3 0.3], 'linewidth', 1);

axis tight


ylim([0 1])
yl = ylim;
xl = xlim;

line([xl(1) Lratio_exampleSession_median], [0.5 0.5], 'linewidth', 0.5, 'linestyle', ':', 'color', [0.3 0.3 1])
line([Lratio_exampleSession_median Lratio_exampleSession_median], [yl(1) 0.5], 'linewidth', 0.5, 'linestyle', ':', 'color', [0.3 0.3 1])
text(Lratio_exampleSession_median, 0.05, sprintf('%.1e', Lratio_exampleSession_median), 'fontsize', fontsize, 'color', [0.3 0.3 1])


line([xl(1) Lratio_pooled_median], [0.5 0.5], 'linewidth', 0.5, 'linestyle', ':', 'color', [1 0.3 0.3])
line([Lratio_pooled_median Lratio_pooled_median], [yl(1) 0.5], 'linewidth', 0.5, 'linestyle', ':', 'color', [1 0.3 0.3])
text(Lratio_exampleSession_median, 0.12, sprintf('%.1e', Lratio_pooled_median), 'fontsize', fontsize, 'color', [1 0.3 0.3])


set(gca,'XScale','log', 'YTick', 0:0.2:1, 'fontsize', 6, 'linewidth', 0.5, 'TickDir', 'out', 'TickLength',[0.01, 0.01])
% set(gca, 'XTick', 10.^linspace(-20, 10, 5))


xtls = get(gca, 'XTickLabel');
xtls{1} = ['<' xtls{1}];

set(gca, 'XTickLabel', xtls)


grid on


xlabel('L-ratio(pairwise)', 'fontsize', fontsize)
ylabel('cumulative probability', 'fontsize', fontsize)

legend([h1, h2], 'example session', 'all sessions', 'location', 'northwest', 'box', 'off', 'fontsize', fontsize)




%% plot the distribution of cross-corrrelation for individual sessions

plotwidth  = 110;
plotheight = 100;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);


hold on
set(gca, 'YScale', 'log')

        
for iSess = 1:nSessions

    if iSess == 6
        cl = [0.3 0.3 1];
    else
        cl = 'k';
    end

    curr_x = (iSess-ceil(nSessions/2))*0.85/nSessions;

    med = nanmedian(Lratio_allSessions{iSess});
    lq  = nanprctile(Lratio_allSessions{iSess}, 25);
    hq  = nanprctile(Lratio_allSessions{iSess}, 75);
    data_min = min(Lratio_allSessions{iSess}); data_min = max(data_min, 1e-20);
    data_max = max(Lratio_allSessions{iSess});


    % whiskers
    line([curr_x curr_x], [data_min lq], 'linewidth', 0.25, 'color', cl)
    line([curr_x curr_x], [data_max hq], 'linewidth', 0.25, 'color', cl)

    % box
    patch([curr_x-0.02 curr_x+0.02 curr_x+0.02 curr_x-0.02], [lq lq hq hq], cl, 'EdgeColor', 'none', 'FaceAlpha', 0.8)

    % median
    scatter(curr_x, med, 5, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'none')

end



pooledData = cell2mat(Lratio_allSessions);

med = nanmedian(pooledData);
lq  = nanprctile(pooledData, 25);
hq  = nanprctile(pooledData, 75);
data_min = min(pooledData); data_min = max(data_min, 1e-20);
data_max = max(pooledData);

curr_x = (iSess+1-ceil(nSessions/2))*0.85/nSessions;


% whiskers
line([curr_x curr_x], [data_min lq], 'linewidth', 0.25, 'color', [1 0.3 0.3])
line([curr_x curr_x], [data_max hq], 'linewidth', 0.25, 'color', [1 0.3 0.3])

% box
patch([curr_x-0.02 curr_x+0.02 curr_x+0.02 curr_x-0.02], [lq lq hq hq], [1 0.3 0.3], 'EdgeColor', 'none', 'FaceAlpha', 0.8)

% median
scatter(curr_x, med, 5, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'none')



set(gca , 'YScale', 'log', 'box', 'off', 'linewidth', 1, 'fontsize', fontsize, 'TickDir', 'out','TickLength',[0.01, 0.01])
ylabel('L-ratio(pairwise)', 'fontsize', fontsize)
ax = gca;
ax.YGrid = 'on';
ax.XAxis.Visible = 'off';

% ytls = get(gca, 'YTickLabel');
% ytls{1} = ['<' ytls{1}];

% set(gca, 'YTickLabel', ytls)



%% functions

function output = nanprctile(data, percentile)

data = data(~isnan(data));
output = prctile(data, percentile);

end

