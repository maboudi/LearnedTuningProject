%% post period

postIdx = HMMprctile_post > 0.95;
prepostIdx = HMMprctile_pre_post > 0.95;
runpostIdx = HMMprctile_run_post > 0.95;

noncong = [length(find(~postIdx & ~prepostIdx & ~runpostIdx)); ...
           length(find(~postIdx & prepostIdx & ~runpostIdx)); ...
           length(find(~postIdx & ~prepostIdx & runpostIdx));...
           length(find(~postIdx & prepostIdx & runpostIdx))];

       
cong = [length(find(postIdx & ~prepostIdx & ~runpostIdx)); ...
           length(find(postIdx & prepostIdx & ~runpostIdx)); ...
           length(find(postIdx & ~prepostIdx & runpostIdx));...
           length(find(postIdx & prepostIdx & runpostIdx))];
       
       
noncong_cumcount = cumsum(noncong);
noncong_cumcount = noncong_cumcount (end:-1:1);

cong_cumcount = cumsum(cong);
cong_cumcount = cong_cumcount(end:-1:1);

thecolors = [100 255 180; 255 100 100; 50 150 255 ;255 255 150]/255;

% {'g', 'r', 'b', 'y'};

h = figure;
hold on
for ii = 1:4
    bar(1:2, [noncong_cumcount(ii) cong_cumcount(ii)], 'FaceColor', thecolors(ii, :))
end
hold off

legend('both pre and run', 'run excl', 'pre excl', 'nor pre nor run')
set(gca, 'xlim', [0.5 2.5],'XTick', [1 2], 'XTickLabel',  {'non-cong.', 'congruent'}, 'fontsize', 16)


xlabel('Congruence with Post BPE model', 'fontsize', 16)
ylabel('number of PBEs', 'fontsize', 20)
title('TED - maze 1', 'fontsize', 20)

print(gcf, ['Ted-maze1_' 'withinVSbetweenEvents'], '-dpng')


%% run period

runIdx = HMMprctile_run > 0.95;
prerunIdx = HMMprctile_pre_run > 0.95;

noncong = [length(find(~runIdx & ~prerunIdx));
           length(find(~runIdx & prerunIdx))];
       
cong = [length(find(runIdx & ~prerunIdx));
        length(find(runIdx & prerunIdx))];      

       
noncong_cumcount = cumsum(noncong);
noncong_cumcount = noncong_cumcount (end:-1:1);

cong_cumcount = cumsum(cong);
cong_cumcount = cong_cumcount(end:-1:1);

thecolors =  [50 150 255 ;255 255 150]/255;

% {'g', 'r', 'b', 'y'};

h = figure;
hold on
for ii = 1:2
    bar(1:2, [noncong_cumcount(ii) cong_cumcount(ii)], 'FaceColor', thecolors(ii, :))
end
hold off

legend('explained by pre', 'not explained by pre')
set(gca, 'xlim', [0.5 2.5],'XTick', [1 2], 'XTickLabel',  {'non-cong.', 'congruent'}, 'fontsize', 16)


xlabel('Congruence with Run BPE model', 'fontsize', 16)
ylabel('number of PBEs', 'fontsize', 20)
title('TED - maze 1', 'fontsize', 20)

print(gcf, ['Ted-maze1_' 'RunPBEs_crossTest_histogram'], '-dpng')


% plot 2D histogram

data = [HMMprctile_pre_run HMMprctile_run];
h = hist3(data, [10 10]);

figure;
imagesc(h); set(gca, 'YDir', 'normal'); 
colormap('hot')
% colormap(flipud(colormap('hot')))
set(gca, 'fontsize', 14)

set(gca, 'XTick', [0.5 10.5], 'XTickLabel', {'0', '100'}, 'YTick', [0.5 10.5], 'YTickLabel', {'0', '100'})
xlabel({'conguence of Run PBEs'; 'with Run model (percentile)'}, 'fontsize', 14)
ylabel({'conguence of Run PBEs'; 'with Pre model (percentile)'}, 'fontsize', 14)

title('TED - maze 1', 'fontsize', 20)
axis square

print(gcf, ['Ted-maze1_' 'RunPBEs_crossTest_2Dhistogram'], '-dpng')





