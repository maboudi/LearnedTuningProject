

parentFolder = '/home/kouroshmaboudi/Documents/NCMLproject/MUAvsonlyPyr';


temp = dir(parentFolder);
sessionNames = {temp(3:end).name};
nSessions = numel(sessionNames);

fraction_of_ripples = cell(nSessions, 1);
fraction_of_PBEs    = cell(nSessions, 1);

for sess = 1:nSessions
    
    sessionName = sessionNames{sess};
    basePath = fullfile(parentFolder, sessionName);
    
    load(fullfile(basePath, 'PopulationBurstEvents', [sessionName '_rippleOverlap_fraction_of_PBEs.mat']), 'rippleOverlap_ratio_of_PBEs')
    load(fullfile(basePath, 'PopulationBurstEvents', [sessionName '_rippleOverlap_fraction_of_ripples.mat']), 'rippleOverlap_ratio_of_ripples')
    
    fraction_of_ripples{sess} = rippleOverlap_ratio_of_ripples;
    fraction_of_PBEs{sess} = rippleOverlap_ratio_of_PBEs;
    
end
    
  

%%

figure; 

tranL = 0.5;

% ripple overlap as a ratio of total number of ripples
% all ripples: ripple power > 2
subplot(2,2, 1)

tt = 0;
for sess = 1:numel(sessionNames)
    
    tt = (sess-1)*4;
    
    plot(tt+1:tt+3, fraction_of_ripples{sess}(:, 1, 1), 'color', [1 0 0 tranL], 'marker', 's', 'linewidth', 2, 'markersize', 5); 
    hold on
    plot(tt+1:tt+3, fraction_of_ripples{sess}(:, 2, 1), 'color', [0 1 0 tranL],  'marker', 's', 'linewidth', 2, 'markersize', 5);
    hold on
    plot(tt+1:tt+3, fraction_of_ripples{sess}(:, 3, 1), 'color', [0 0 1 tranL],  'marker', 's', 'linewidth', 2, 'markersize', 5);
    
    rr = sessionNames{sess};
    rr(rr == '_') = '-';
    text(tt+1, 1, rr)
end

set(gca, 'box', 'off')
xlabel('PBE thresholds(z)')
ylabel({'fraction of ripples'; 'overlapping w PBEs'})
title('all ripples(>2z)')
% xlim([0 4])
ylim([0 1])
xticks([1 2 3 5 6 7 9 10 11 13 14 15 17 18 19])
xticklabels({'1';'2';'3'; '1'; '2'; '3'; '1'; '2'; '3'; '1'; '2'; '3'; '1'; '2'; '3'})
legend('MUA', 'pyr+int', 'only pyr', 'location', 'best')


% high amplitude ripples : ripple power > 5
subplot(2,2, 3)

tt = 0;
for sess = 1:numel(sessionNames)
    
    tt = (sess-1)*4;
    
    plot(tt+1:tt+3,fraction_of_ripples{sess}(:, 1, 2), 'color', [1 0 0 tranL],  'marker', 's', 'linewidth', 2, 'markersize', 5)
    hold on
    plot(tt+1:tt+3,fraction_of_ripples{sess}(:, 2, 2), 'color', [0 1 0 tranL],  'marker', 's', 'linewidth', 2, 'markersize', 5)
    hold on
    plot(tt+1:tt+3,fraction_of_ripples{sess}(:, 3, 2), 'color', [0 0 1 tranL],  'marker', 's', 'linewidth', 2, 'markersize', 5)
    
    rr = sessionNames{sess};
    rr(rr == '_') = '-';
    text(tt+1, 1, rr)
end
  
set(gca, 'box', 'off')
xlabel('PBE thresholds(z)')
ylabel({'fraction of ripples'; 'overlapping w PBEs'})
title('stronger ripples(>5z)')
% xlim([0 4])
ylim([0 1])
xticks([1 2 3 5 6 7 9 10 11 13 14 15 17 18 19])
xticklabels({'1';'2';'3'; '1'; '2'; '3'; '1'; '2'; '3'; '1'; '2'; '3'; '1'; '2'; '3'})


% ripple overlap as a ratio of total number of PBEs
% ripple power  > 2
subplot(2,2, 2)

for sess = 1:numel(sessionNames)
    
    tt = (sess-1)*4;
    
    plot(tt+1:tt+3,fraction_of_PBEs{sess}(:, 1, 1), 'color', [1 0 0 tranL],  'marker', 's', 'linewidth', 2, 'markersize', 5)
    hold on
    plot(tt+1:tt+3,fraction_of_PBEs{sess}(:, 2, 1), 'color', [0 1 0 tranL], 'marker', 's', 'linewidth', 2, 'markersize', 5)
    hold on
    plot(tt+1:tt+3,fraction_of_PBEs{sess}(:, 3, 1), 'color', [0 0 1 tranL],  'marker', 's', 'linewidth', 2, 'markersize', 5)
    
    rr = sessionNames{sess};
    rr(rr == '_') = '-';
    text(tt+1, 1, rr)
end

set(gca, 'box', 'off')
xlabel('PBE thresholds(z)')
ylabel({'fraction of PBEs'; 'overlapping w ripples'})
title('all ripples(>2z)')
% xlim([0 4])
ylim([0 1])
xticks([1 2 3 5 6 7 9 10 11 13 14 15 17 18 19])
xticklabels({'1';'2';'3'; '1'; '2'; '3'; '1'; '2'; '3'; '1'; '2'; '3'; '1'; '2'; '3'})


%ripple power > 5
subplot(2,2, 4)

for sess = 1:numel(sessionNames)
    
    tt = (sess-1)*4;
    
    plot(tt+1:tt+3,fraction_of_PBEs{sess}(:, 1, 2), 'color', [1 0 0 tranL],  'marker', 's', 'linewidth', 2, 'markersize', 5)
    hold on
    plot(tt+1:tt+3,fraction_of_PBEs{sess}(:, 2, 2), 'color', [0 1 0 tranL],  'marker', 's', 'linewidth', 2, 'markersize', 5)
    hold on
    plot(tt+1:tt+3,fraction_of_PBEs{sess}(:, 3, 2), 'color', [0 0 1 tranL],  'marker', 's', 'linewidth', 2, 'markersize', 5)
    
    rr = sessionNames{sess};
    rr(rr == '_') = '-';
    text(tt+1, 1, rr)
end

set(gca, 'box', 'off')
xlabel('PBE thresholds(z)')
ylabel({'fraction of PBEs'; 'overlapping w ripples'})
title('stronger ripples(>5z)')
% xlim([0 4])
ylim([0 1])
xticks([1 2 3 5 6 7 9 10 11 13 14 15 17 18 19])
xticklabels({'1';'2';'3'; '1'; '2'; '3'; '1'; '2'; '3'; '1'; '2'; '3'; '1'; '2'; '3'})

