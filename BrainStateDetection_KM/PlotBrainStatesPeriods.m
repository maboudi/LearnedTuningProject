

figure;
hold on

for ii = 1:length(behavior.Wake)
    h1 = patch([behavior.Wake(ii,1) behavior.Wake(ii,2) behavior.Wake(ii,2) behavior.Wake(ii,1)], [0 0 110 110], ...
    'k', 'EdgeColor', 'none');
    set(h1,'FaceAlpha', 0.25)
end

for ii = 1:length(behavior.Drowsy)
    h2 = patch([behavior.Drowsy(ii,1) behavior.Drowsy(ii,2) behavior.Drowsy(ii,2) behavior.Drowsy(ii,1)], [0 0 110 110], ...
    'c', 'EdgeColor', 'none');
    set(h2,'FaceAlpha', 0.25)
end

for ii = 1:length(behavior.NREM)
    h3 = patch([behavior.NREM(ii,1) behavior.NREM(ii,2) behavior.NREM(ii,2) behavior.NREM(ii,1)], [0 0 110 110], ...
    'b', 'EdgeColor', 'none');
    set(h3,'FaceAlpha', 0.25)
end

for ii = 1:length(behavior.Intermediate)
    h4 = patch([behavior.Intermediate(ii,1) behavior.Intermediate(ii,2) behavior.Intermediate(ii,2) behavior.Intermediate(ii,1)], [0 0 110 110], ...
    'y', 'EdgeColor', 'none');
    set(h4,'FaceAlpha', 0.25)
end

for ii = 1:length(behavior.REM)
    h5 = patch([behavior.REM(ii,1) behavior.REM(ii,2) behavior.REM(ii,2) behavior.REM(ii,1)], [0 0 110 110], ...
    'g', 'EdgeColor', 'none');
    set(h5,'FaceAlpha', 0.25)
end

plot(speed.t, speed.v, 'linewidth', 1, 'color', 'k')

ylim([-0.5 110])

line(behavior.PREEpoch, [1 1], 'linewidth', 4, 'color', 'k')
line(behavior.POSTEpoch, [1 1], 'linewidth', 4, 'color', 'k')

legend([h1, h2, h3, h4, h5], 'Wake', 'Drowsy', 'NREM', 'Intemediate', 'REM')