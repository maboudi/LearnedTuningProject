clear; clc; close all

currDir = '/home/kouroshmaboudi/Documents/HMM project';
cd(currDir)


%% loading session data

rats = {'Roy','Ted', 'Kevin'};
allsessionNumbers = {[1], [1 2 3], 1};


cumPRE = [];
cumPREshuffle = [];

cumRUN = [];
cumRUNshuffle = [];

cumPOST = [];
cumPOSTshuffle = [];



for rr = 1:length(rats) % loop for all the rats and sessions
    
rat = rats{rr};
sessionNumbers = allsessionNumbers{rr};

    
for sessionNumber = sessionNumbers

% 
% rat = input('\n Enter the animal name (in quotes) \n');
% rat(1) = upper(rat(1));
% 
% 
% sessionNumber = input('\n Enter the session number \n');


sessionName = [rat 'Maze' num2str(sessionNumber)];
sessionName2 = [rat '-maze' num2str(sessionNumber)];

VarList = {'spikes','behavior','position','speed','basics','ripple'};

for var = 1 : length(VarList)
    load([currDir '/pooled_includingKevin/wake-' VarList{var} '.mat'])
end

spikes = eval(['spikes.' sessionName]);

behavior = eval(['behavior.' sessionName]);


if strcmp(rat, 'Kevin') % there are two wake periods, we need to concatenate them
   behavior.time(2,2) = behavior.time(3,2);
   behavior.time(3,1) = behavior.time(4,1);
   behavior.time(3,2) = behavior.time(4,2);
   behavior.time(4,:) = [];
end
    
position = eval(['position.' sessionName]);
speed = eval(['speed.' sessionName]);
basics = eval(['basics.' sessionName]);
ripple = eval(['ripple.' sessionName]);


%%



load(['/home/kouroshmaboudi/Documents/HMM project/Results_02_19_2018/' sessionName2 '/HMM/lsPFs/' sessionName2 '_decodingError.mat'])

cumPRE = [cumPRE; PosDecErrorPRE];
cumPREshuffle = [cumPREshuffle; PosDecErrorPRE_shuffledGamma];

cumRUN = [cumRUN; PosDecErrorRUN];
cumRUNshuffle = [cumRUNshuffle; PosDecErrorRUN_shuffledGamma];

cumPOST = [cumPOST; PosDecErrorPOST];
cumPOSTshuffle = [cumPOSTshuffle; PosDecErrorPOST_shuffledGamma];


end
end

load('/home/kouroshmaboudi/Documents/HMM project/Results_02_19_2018/Roy-maze1/HMM/lsPFs/Roy-maze1_decodingError.mat')

% plot decoding error

figure;
set(gcf, 'Units', 'pixels', 'position', [814 399 1648 606])

subplot(1,3,1)
hold on
plotErrorCdf(cumPRE, cumPREshuffle, PosDecErrorPRE, PosDecErrorPRE_shuffledGamma, 1)

title('PRE', 'fontsize', 16)


subplot(1,3,2)
plotErrorCdf(cumRUN, cumRUNshuffle, PosDecErrorRUN, PosDecErrorRUN_shuffledGamma, 0)
title('RUN', 'fontsize', 16)


subplot(1,3,3)
plotErrorCdf(cumPOST, cumPOSTshuffle, PosDecErrorPOST, PosDecErrorPOST_shuffledGamma, 0)
title('POST', 'fontsize', 16)


savepdf(gcf, '_decodingError')
print(gcf, '_decodingError', '-dpng')




%% functions

function plotErrorCdf(cum, cumShuffle, example, exampleshuffle, ylabelneeded)

hold on

% example

[h1, bins] = hist(example, 100);
h1c = cumsum(h1)/sum(h1);

exampleMed = median(example);
 
curve1 = plot(bins, h1c, 'linewidth', 4, 'color', [150, 150, 255]/255);
line([exampleMed exampleMed],[0 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', [150, 150, 255]/255)
line([0 exampleMed],[0.5 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', [150, 150, 255]/255)
text(exampleMed, 0.03, sprintf('%.1f cm', exampleMed), 'color', [150, 150, 255]/255, 'fontsize', 12, 'FontWeight', 'Bold')

% example shuffle

[h2, bins] = hist(exampleshuffle, 100);
h2c = cumsum(h2)/sum(h2);

exampleshuffleMed = median(exampleshuffle);
 
curve2 = plot(bins, h2c, 'linewidth', 4, 'color', 'k');
line([exampleshuffleMed exampleshuffleMed],[0 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', 'k')
line([0 exampleshuffleMed],[0.5 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', 'k')
text(exampleshuffleMed, 0.03, sprintf('%.1f cm', exampleshuffleMed), 'color', 'k', 'fontsize', 12, 'FontWeight', 'Bold')



% cumulative data

[h3, bins] = hist(cum, 100);
h3c = cumsum(h3)/sum(h3);

cumMed = median(cum);

curve3 = plot(bins, h3c, 'linewidth', 4, 'color', [160, 160, 160]/255);
line([cumMed cumMed],[0 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', [160, 160, 160]/255)
line([0 cumMed],[0.5 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', [160, 160, 160]/255)
text(cumMed, 0.1, sprintf('%.1f cm', cumMed), 'color', [160, 160, 160]/255, 'fontsize', 12, 'FontWeight', 'Bold')



% cumulative shuffle


[h4, bins] = hist(cumShuffle, 100);
h4c = cumsum(h4)/sum(h4);

cumShuffleMed = median(cumShuffle);

curve4 = plot(bins, h4c, 'linewidth', 4, 'color', 'r');
line([cumShuffleMed cumShuffleMed],[0 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', 'r')
line([0 cumShuffleMed],[0.5 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', 'r')
text(cumShuffleMed, 0.1, sprintf('%.1f cm', cumShuffle), 'color', 'r', 'fontsize', 12, 'FontWeight', 'Bold')





xlim([0 max([cum; cumShuffle; example; exampleshuffle])])
set(gca, 'fontsize', 14, 'linewidth', 3, 'TickDir', 'out', 'TickLength',[0.02, 0.01])


xlabel('decoding drror(cm)', 'fontsize', 16)

if ylabelneeded
    ylabel('cumulative ratio', 'fontsize', 24)
    
%     legend('actual', 'shuffle state probability', 'shuffle position index', 'Location', 'northwest'); 
    legend([curve1, curve2, curve3, curve4], 'example', 'example shuffle', 'all sessions', 'all sessions shuffle', 'Location', 'northwest'); 

    legend boxoff 
end

axis square
end




function savepdf(gcf, pdfname)

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

print(gcf, pdfname,'-depsc','-r0')

end






