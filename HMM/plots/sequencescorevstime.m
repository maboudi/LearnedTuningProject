
clear; clc; close all


currDir = '/home/kouroshmaboudi/Documents/HMM_project';
cd(currDir)


%% loading session data

rats = {'Roy','Ted', 'Kevin'};
allsessionNumbers = {[1 ], [1 2 3], 1};


for rr = 1%:length(rats) % loop for all the rats and sessions
    
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



% lfpSampleRate = basics.lfpSampleRate;
lfpSampleRate = 1e6; % agian we use the time unit after conversion
Fs = 3600*1e6; % to hour
% 
% nCh = basics.nChannels;

fileinfo = struct('name', sessionName2, 'animal', rat, 'xyt', [], 'tbegin', [], 'tend', [], 'Fs', Fs, 'lfpSampleRate', lfpSampleRate, 'pix2cm', 0.3861); 
fileinfo.pix2cm = 0.3861; % in cm

 


FileBase = [currDir '/' fileinfo.name];
mkdir(FileBase)
cd(FileBase)

load(['/home/kouroshmaboudi/Documents/HMM_project/' sessionName2 '/HMM/sequenceDetection/congruence(wi_and_bw)/' sessionName2 '_congruence.mat'])




bvrTimeList = behavior.list(:,[1 2])/Fs;
bvrState = behavior.list(:,3);


figure;
set(gca, 'position', [995   783   796   538])
hold on

%% time courses during RUN (for PRE and RUN(CV))

subplot(2,1,1)

load(['/home/kouroshmaboudi/Documents/HMM_project/' sessionName2 '/PopulationBurstEvents/RUN/' sessionName2 '-PBEs.mat'])

tbeginRUN = behavior.time(2,1)/Fs;
tendRUN = behavior.time(2,2)/Fs;

PBEtimesRUN = secondaryPBEs(:,1)/Fs;



%QW periods during RUN

qw = bvrTimeList(bvrState == 3,:); % nrem=1, rem=2, quiet=3, wake=4
qwRUN = qw(qw(:,1) > tbeginRUN & qw(:,2) < tendRUN, :);
qwRUN = qwRUN - tbeginRUN;


PlotSeqFreqvsTime(PBEtimesRUN, HMMprctileRUN, HMMprctile_PRE_RUN, 'RUN', 'PRE', 'RUN', qwRUN, 12, tbeginRUN, tendRUN)



%% time courses during RUN (for PRE and RUN(CV))

subplot(2,1,2)

load(['/home/kouroshmaboudi/Documents/HMM_project/' sessionName2 '/PopulationBurstEvents/POST/' sessionName2 '-PBEs.mat'])

tbeginPOST = behavior.time(3,1)/Fs;
tendPOST = behavior.time(3,2)/Fs;

PBEtimesPOST = secondaryPBEs(:,1)/Fs;

%nrem periods during POST

nrem = bvrTimeList(ismember(bvrState, [1 3]),:); % nrem=1, rem=2, quiet=3, wake=4
nremPOST = nrem(nrem(:,1) > tbeginPOST & nrem(:,2) < tendPOST, :);
nremPOST = nremPOST - tbeginPOST;


PlotSeqFreqvsTime(PBEtimesPOST, HMMprctile_RUN_POST, HMMprctile_PRE_POST, 'RUN', 'PRE', 'POST', nremPOST, 12, tbeginPOST, tendPOST)


savepdf(gcf, [FileBase '/' fileinfo.name '_seqvstime'])
print(gcf, [FileBase '/' fileinfo.name  '_seqvstime'], '-dpng')



end
end



function PlotSeqFreqvsTime(PBEtimes, PBEscores1, PBEscores2, trainPeriod1, trainPeriod2, testPeriod, offPeriods, noEpochs, tbegin, tend)

% noEpochs = 12; % epochs with duration of ~ 15 min

epochs = linspace(tbegin, tend, noEpochs+1);
epochDur = epochs(2) - epochs(1);

epochCenters = epochs(1:end-1)+epochDur/2 - tbegin;

[~, epochIdx] = histc(PBEtimes, epochs);


% model #1
epochMean1 = zeros(noEpochs, 1);
epochErr1 = zeros(noEpochs, 1);

for epoch = 1 : noEpochs
    epochPrctls = PBEscores1(epochIdx == epoch);
    epochMean1(epoch) = mean(epochPrctls);
    epochErr1(epoch) = std(epochPrctls)/length(epochPrctls);
end





% model #2
epochMean2 = zeros(noEpochs, 1);
epochErr2 = zeros(noEpochs, 1);

for epoch = 1 : noEpochs
    epochPrctls = PBEscores2(epochIdx == epoch);
    epochMean2(epoch) = mean(epochPrctls);
    epochErr2(epoch) = std(epochPrctls)/length(epochPrctls);
end





hold on

minY = floor(min([epochMean1; epochMean2])*10)/10;
maxY = ceil(max([epochMean1; epochMean2])*10)/10;

for ii = 1:length(offPeriods)
    patch([offPeriods(ii,1) offPeriods(ii,2) offPeriods(ii,2) offPeriods(ii,1)], [minY minY maxY maxY], ...
           [204 229 255]/255, 'EdgeColor', 'none');
%         set(p,'FaceAlpha', 0.5)
end

h1 = errorbar(epochCenters, epochMean1, epochErr1, 'color', 'k', 'linewidth', 2);

h2 = errorbar(epochCenters, epochMean2, epochErr2, 'color', 'r', 'linewidth', 2);




legend([h1, h2], ['train on ' trainPeriod1], ['train on ' trainPeriod2], 'Location','southeast')
legend boxoff 

set(gca, 'fontsize', 12, 'linewidth', 2, 'TickDir', 'out', 'TickLength',[0.02, 0.01], 'box', 'off', 'Layer', 'Top')
xlim([0 3])
ylim()

xlabel(['Time in ' testPeriod '(hr)'], 'fontsize', 14)
ylabel('PBE sequence score', 'fontsize', 12)



end


function savepdf(gcf, pdfname)

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

print(gcf, pdfname,'-depsc','-r0')

end



