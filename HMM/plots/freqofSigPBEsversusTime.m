clear; clc; close all

currDir = '/home/kouroshmaboudi/Documents/HMM project';
cd(currDir)


%% loading session data

rats = {'Roy','Ted', 'Kevin'};
allsessionNumbers = {[1 2 3], [1 2 3], 1};

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

load(['/home/kouroshmaboudi/Documents/HMM project/Results_02_19_2018/' sessionName2 '/HMM/sequenceDetection/congruence(wi_and_bw)/' sessionName2 '_congruence.mat'])




bvrTimeList = behavior.list(:,[1 2])/Fs;
bvrState = behavior.list(:,3);


figure;
set(gcf, 'position', [ 1188 482 965 554])


%% time courses during RUN (for PRE and RUN(CV))

subplot(2,1,1)

load(['/home/kouroshmaboudi/Documents/HMM project/Results_02_19_2018/' sessionName2 '/PopulationBurstEvents/RUN/' sessionName2 '-PBEs.mat'])

tbeginRUN = behavior.time(2,1)/Fs;
tendRUN = behavior.time(2,2)/Fs;

PBEtimesRUN = secondaryPBEs(:,1)/Fs;



%QW periods during RUN

qw = bvrTimeList(bvrState == 3,:); % nrem=1, rem=2, quiet=3, wake=4
qwRUN = qw(qw(:,1) > tbeginRUN & qw(:,2) < tendRUN, :);
qwRUN = qwRUN - tbeginRUN;


PlotSeqFreqvsTime(PBEtimesRUN, HMMprctile_run, HMMprctile_pre_run, 'exc. RUN', 'exc. PRE', 'RUN', qwRUN, 12, tbeginRUN, tendRUN)



%% time courses during RUN (for PRE and RUN(CV))

subplot(2,1,2)

load(['/home/kouroshmaboudi/Documents/HMM project/Results_02_19_2018/' sessionName2 '/PopulationBurstEvents/POST/' sessionName2 '-PBEs.mat'])

tbeginPOST = behavior.time(3,1)/Fs;
tendPOST = behavior.time(3,2)/Fs;

PBEtimesPOST = secondaryPBEs(:,1)/Fs;

%nrem periods during RUN

nrem = bvrTimeList(bvrState == 1,:); % nrem=1, rem=2, quiet=3, wake=4
nremPOST = nrem(nrem(:,1) > tbeginPOST & nrem(:,2) < tendPOST, :);
nremPOST = nremPOST - tbeginPOST;


PlotSeqFreqvsTime(PBEtimesPOST, HMMprctile_run_post, HMMprctile_pre_post, 'exc. RUN', 'exc. PRE', 'POST', nremPOST, 12, tbeginPOST, tendPOST)


savepdf(gcf, [FileBase '/' fileinfo.name '_freqSigPBEsvstime'])
print(gcf, [FileBase '/' fileinfo.name  '_freqSigPBEsvstime'], '-dpng')



end
end



function PlotSeqFreqvsTime(PBEtimes, PBEscores1, PBEscores2, trainPeriod1, trainPeriod2, testPeriod, offPeriods, noEpochs, tbegin, tend)


hold on

epochs = linspace(tbegin, tend, noEpochs+1);
epochDur = epochs(2) - epochs(1);

epochCenters = epochs(1:end-1)+epochDur/2 - tbegin;




% model #1
count1 = zeros(noEpochs, 1);

PBEtimes1 = PBEtimes(PBEscores1 > 0.8 & PBEscores2 < 0.8);
PBEscores11 = PBEscores1(PBEscores1 > 0.8 & PBEscores2 < 0.8);


[~, epochIdx1] = histc(PBEtimes, epochs);
[~, epochIdx11] = histc(PBEtimes1, epochs);

for epoch = 1 : noEpochs
    epochPrctls = PBEscores11(epochIdx11 == epoch);
    count1(epoch) = length(epochPrctls);
    
    normTemp = length(PBEscores1(epochIdx1 == epoch));
    count1(epoch) = count1(epoch)/normTemp;
end

h1 = plot(epochCenters, count1, 'color', 'b', 'linewidth', 2);


% model #2
count2 = zeros(noEpochs, 1);


PBEtimes2 = PBEtimes(PBEscores2 > 0.8 & PBEscores1 < 0.8);
PBEscores22 = PBEscores2(PBEscores2 > 0.8 & PBEscores1 < 0.8);

[~, epochIdx2] = histc(PBEtimes, epochs);
[~, epochIdx22] = histc(PBEtimes2, epochs);

for epoch = 1 : noEpochs
    
    epochPrctls = PBEscores22(epochIdx22 == epoch);
    count2(epoch) = length(epochPrctls);
    
    normTemp = length(PBEscores1(epochIdx2 == epoch));
    count2(epoch) = count2(epoch)/normTemp;
    
end

h2 = plot(epochCenters, count2, 'color', 'r', 'linewidth', 2);

% h3 = plot(epochCenters, (count1 - count2), 'color', 'k', 'linewidth', 2);





% noEpochs = 12; % epochs with duration of ~ 15 min

ydim = ylim;

for ii = 1:length(offPeriods)
    p = patch([offPeriods(ii,1) offPeriods(ii,2) offPeriods(ii,2) offPeriods(ii,1)], [ydim(1) ydim(1) ydim(2) ydim(2)], ...
           [204 229 255]/255, 'EdgeColor', 'none');
%         set(p,'FaceAlpha', 0.5)
end




uistack(h2,'top')
uistack(h1,'top')


% legend([h1, h2, h3], ['train on ' trainPeriod1], ['train on ' trainPeriod2], [trainPeriod1 '/' trainPeriod2], 'Location','northeast')


% legend(h1, ['[' trainPeriod1 '-' trainPeriod2 ']* percentage of total'], 'Location','southeast')
% legend(h1, ['[' trainPeriod1  ']* percentage of total'], 'Location','southeast')

legend([h1, h2], ['congruent with' trainPeriod1], ['congruent with' trainPeriod2], 'Location','southeast')
legend boxoff 



set(gca, 'fontsize', 12, 'linewidth', 2, 'TickDir', 'out', 'TickLength',[0.02, 0.01], 'box', 'off', 'Layer', 'Top')
xlim([0 3])

xlabel(['Time in ' testPeriod '(hr)'], 'fontsize', 14)
ylabel('Ratio of significant PBEs', 'fontsize', 12)



end


function savepdf(gcf, pdfname)

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

print(gcf, pdfname,'-depsc','-r0')

end



