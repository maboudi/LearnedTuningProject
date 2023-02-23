% clear; clc;

% currDir = '/Volumes/TempHDD/';
% DatasetName = 'CA1Sleep+Acute-NintendoDS-20170125-01-nlx';
% 
% FileBase = [currDir DatasetName '/' DatasetName];
% 
% varList = {'spikes';... %% spiking time stamps and quality of each unit (1 -> excitatory, 8 -> interneurons, 3 -> noisy, 4 -> MUA)
%             'ure';... %% 
%             'basics'...
%             };
% 
% for var = 1: length(varList)
%     load([FileBase '-' varList{var}])
% end
% 
% 
% Fs = 1e6;
% lfpSampleRate = basics.lfpSampleRate;
% nCh = basics.nChannels;
% t0 = basics.period.time(1);
% tend = basics.period.time(end);
% 
% 
% fileinfo.name = DatasetName;
% fileinfo.CA = [1 1 1 1];
% fileinfo.CA1thetach = [];
% fileinfo.badch = [];
% 
% %%% finding the best ripple channels
% 
% [best_channels,ch_medians] = BigRippleChannels(fileinfo,1);
% 
% 
% 
% beforeLaserPrd = 45; %min
% beforeLaserPrd = beforeLaserPrd *60 * 1e6;
% 
% firstLaserTime = evts_urethane.time(1,1);
% analysisPrd = [firstLaserTime-beforeLaserPrd firstLaserTime]; % firstLaserTime


[SWR, HFE]=SWRdetectionHilbert(lfpSampleRate,nCh,FileBase,best_channels, (analysisPrd(1)-t0)/1e6, (analysisPrd(2)-t0)/1e6);

