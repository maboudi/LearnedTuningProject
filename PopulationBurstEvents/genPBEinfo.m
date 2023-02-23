
clear 
clc

currDir = '~/Documents/HMM_project/GreatLakes_March2021_addingMoreShuffles_CircularPF_CircularPosterior';
cd(currDir)
rr = dir;


for sessionNumber = 1:12

    sessionNumber
    
    
sessionName = rr(sessionNumber+2).name;

currSessionFolder = fullfile(currDir, sessionName);

cd(currSessionFolder)


load('Bayesian.mat', 'behavior', 'secondaryPBEs', 'RUNbinnedPBEs', 'PREbinnedPBEs', 'POSTbinnedPBEs', 'begPosition', 'endPosition', 'BDseqscore', '')


currPBEinfo = struct('sessionName', [], 'PBEno', [], ...
                    'startT', [], 'endT', [], 'peakT', [], 'peakMUA', [], ...
                    'nFiringUnits', [], 'duration', [], ...
                    'rippleOverlap', [], ...
                    'NREM', [], 'REM', [], 'QW', [], 'AW', [], ...
                    'PRE', [], 'RUN', [], 'POST', [], ...
                    'firing_1msbin', [], 'firing_20msbin', [], ...
                    'posteriorProbMat', [], ...
                    'rt_ts', [], 'rt_ui', [], 'rt_pf', [], 'rt_ds', [], ...
                    'wc_ts', [], 'wc_ui', [], 'wc_pf', [], 'wc_ds', [], ...
                    'begPosition', [], 'endPosition', [], ...
                    'radonIntegral', [], 'weightedCorr', []);
                
                
PREidx = find(secondaryPBEs(:, 1) > behavior.time(1,1) & secondaryPBEs(:, 2) < behavior.time(1,2));
RUNidx = find(secondaryPBEs(:, 1) > behavior.time(2,1) & secondaryPBEs(:, 2) < behavior.time(2,2));
POSTidx = find(secondaryPBEs(:, 1) > behavior.time(3,1) & secondaryPBEs(:, 2) < behavior.time(3,2));
idx = [PREidx; RUNidx; POSTidx];



% % PRE
                
nPBEs_PRE = size(PREbinnedPBEs, 1);

for ipbe = 1:nPBEs_PRE
    
    PBEno = PREidx(ipbe);
    
    currPBEinfo(PBEno).sessionName = sessionName;
    
    currPBEinfo(PBEno).PBEno   = PBEno;
    currPBEinfo(PBEno).startT  = secondaryPBEs(PBEno, 1);
    currPBEinfo(PBEno).endT    = secondaryPBEs(PBEno, 2);
    currPBEinfo(PBEno).peakT   = secondaryPBEs(PBEno, 3);
    currPBEinfo(PBEno).peakMUA = secondaryPBEs(PBEno, 4);
    
    currPBEinfo(PBEno).rippleOverlap = secondaryPBEs(PBEno, 5);
    
    currPBEinfo(PBEno).NREM    = secondaryPBEs(PBEno, 6);
    currPBEinfo(PBEno).REM     = secondaryPBEs(PBEno, 7);
    currPBEinfo(PBEno).QW      = secondaryPBEs(PBEno, 8);
    currPBEinfo(PBEno).AW      = secondaryPBEs(PBEno, 9);
    
    currPBEinfo(PBEno).PRE     = 1;
    currPBEinfo(PBEno).RUN     = 0;
    currPBEinfo(PBEno).POST    = 0;
    
    currPBEinfo(PBEno).firing_1msbin  = PREbinnedPBEs{ipbe, 1};
    currPBEinfo(PBEno).firing_20msbin = PREbinnedPBEs{ipbe, 2};
    
    currPBEinfo(PBEno).nFiringUnits   = numel(find(sum(PREbinnedPBEs{ipbe, 2}, 2)));
    currPBEinfo(PBEno).duration       = size(PREbinnedPBEs{ipbe, 1}, 2);
    
    
    currPBEinfo(PBEno).posteriorProbMat = posteriorProbMatrix.PRE.data{ipbe};
    
    currPBEinfo(PBEno).rt_ts   = BDseqscore.PRE.data.wPBEtimeswap.RadonTF.prctilescore(ipbe);
    currPBEinfo(PBEno).rt_ui   = BDseqscore.PRE.data.unitIDshuffle.RadonTF.prctilescore(ipbe);
    currPBEinfo(PBEno).rt_pf   = BDseqscore.PRE.data.circularPFshuffle.RadonTF.prctilescore(ipbe);
    currPBEinfo(PBEno).rt_ds   = BDseqscore.PRE.data.circularPosteriorShuffle.RadonTF.prctilescore(ipbe);
    
    currPBEinfo(PBEno).wc_ts   = BDseqscore.PRE.data.wPBEtimeswap.weightedCorr.prctilescore(ipbe);
    currPBEinfo(PBEno).wc_ui   = BDseqscore.PRE.data.unitIDshuffle.weightedCorr.prctilescore(ipbe);
    currPBEinfo(PBEno).wc_pf   = BDseqscore.PRE.data.circularPFshuffle.weightedCorr.prctilescore(ipbe);
    currPBEinfo(PBEno).wc_ds   = BDseqscore.PRE.data.circularPosteriorShuffle.weightedCorr.prctilescore(ipbe);
    
    
    currPBEinfo(PBEno).begPosition = begPosition.PRE.data(ipbe, :);
    currPBEinfo(PBEno).endPosition = endPosition.PRE.data(ipbe, :);

    currPBEinfo(PBEno).radonIntegral = RadonTF.PRE.data(ipbe);
    currPBEinfo(PBEno).weightedCorr  = weightedCorr.PRE.data(ipbe);
    
end



% % RUN

nPBEs_RUN = size(RUNbinnedPBEs, 1);

for ipbe = 1:nPBEs_RUN
    
    PBEno = nPBEs_PRE+ipbe;
    
    currPBEinfo(PBEno).sessionName = sessionName;
    
    currPBEinfo(PBEno).PBEno   = PBEno;
    
    currPBEinfo(PBEno).startT  = secondaryPBEs(PBEno, 1);
    currPBEinfo(PBEno).endT    = secondaryPBEs(PBEno, 2);
    currPBEinfo(PBEno).peakT   = secondaryPBEs(PBEno, 3);
    currPBEinfo(PBEno).peakMUA = secondaryPBEs(PBEno, 4);
    
    currPBEinfo(PBEno).rippleOverlap = secondaryPBEs(PBEno, 5);
    
    currPBEinfo(PBEno).NREM    = secondaryPBEs(PBEno, 6);
    currPBEinfo(PBEno).REM     = secondaryPBEs(PBEno, 7);
    currPBEinfo(PBEno).QW      = secondaryPBEs(PBEno, 8);
    currPBEinfo(PBEno).AW      = secondaryPBEs(PBEno, 9);
    
    currPBEinfo(PBEno).PRE     = 0;
    currPBEinfo(PBEno).RUN     = 1;
    currPBEinfo(PBEno).POST    = 0;
    
    currPBEinfo(PBEno).firing_1msbin  = RUNbinnedPBEs{ipbe, 1};
    currPBEinfo(PBEno).firing_20msbin = RUNbinnedPBEs{ipbe, 2};
    
    currPBEinfo(PBEno).nFiringUnits   = numel(find(sum(RUNbinnedPBEs{ipbe, 2}, 2)));
    currPBEinfo(PBEno).duration       = size(RUNbinnedPBEs{ipbe, 1}, 2);
    
    currPBEinfo(PBEno).posteriorProbMat = posteriorProbMatrix.RUN.data{ipbe};
    
    currPBEinfo(PBEno).rt_ts   = BDseqscore.RUN.data.wPBEtimeswap.RadonTF.prctilescore(ipbe);
    currPBEinfo(PBEno).rt_ui   = BDseqscore.RUN.data.unitIDshuffle.RadonTF.prctilescore(ipbe);
    currPBEinfo(PBEno).rt_pf   = BDseqscore.RUN.data.circularPFshuffle.RadonTF.prctilescore(ipbe);
    currPBEinfo(PBEno).rt_ds   = BDseqscore.RUN.data.circularPosteriorShuffle.RadonTF.prctilescore(ipbe);
    
    currPBEinfo(PBEno).wc_ts   = BDseqscore.RUN.data.wPBEtimeswap.weightedCorr.prctilescore(ipbe);
    currPBEinfo(PBEno).wc_ui   = BDseqscore.RUN.data.unitIDshuffle.weightedCorr.prctilescore(ipbe);
    currPBEinfo(PBEno).wc_pf   = BDseqscore.RUN.data.circularPFshuffle.weightedCorr.prctilescore(ipbe);
    currPBEinfo(PBEno).wc_ds   = BDseqscore.RUN.data.circularPosteriorShuffle.weightedCorr.prctilescore(ipbe);
    
    
    currPBEinfo(PBEno).begPosition = begPosition.RUN.data(ipbe, :);
    currPBEinfo(PBEno).endPosition = endPosition.RUN.data(ipbe, :);

    currPBEinfo(PBEno).radonIntegral = RadonTF.RUN.data(ipbe);
    currPBEinfo(PBEno).weightedCorr  = weightedCorr.RUN.data(ipbe);
    
end



% % POST

nPBEs_POST = size(POSTbinnedPBEs, 1);

for ipbe = 1:nPBEs_POST
    
    PBEno = nPBEs_PRE+nPBEs_RUN+ipbe;
    
    currPBEinfo(PBEno).sessionName = sessionName;
    
    currPBEinfo(PBEno).PBEno   = PBEno;
    
    currPBEinfo(PBEno).startT  = secondaryPBEs(PBEno, 1);
    currPBEinfo(PBEno).endT    = secondaryPBEs(PBEno, 2);
    currPBEinfo(PBEno).peakT   = secondaryPBEs(PBEno, 3);
    currPBEinfo(PBEno).peakMUA = secondaryPBEs(PBEno, 4);
    
    currPBEinfo(PBEno).rippleOverlap = secondaryPBEs(PBEno, 5);
    
    currPBEinfo(PBEno).NREM    = secondaryPBEs(PBEno, 6);
    currPBEinfo(PBEno).REM     = secondaryPBEs(PBEno, 7);
    currPBEinfo(PBEno).QW      = secondaryPBEs(PBEno, 8);
    currPBEinfo(PBEno).AW      = secondaryPBEs(PBEno, 9);
    
    currPBEinfo(PBEno).PRE     = 0;
    currPBEinfo(PBEno).RUN     = 0;
    currPBEinfo(PBEno).POST    = 1;
    
    currPBEinfo(PBEno).firing_1msbin  = POSTbinnedPBEs{ipbe, 1};
    currPBEinfo(PBEno).firing_20msbin = POSTbinnedPBEs{ipbe, 2};
    
    currPBEinfo(PBEno).nFiringUnits   = numel(find(sum(POSTbinnedPBEs{ipbe, 2}, 2)));
    currPBEinfo(PBEno).duration       = size(POSTbinnedPBEs{ipbe, 1}, 2);
    
    currPBEinfo(PBEno).posteriorProbMat = posteriorProbMatrix.POST.data{ipbe};
    
    currPBEinfo(PBEno).rt_ts   = BDseqscore.POST.data.wPBEtimeswap.RadonTF.prctilescore(ipbe);
    currPBEinfo(PBEno).rt_ui   = BDseqscore.POST.data.unitIDshuffle.RadonTF.prctilescore(ipbe);
    currPBEinfo(PBEno).rt_pf   = BDseqscore.POST.data.circularPFshuffle.RadonTF.prctilescore(ipbe);
    currPBEinfo(PBEno).rt_ds   = BDseqscore.POST.data.circularPosteriorShuffle.RadonTF.prctilescore(ipbe);
    
    currPBEinfo(PBEno).wc_ts   = BDseqscore.POST.data.wPBEtimeswap.weightedCorr.prctilescore(ipbe);
    currPBEinfo(PBEno).wc_ui   = BDseqscore.POST.data.unitIDshuffle.weightedCorr.prctilescore(ipbe);
    currPBEinfo(PBEno).wc_pf   = BDseqscore.POST.data.circularPFshuffle.weightedCorr.prctilescore(ipbe);
    currPBEinfo(PBEno).wc_ds   = BDseqscore.POST.data.circularPosteriorShuffle.weightedCorr.prctilescore(ipbe);
    
    
    currPBEinfo(PBEno).begPosition = begPosition.POST.data(ipbe, :);
    currPBEinfo(PBEno).endPosition = endPosition.POST.data(ipbe, :);

    currPBEinfo(PBEno).radonIntegral = RadonTF.POST.data(ipbe);
    currPBEinfo(PBEno).weightedCorr  = weightedCorr.POST.data(ipbe);
    
end

sessionName(sessionName == '-') = '_';

PBEinfo_pooled.(sessionName) = currPBEinfo;


end