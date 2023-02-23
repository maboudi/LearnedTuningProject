

detectParams.time_resolution    = 0.001;
detectParams.exclude            = [];
detectParams.smoothingSigma     = 0.01;
detectParams.smoothinghalfwidth = 0.03;
detectParams.minDur             = 0.04; % number of bins
detectParams.maxDur             = 0.5;
detectParams.maxVelocity        = 5; 
detectParams.thRatioThresh      = 1;
detectParams.rippleAThresh      = 3;



thresholds = [1 2 3];
detectionUnits = {'MUA'; 'PYR_INT';'PYR'};
rippleAThreshs = [2 5]; 


rippleOverlap_ratio_of_ripples = zeros(numel(thresholds), numel(detectionUnits), numel(rippleAThreshs));
rippleOverlap_ratio_of_PBEs    = zeros(numel(thresholds), numel(detectionUnits), numel(rippleAThreshs));


for ithresh = 1:numel(thresholds)
    ithresh
    for idu = 1:numel(detectionUnits)


        detectParams.threshold = thresholds(ithresh);
        detectParams.MUAorPyr  = detectionUnits{idu};

        
        if strcmp(detectParams.MUAorPyr, 'MUA') 
            [primaryPBEs, sdat] = PBPeriods(MUA, detectParams, fileInfo);
        elseif strcmp(detectParams.MUAorPyr, 'PYR_INT')
            [primaryPBEs, sdat] = PBPeriods(PYR_INT, detectParams, fileInfo);
        elseif strcmp(detectParams.MUAorPyr, 'PYR')
            [primaryPBEs, sdat] = PBPeriods(PYR, detectParams, fileInfo);
        end


        nPBEs = size(primaryPBEs, 1);
        PBEInfo = struct('sessionName', [], 'PBEno', [], ...
                            'startT', [], 'endT', [], 'peakT', [], ...
                            'startT_adj', [], 'endT_adj', [], 'peakT_adj', [], ...
                            'peakMUA', [], 'MUAwave', [], ...
                            'nFiringUnits', [], 'duration', [], ...
                            'peakRippleA', [], 'rippleWave', [], ...
                            'thetaRatio', [], 'SWA', [], 'emg', [], ...
                            'epoch', [], 'brainState', [], ...
                            'fr_1msbin', [], 'fr_20msbin', []); 
        
        
        binDur = 0.02; % in sec

        PBEs_all = [];

        for ipbe = 1:nPBEs

            PBE_primary = primaryPBEs(ipbe, :);
            
            [~, nFiringUnits] = timeBinning(PBE_primary, spikes_pyr, binDur, fileInfo);


            PBEs_all = [PBEs_all; PBE_primary];
                     
            PBEInfo(ipbe).startT      = PBE_primary(1);
            PBEInfo(ipbe).endT        = PBE_primary(2);
            PBEInfo(ipbe).peakT       = PBE_primary(3);
            
            PBEInfo(ipbe).duration     = PBE_primary(2) - PBE_primary(1);  
            PBEInfo(ipbe).nFiringUnits = nFiringUnits;

            PBEInfo(ipbe).peakRippleA = max(ripplePower_adjust(floor(PBE_primary(1)*1e3): floor(PBE_primary(2)*1e3)));

            PBEInfo(ipbe).thetaRatio  = interp1(timePnts, theratio, PBE_primary(3));
            PBEInfo(ipbe).SWA         = interp1(sw_timePnts, slowWave, PBE_primary(3)); % what exactly is the SWA?
            PBEInfo(ipbe).emg         = interp1(sw_timePnts, emg, PBE_primary(3));

            PBEInfo(ipbe).velocity    = interp1(fileInfo.speed.t, fileInfo.speed.v, PBE_primary(3));
            PBEInfo(ipbe).velocity(isnan(PBEInfo(ipbe).velocity)) = -1;
            
            
            if PBEInfo(ipbe).startT > behavior.time(1,1) && PBEInfo(ipbe).startT < behavior.time(1,2)
               PBEInfo(ipbe).epoch = 'pre';
            elseif PBEInfo(ipbe).startT > behavior.time(2,1) && PBEInfo(ipbe).startT < behavior.time(2,2)
               PBEInfo(ipbe).epoch = 'run';
            elseif PBEInfo(ipbe).startT > behavior.time(3,1) && PBEInfo(ipbe).startT < behavior.time(3,2)
               PBEInfo(ipbe).epoch = 'post';
            end

            
            if PBEInfo(ipbe).SWA > swthresh && ismember(PBEInfo(ipbe).epoch, {'pre'; 'post'})
                PBEInfo(ipbe).brainState = 'NREM';
            elseif PBEInfo(ipbe).thetaRatio < detectParams.thRatioThresh || (PBEInfo(ipbe).peakRippleA > detectParams.rippleAThresh && PBEInfo(ipbe).peakRippleA < 15) 
                PBEInfo(ipbe).brainState = 'QW';
            elseif PBEInfo(ipbe).emg < emgThresh
                PBEInfo(ipbe).brainState = 'REM';
            else
                PBEInfo(ipbe).brainState = 'WAKE';
            end

            
            if ipbe == 1
                fprintf('\n dividing the PBEs into 20 ms time bins ..')
            end

            
            if mod(ipbe, 1000) == 0
                fprintf('.')
            end

        end
               
        PBEs_accepted = PBEs_all(ismember({PBEInfo.brainState}, {'NREM'; 'QW'}) & ...
                                [PBEInfo.duration] >= detectParams.minDur & [PBEInfo.duration] <= detectParams.maxDur & ...
                                [PBEInfo.velocity] < detectParams.maxVelocity, :);


        % find the overlap between ripple events and detected PBEs using any of the
        % methods

        nPBE = size(PBEs_accepted, 1);

        for irt = 1:numel(rippleAThreshs)

            currentRippleThresh = rippleAThreshs(irt);

            highAmpRipples = rippleInfo([rippleInfo.peakRippleA] > currentRippleThresh);

            nRipples = numel(highAmpRipples);


            rippleIncluded = zeros(nRipples, 1);
            for ir = 1:nRipples 

                temp = find(PBEs_accepted(:,1) < highAmpRipples(ir).peakT & PBEs_accepted(:,2) > highAmpRipples(ir).peakT, 1, 'first');

                if ~isempty(temp)
                    rippleIncluded(ir) = 1;
                end

            end

            rippleOverlap_ratio_of_ripples(ithresh, idu, irt) = length(find(rippleIncluded))/nRipples;



            highAmpRipples_peakT = [highAmpRipples.peakT];

            PBEoverlapped = zeros(nPBE, 1);
            for ip = 1:nPBE 

                temp = find(highAmpRipples_peakT > PBEs_accepted(ip, 1) & highAmpRipples_peakT < PBEs_accepted(ip, 2), 1, 'first');

                if ~isempty(temp)
                    PBEoverlapped(ip) = 1;
                end

            end

            rippleOverlap_ratio_of_PBEs(ithresh, idu, irt)    = length(find(PBEoverlapped))/size(PBEs_accepted, 1);

        end
    end
end


folderName = fullfile(storageDir,  'PopulationBurstEvents');
if ~exist(folderName, 'dir'); mkdir(folderName); end


fileName = [sessionName '_rippleOverlap_fraction_of_ripples.mat'];
save(fullfile(folderName, fileName), 'rippleOverlap_ratio_of_ripples', 'thresholds', 'detectionUnits', 'rippleAThreshs')

fileName = [sessionName '_rippleOverlap_fraction_of_PBEs.mat'];
save(fullfile(folderName, fileName), 'rippleOverlap_ratio_of_PBEs', 'thresholds', 'detectionUnits', 'rippleAThreshs')

clear detectParams

figure; 


% ripple overlap as a ratio of total number of ripples
% all ripples: ripple power > 2
subplot(2,2, 1)

plot(rippleOverlap_ratio_of_ripples(:, 1, 1), 'r', 'marker', 's', 'linewidth', 2, 'markersize', 5)
hold on
plot(rippleOverlap_ratio_of_ripples(:, 2, 1), 'g',  'marker', 's', 'linewidth', 2, 'markersize', 5)
hold on
plot(rippleOverlap_ratio_of_ripples(:, 3, 1), 'b',  'marker', 's', 'linewidth', 2, 'markersize', 5)

xlabel('PBE thresholds(z)')
ylabel({'fraction of ripples'; 'overlapping w PBEs'})
title('all ripples(>2z)')
xlim([0 4])
xticks([1 2 3])

legend('MUA', 'pyr+int', 'only pyr', 'location', 'best')


% high amplitude ripples : ripple power > 5
subplot(2,2, 3)

plot(rippleOverlap_ratio_of_ripples(:, 1, 2), 'r',  'marker', 's', 'linewidth', 2, 'markersize', 5)
hold on
plot(rippleOverlap_ratio_of_ripples(:, 2, 2), 'g',  'marker', 's', 'linewidth', 2, 'markersize', 5)
hold on
plot(rippleOverlap_ratio_of_ripples(:, 3, 2), 'b',  'marker', 's', 'linewidth', 2, 'markersize', 5)

xlabel('PBE thresholds(z)')
ylabel({'fraction of ripples'; 'overlapping w PBEs'})
title('stronger ripples(>5z)')
xlim([0 4])
xticks([1 2 3])


% ripple overlap as a ratio of total number of PBEs
% ripple power  > 2
subplot(2,2, 2)

plot(rippleOverlap_ratio_of_PBEs(:, 1, 1), 'r',  'marker', 's', 'linewidth', 2, 'markersize', 5)
hold on
plot(rippleOverlap_ratio_of_PBEs(:, 2, 1), 'g', 'marker', 's', 'linewidth', 2, 'markersize', 5)
hold on
plot(rippleOverlap_ratio_of_PBEs(:, 3, 1), 'b',  'marker', 's', 'linewidth', 2, 'markersize', 5)


xlabel('PBE thresholds(z)')
ylabel({'fraction of PBEs'; 'overlapping w ripples'})
title('all ripples(>2z)')
xlim([0 4])
xticks([1 2 3])


%ripple power > 5
subplot(2,2, 4)

plot(rippleOverlap_ratio_of_PBEs(:, 1, 2), 'r',  'marker', 's', 'linewidth', 2, 'markersize', 5)
hold on
plot(rippleOverlap_ratio_of_PBEs(:, 2, 2), 'g',  'marker', 's', 'linewidth', 2, 'markersize', 5)
hold on
plot(rippleOverlap_ratio_of_PBEs(:, 3, 2), 'b',  'marker', 's', 'linewidth', 2, 'markersize', 5)

xlabel('PBE thresholds(z)')
ylabel({'fraction of PBEs'; 'overlapping w ripples'})
title('stronger ripples(>5z)')
xlim([0 4])
xticks([1 2 3])


suptitle(sessionName)

folderName = fullfile(storageDir,  'PopulationBurstEvents');
if ~exist(folderName, 'dir'); mkdir(folderName); end


fileName = [sessionName '_rippleOverlap_MUAvsPyr'];
savepdf(gcf, fullfile(folderName, fileName), '-dpdf')

