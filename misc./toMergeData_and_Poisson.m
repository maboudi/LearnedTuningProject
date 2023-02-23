
currDir = '/home/kouroshmaboudi/Documents/HMM_project/GreatLakes_firstRun_Nov2020/Grosmark_originalClusters'; 

cd(currDir)


rr = dir;

sessionNames = fieldnames(rr);
periodNames  = {'PRE';'RUN';'POST'};

for sessionNumber = [2 5 8]
    
    clearvars -except currDir sessionNames periodNames sessionNumber rr
    
    % loading poisson Data
    sessionName = rr(sessionNumber+2).name;
    
    currSessionFolder = fullfile(currDir, 'Poisson', sessionName, 'BayesianDecoding');
    ll = dir(currSessionFolder);
    
    load(fullfile(currSessionFolder, ll(3).name))
    
    
    for iperiod = 1:3
                
        currPeriod = periodNames{iperiod};
        
        BDseqscore_p.(currPeriod)          = BDseqscore.(currPeriod).p;
        begPosition_p.(currPeriod)         = begPosition.(currPeriod).p;
        endPosition_p.(currPeriod)         = endPosition.(currPeriod).p;
        coveredLen_p.(currPeriod)          = coveredLen.(currPeriod).p;
        coveredLen_null_p.(currPeriod)     = coveredLen_null.(currPeriod).p;
        jumpDist_p.(currPeriod)            = jumpDist.(currPeriod).p;
        onlyLineElements_p.(currPeriod)    = onlyLineElements.(currPeriod).p;
        posteriorProbMatrix_p.(currPeriod) = posteriorProbMatrix.(currPeriod).p;
        replayScore_p.(currPeriod)         = replayScore.(currPeriod).p;
        replayScore_null_p.(currPeriod)    = replayScore_null.(currPeriod).p;
        sigMatrix_p.(currPeriod)           = sigMatrix.(currPeriod).p;
        weightedCorr_p.(currPeriod)        = weightedCorr.(currPeriod).p;
        weightedCorr_null_p.(currPeriod)   = weightedCorr_null.(currPeriod).p;
        PREbinnedPBEs_p                    = PREbinnedPBEs;
        RUNbinnedPBEs_p                    = RUNbinnedPBEs;
        POSTbinnedPBEs_p                   = POSTbinnedPBEs;
        
    end
    
    clear BDseqscore begPosition endPosition coveredLen coveredLen_null jumpDist onlyLineElements posteriorProbMatrix ...
        replayScore replayScore_null sigMatrix weightedCorr weightedCorr_null PREbinnedPBEs RUNbinnedPBEs POSTbinnedPBEs
    
    
    currSessionFolder = fullfile(currDir, sessionName, 'BayesianDecoding');
    ll = dir(currSessionFolder);
    
    load(fullfile(currSessionFolder, ll(3).name))
    
    for iperiod = 1:3
                
        currPeriod = periodNames{iperiod};
        
        BDseqscore.(currPeriod).p          = BDseqscore_p.(currPeriod);
        begPosition.(currPeriod).p         = begPosition_p.(currPeriod);
        endPosition.(currPeriod).p         = endPosition_p.(currPeriod);
        coveredLen.(currPeriod).p          = coveredLen_p.(currPeriod);
        coveredLen_null.(currPeriod).p     = coveredLen_null_p.(currPeriod);
        jumpDist.(currPeriod).p            = jumpDist_p.(currPeriod);
        onlyLineElements.(currPeriod).p    = onlyLineElements_p.(currPeriod);
        posteriorProbMatrix.(currPeriod).p = posteriorProbMatrix_p.(currPeriod);
        replayScore.(currPeriod).p         = replayScore_p.(currPeriod);
        replayScore_null.(currPeriod).p    = replayScore_null_p.(currPeriod);
        sigMatrix.(currPeriod).p           = sigMatrix_p.(currPeriod);
        weightedCorr.(currPeriod).p        = weightedCorr_p.(currPeriod);
        weightedCorr_null.(currPeriod).p   = weightedCorr_null_p.(currPeriod);
        
    end
    
    
    save(fullfile(currSessionFolder, ['Bayesian_merged.mat']), ...
    'POSTbinnedPBEs', ...
    'PREbinnedPBEs', ...      
    'RUNbinnedPBEs', ...      
    'POSTbinnedPBEs_p', ...
    'PREbinnedPBEs_p', ...      
    'RUNbinnedPBEs_p', ... 
    'BDseqscore', ...
    'begPosition', ...        
    'behavior', ...           
    'coveredLen', ...         
    'coveredLen_null', ...   
    'endPosition', ...        
    'fileinfo', ...           
    'jumpDist', ...           
    'onlyLineElements', ...   
    'posteriorProbMatrix', ...
    'replayScore', ...        
    'replayScore_null', ...   
    'runTemplate', ...     
    'secondaryPBEs', ...      
    'sigMatrix', ...          
    'spatialInfo', ...          
    'spatialTunings', ...    
    'weightedCorr', ...       
    'weightedCorr_null')
end
    