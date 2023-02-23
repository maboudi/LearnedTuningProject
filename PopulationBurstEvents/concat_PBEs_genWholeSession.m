parentDir = '/data/GreatLakes_datasets_temp/RatU_Day2/';


storageDir = '/data/GreatLakes_datasets_temp/RatU_Day2';

files = dir(parentDir);
isdir = [files.isdir];
subfolders = files(isdir);

folderNames = subfolders(3:end);

nSubSessions = numel({folderNames.name});
sessionNames = cell(nSubSessions, 1);
for ii = 1:nSubSessions
    sessionNames{ii} = folderNames(ii).name(1:find(folderNames(ii).name == '_', 1, 'last')-1);
end


uniqSessionNames = unique(sessionNames);


for ii = 1:numel(uniqSessionNames)

    sessionDir = fullfile(storageDir, uniqSessionNames{ii});
    mkdir2(sessionDir)
    
    curr_subSession_nums = find(strcmp(sessionNames, uniqSessionNames{ii}));
    PBEInfo_replayScores_all = [];
    PBEInfo_null_all         = struct('radonIntegral', [], 'coveredLen', [], 'weightedCorr', [], 'jumpDist', []);
    
    
    for jj = 1:numel(curr_subSession_nums)
        
        
        subSessionName = folderNames(curr_subSession_nums(jj)).name;
        
        subSessionDir = fullfile(parentDir, subSessionName);
        
        
        % BD replay scores
        load(fullfile(subSessionDir, 'BayesianDecoding', [subSessionName '.PBEInfo_replayScores.mat']), 'PBEInfo_replayScores', 'PBEInfo_null')
        
        PBEInfo_replayScores_all = [PBEInfo_replayScores_all PBEInfo_replayScores];
        
        PBEInfo_null_all.radonIntegral = cat(1, PBEInfo_null_all.radonIntegral, PBEInfo_null.radonIntegral);
        PBEInfo_null_all.coveredLen    = cat(1, PBEInfo_null_all.coveredLen, PBEInfo_null.coveredLen);
        PBEInfo_null_all.weightedCorr  = cat(1, PBEInfo_null_all.weightedCorr, PBEInfo_null.weightedCorr);
        
        PBEInfo_null_all.jumpDist      = [PBEInfo_null_all.jumpDist PBEInfo_null.jumpDist];
        
%         if jj == 1
%            load(fullfile(subSessionDir, 'spikes', [subSessionName '.spikes.mat']), 'fileInfo', 'spikes', 'spikes_pyr') 
%         end
       
    end
  
    PBEInfo_null         = PBEInfo_null_all;
    PBEInfo_replayScores = PBEInfo_replayScores_all;
    
    temp = fullfile(sessionDir, 'BayesianDecoding');
    mkdir2(temp)
    save(fullfile(temp, [uniqSessionNames{ii} '.PBEInfo_replayScores.mat']), 'PBEInfo_replayScores', 'PBEInfo_null', '-v7.3')
    

%     temp = fullfile(sessionDir, 'spikes');
%     mkdir2(temp)
%     save(fullfile(temp, [uniqSessionNames{ii} '.spikes.mat']), 'fileInfo', 'spikes', 'spikes_pyr')

end


function mkdir2(dirName)

if ~exist(dirName, 'dir'); mkdir(dirName); end

end
        
