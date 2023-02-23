parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/sessions_calculated_PBEinfo3';


storageDir = '/data/GreatLakes_datasets_temp';

files = dir(parentDir);
isdir = [files.isdir];
subfolders = files(isdir);

sessionNames = subfolders(3:end);

maxN = 1000;

for sessionNumber = 9%1:numel(sessionNames)
    
    sessionName = sessionNames(sessionNumber).name;
    basePath = fullfile(parentDir, sessionName);


    load(fullfile(basePath, 'PopulationBurstEvents', [sessionName '.PBEInfo.mat']), 'PBEInfo_Bayesian')
    
    
    
    PBEInfo_Bayesian_all = PBEInfo_Bayesian;
    
    nT_PBEs = numel(PBEInfo_Bayesian_all);
    
    % if number of PBEs was higher than maxN, we parse them to smaller
    % datasets

   nSets = ceil(nT_PBEs/maxN);
       
   for ii = 1:nSets

       if ii < nSets
           PBEInfo_Bayesian = PBEInfo_Bayesian_all((ii-1)*maxN+1 : ii*maxN);
       else
           PBEInfo_Bayesian = PBEInfo_Bayesian_all((ii-1)*maxN+1 : nT_PBEs);
       end


       % save the subset of PBEs as a dataset
       sessionNames_subset = [sessionName '_' num2str(ii)];
       session_subset_folder = fullfile(storageDir, sessionNames_subset);
       if ~exist(session_subset_folder, 'dir')
        mkdir(session_subset_folder)
       end

       folderName = fullfile(session_subset_folder,  'PopulationBurstEvents');

       if ~exist(folderName, 'dir')
          mkdir(folderName)
       end
       save(fullfile(folderName, [sessionNames_subset '.PBEInfo.mat']), 'PBEInfo_Bayesian', '-v7.3')


       % copy spikes -- I couldn't figure out how to work with copyfile
       % function

       load(fullfile(basePath, 'spikes', [sessionName '.spikes.mat']));

       folderName = fullfile(session_subset_folder, 'spikes');

       if ~exist(folderName, 'dir')
            mkdir(folderName)
       end

       save(fullfile(folderName, [sessionNames_subset '.spikes.mat']), 'spikes', 'spikes_pyr', 'fileInfo', '-v7.3')

   end
       
    
end
    
    
    
    