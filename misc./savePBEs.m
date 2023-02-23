function savePBEs(PBEs_primary, PBEs_secondary, fileinfo, folderName) 


primaryPBEs = struct('interval', [], 'peakTime', [], 'peakMUA', []);
for ii = 1:length(PBEs_primary)
   primaryPBEs(ii).interval = PBEs_primary(ii, 1:2);
   primaryPBEs(ii).peakTime = PBEs_primary(ii, 3);
   primaryPBEs(ii).peakMUA  = PBEs_primary(ii, 4);
end


secondaryPBEs = struct('interval', [], 'peakTime', [], 'peakMUA', [], 'containsRipple', [], 'nFiringUnits', []);
for ii = 1:length(PBEs_secondary)
   secondaryPBEs(ii).interval = PBEs_secondary(ii, 1:2);
   secondaryPBEs(ii).peakTime = PBEs_secondary(ii, 3);
   secondaryPBEs(ii).peakMUA  = PBEs_secondary(ii, 4);
   
%    secondaryPBEs(ii).containsRipple = PBErippleIdx(ii);
%    secondaryPBEs(ii).nFiringUnits   = nFiringUnits(ii);
end



save(fullfile(folderName,  [fileinfo.name  '-' 'PBEs2.mat']), 'primaryPBEs', 'secondaryPBEs') 

% Filename = fullfile(folderName, [fileinfo.name '.pr3.evt']);
% MakeEvtFile(PBEs_primary(:, 1:3), Filename, {'beg', 'end', 'peak'}, fileinfo.Fs, 1)

Filename = fullfile(folderName, [fileinfo.name '.pb3.evt']);
MakeEvtFile(PBEs_secondary(:, 1:3), Filename, {'beg', 'end', 'peak'}, fileinfo.Fs, 1)


end