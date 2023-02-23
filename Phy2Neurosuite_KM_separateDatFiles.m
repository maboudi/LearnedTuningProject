function Phy2Neurosuite_KM_separateDatFiles(fileBase)

% this function extract the wavforms and computes PCs, if we already have
% the spike times, .xml and .dat files.



% badChannels
% in Utku's dataset: AG_2019_12_23_NSD there was no bad channel

% load([fileBase '.basics.mat'], 'badchans', 'auxchans')
% badChannels = [badchans auxchans]+1; %add one if indices are based on Neuroscope
badChannels = []; 


% parameters
par = LoadXml([fileBase '.shx.xml']);

params.ForceMaxRAMforDat = 15000000000;
params.Nchan             = numel(par.AnatGrps.Channels);
params.ntbuff            = 64;
params.fshigh            = 500;
params.fs                = par.SampleRate;
params.nt0               = round(1.6*params.fs/1000); % number of samples around each spikes (at fs = 30k, it's equal to 48 sample)
params.GPU               = 0;
params.NT                = 32*1028+ params.ntbuff;
params.nTotCh            = par.nChannels;
params.badChannels       = badChannels;


% load spikes
load([fileBase '.spikes.mat'], 'spikes');


for iUnit = 1:numel(spikes)
    spikes(iUnit).shank = spikes(iUnit).id(1); 
end

uniqShanks = unique([spikes.shank]);
uniqShanks = sort(uniqShanks, 'ascend');


for sh = 2:numel(uniqShanks)
    
    
    shankNumber = uniqShanks(sh);
    
    % .dat file
    datFile = sprintf([fileBase '.sh%d.dat'], shankNumber);
    
    
    units = find([spikes.shank] == shankNumber);
    
    nUnits = numel(units); 
    
%     channels   = par.AnatGrps(shankNumber).Channels+1;
    channels   = 1:params.Nchan;
    channels(ismember(channels, params.badChannels)) = [];
    
    
    tspktimes  = cell(nUnits, 1);
    tclus      = cell(nUnits, 1);
    nSpikes    = zeros(nUnits, 1);
    
    for iUnit = 1:nUnits
        
        currUnit = units(iUnit);
        
        tspktimes{iUnit} = spikes(currUnit).time * params.fs;
        nSpikes(iUnit)   = numel(tspktimes{iUnit});
        tclus{iUnit}     = spikes(currUnit).id(2)*ones(nSpikes(iUnit), 1);
        
    end
    
    
    tspktimes = cell2mat(tspktimes);
    [tspktimes, sortIdx]  = sort(tspktimes, 'ascend');
    
    tclus     = cell2mat(tclus);
    tclus     = tclus(sortIdx);

    tclus     = cat(1, numel(unique(tclus)), tclus);
    
    
    % For rat N, we need to exluded 68 to 172 minutes from the dat file or
    % accrodingly add to the spike times
    
    tspktimes2 = tspktimes; 
%     adjustIdx = tspktimes2 > 168*60*params.fs;
%     tspktimes2(adjustIdx) = tspktimes2(adjustIdx) + 4*60*params.fs;
    
    tspkwvs = extractWaveforms(datFile, tspktimes2, channels, params);
%     tspkwvs = extractWaveforms(datFile, tspktimes, channels, params);
    
    
    
    fprintf('\n')
    fprintf('shank %d ..', shankNumber)
        
    fprintf('\nSaving .clu files to disk (cluster indexes)')
    cluname = fullfile([fileBase '.clu.' num2str(shankNumber)]);
    fid=fopen(cluname,'w');
    fprintf(fid,'%.0f\n',tclus);
    fclose(fid);
    clear fid
    fprintf('\n')

    
    fprintf('\nSaving .res files to disk (spike times)')
    resname = fullfile([fileBase '.res.' num2str(shankNumber)]);
    fid=fopen(resname,'w');
    fprintf(fid,'%.0f\n',tspktimes);
    fclose(fid);
    clear fid

    fprintf('\n');

    fprintf('\nSaving .spk files to disk (waveforms)')

    fid=fopen([fileBase,'.spk.',num2str(shankNumber)],'w');
    fwrite(fid, tspkwvs,'int16');
    fclose(fid);

    fprintf('\n'); 
    
    
    
    fprintf('\nComputing PCAs')

    PCAs_global = zeros(3, size(tspkwvs, 1), size(tspkwvs, 3));
    waveforms = tspkwvs;

    waveforms2 = reshape(waveforms,[size(waveforms,1)*size(waveforms,2),size(waveforms,3)]);
    wranges = int64(range(waveforms2,1));
    
    
    
%     wpowers = int64(sum(waveforms2.^2,1)/size(waveforms2,1)/100); % out of memory sometimes if there are too many spikes
    
    temp = floor(numel(tspktimes)/2);
    
    wpowers1 = int64(sum(waveforms2(:, 1:temp).^2,1)/size(waveforms2(:, 1:temp),1)/100);
    wpowers2 = int64(sum(waveforms2(:, temp+1:end).^2,1)/size(waveforms2(:, temp+1:end),1)/100);
    wpowers  = [wpowers1 wpowers2];
    


    % Calculating PCs 
    for k = 1:size(waveforms,1) % for each channel
        PCAs_global(:,k,:) = pca(zscore(permute(waveforms(k,:,:),[2,3,1]),[],2),'NumComponents',3)';
    end
 
     
    fprintf(['Saving .fet files for group ', num2str(shankNumber),'. '])
    PCAs_global2 = reshape(PCAs_global,size(PCAs_global,1)*size(PCAs_global,2),size(PCAs_global,3));
    factor = (2^15)./max(abs(PCAs_global2'));
    PCAs_global2 = int64(PCAs_global2 .* factor');

    fid=fopen([fileBase,'.fet.',num2str(shankNumber)],'w');
    Fet = double([PCAs_global2; wranges; wpowers; tspktimes']);
    nFeatures = size(Fet, 1);
    formatstring = '%d';
    for ii=2:nFeatures
        formatstring = [formatstring,'\t%d'];
    end
    formatstring = [formatstring,'\n'];

    fprintf(fid, '%d\n', nFeatures);
    fprintf(fid,formatstring,Fet);
    fclose(fid);

    fprintf('\n')
    fprintf('\nComplete!')
    
end
    
    
end

function waveforms = extractWaveforms(datFile, spikeTimes, channels, params)



[b1, a1] = butter(3, params.fshigh/params.fs*2, 'high');

d = dir(datFile);
fid = fopen(datFile, 'r');

NTbuff      = params.NT + 4*params.ntbuff;
Nbatch      = ceil(d.bytes/2/params.nTotCh /(params.NT-params.ntbuff));


waveforms = zeros(numel(channels), params.nt0, numel(spikeTimes));


fprintf('Extraction of waveforms begun \n')
for ibatch = 1:1000%Nbatch
    if mod(ibatch,10)==0
        if ibatch~=10
            fprintf(repmat('\b',[1 length([num2str(round(100*(ibatch-10)/Nbatch)), ' percent complete'])]))
        end
        fprintf('%d percent complete', round(100*ibatch/Nbatch));
    end

    offset = max(0, 2*params.nTotCh*((params.NT - params.ntbuff) * (ibatch-1) - 2*params.ntbuff));
    if ibatch==1
        ioffset = 0;
    else
        ioffset = params.ntbuff;
    end
    fseek(fid, offset, 'bof');
    buff = fread(fid, [params.nTotCh NTbuff], '*int16');

    if isempty(buff)
        break;
    end
    nsampcurr = size(buff,2);
    if nsampcurr<NTbuff
        buff(:, nsampcurr+1:NTbuff) = repmat(buff(:,nsampcurr), 1, NTbuff-nsampcurr);
    end
    if params.GPU
        try%control for if gpu is busy
            dataRAW = gpuArray(buff);
        catch
            dataRAW = buff;
        end
    else
        dataRAW = buff;
    end

    dataRAW = dataRAW';
    dataRAW = single(dataRAW);
    dataRAW = dataRAW-median(dataRAW,2);
    datr = filter(b1, a1, dataRAW);
    datr = flipud(datr);
    datr = filter(b1, a1, datr);
    datr = flipud(datr);
    DATA = gather(int16(datr(ioffset + (1:params.NT),:)));
    dat_offset = offset/params.nTotCh/2+ioffset;

    % Saves the waveforms occuring within each batch

    temp = find(ismember(spikeTimes, (params.nt0/2+1:size(DATA,1)-params.nt0/2)+dat_offset));
    temp2 = spikeTimes(temp)-dat_offset;

    startIndicies = temp2-params.nt0/2+1;
    stopIndicies = temp2+params.nt0/2;
    X = cumsum(accumarray(cumsum([1;stopIndicies(:)-startIndicies(:)+1]),[startIndicies(:);0]-[0;stopIndicies(:)]-1)+1);
    X = X(1:end-1);
    
    waveforms(:,:,temp) = reshape(DATA(X, channels)',numel(channels), params.nt0, []);
        
end
        
end