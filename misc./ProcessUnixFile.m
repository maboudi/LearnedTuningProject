function spikeInf = ProcessUnixFile(InputSession, ElecNum)

MuElcClus = [];
for shank = 1 : length(ElecNum) 
    Clus = load(strcat(InputSession, '.clu.', num2str(ElecNum(shank))));
    Spiketimes = load(strcat(InputSession, '.res.', num2str(ElecNum(shank))));


    %%% calculation of each spike position, which bin each spike resides
    % end1Crd = [400  640];
    % end2Crd = [10 345];
    % 
    % if end2Crd(1) > end1Crd(1)
    %     ax = 1;
    % else 
    %     ax = -1;
    % end
    % 
    % if end2Crd(2) > end2Crd(2)
    %     ay = 1;
    % else
    %     ay = -1;
    % end
    %     
    %     
    % 
    % binsize = 10; %% length of each position bin
    % line_length = ((end2Crd(2) - end1Crd(2))^2 + (end2Crd(1) - end1Crd(1))^2)^0.5;
    % a = (end2Crd(2) - end1Crd(2)) / (end2Crd(1) - end1Crd(1));
    % deltaX = (binsize^2 / (1 + a^2))^0.5;
    % deltaY = a * deltaX;
    % 
    % NoPosBins = floor(line_length/binsize);
    % 
    % for ibin = 1 : (NoPosBins - 1)
    %     pntCrd(ibin , :) = end1Crd + ibin*[ax*deltaX ay*deltaY];
    %     
    %     if ibin == 1
    %         binCord(ibin,:) = [end1Crd pntCrd(ibin , :)];
    %     else
    %         binCord(ibin,:) = [pntCrd(ibin , :) pntCrd(ibin-1 , :)];
    %     end
    % end
    % 
    % binCord(NoPosBins, :) = [pntCrd(end , :) end2Crd];
    % 
    % for ibin = 1 : NoPosBins
    %     
    %     rr(ibin) = ((binCord(ibin, 3) - binCord(ibin, 1))^2 + (binCord(ibin, 4) - binCord(ibin, 2))^2)^0.5;
    % end
    % 
    % 
    % 
    % plot(pntCrd(:,1), pntCrd(:,2), '*')

    %% processing spikes

    spktimes = round(Spiketimes ./ 20); %% spikes of 1 ms  resolution
    
    cur_spktrain_length = spktimes(end);
    cum_NoClusteres = size(MuElcClus , 1);
    cur_NoClusters = Clus(1) - 1;

    prv_spk_lenght = 0; %% just an initiation
    if shank > 1
        
        prv_spk_lenght = size(spkMat, 2);
        if cur_spktrain_length > size(spkMat, 2)
            spkMat = [spkMat zeros(size(spkMat, 1), (cur_spktrain_length - prv_spk_lenght))];
        end  
    end
    spkMat(cum_NoClusteres + 1 : cum_NoClusteres + cur_NoClusters, :) = zeros(cur_NoClusters, max(cur_spktrain_length, prv_spk_lenght));
    
    
    for spk = 1 : length(spktimes)
        spk_clu = Clus(spk + 1);

        if spk_clu > 1
            spkMat((spk_clu - 1) + cum_NoClusteres , spktimes(spk)) = 1;
        end
    end

    MuElcClus = [MuElcClus; [1:cur_NoClusters shank * ones(cur_NoClusters, 1)]];
end
 

    SpkRaster = spkMat(:, 1:100000);

    for unit = 1 : size(SpkRaster, 1)

        timeind = find(SpkRaster(unit, :));
        if ~isempty(timeind)
            plot(timeind, unit, '.', 'color', 'k', 'markersize', 8)
            hold on
        end
        ylim([0 (size(SpkRaster, 1) + 1)])
    end
%     

spikeInf = [];






