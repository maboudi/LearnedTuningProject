function [Spk, nCh] = readSpk(FileName, nCh, smpPerSpk)

    fh = fopen(FileName);
    Spk = fread(fh, [nCh, Inf],'int16');
    fclose(fh);
    
    nSpikes = size(Spk, 2)/smpPerSpk;

    if floor(nSpikes)~= nSpikes
        
        clear Spk
        warning('The number of channels do not seem correct!')
        
%         for nCh2 = [nCh+1 nCh+2 nCh+3 nCh+4]
%         
%             fh = fopen(FileName);
%             Spk = fread(fh, [nCh2, Inf],'int16');
%             fclose(fh);
%             
%             nSpikes = size(Spk, 2)/smpPerSpk;
%             
%             if ~isfloat(nSpikes)
%                 break
%             end
%         end
        fh = fopen(FileName);
        Spk = fread(fh,'int16');
        fclose(fh);
        
        nSpikeChann = numel(Spk)/smpPerSpk;
        
        for nCh2 = setdiff(nCh-3:nCh+3, nCh)
            nSpikes = nSpikeChann/nCh2;
            if floor(nSpikes)== nSpikes
                nCh = nCh2;
                break
            end
        end
    end
    clear Spk

    fh = fopen(FileName);
    Spk = fread(fh, [nCh, Inf],'int16');
    fclose(fh);
    
    nSpikes = size(Spk, 2)/smpPerSpk;
    
    Spk = reshape(Spk, [nCh, smpPerSpk, nSpikes]);
end

