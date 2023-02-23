function [Fs, lfpSampleRate, nCh] = getParameters(file)


fileID = fopen(file);

flag1 = 1;
flag2 = 1;
flag3 = 1;

while 1
    line = fgetl(fileID);
    
    if isempty(line)
        break
    else
    
        out1 = sscanf(line, ['\n<samplingRate>' '%d']);
        out2 = sscanf(line, ['\n<nChannels>' '%d']);
        out3 = sscanf(line, ['\n<lfpSamplingRate>' '%d']);

        
        if ~isempty(out1) && flag1 ~= 0
            Fs = out1;
            flag1 = 0;
        end
        
        if ~isempty(out2) && flag2 ~= 0
            nCh = out2;
            flag2 = 0;
        end
        
        if ~isempty(out3) && flag3 ~= 0
            lfpSampleRate = out3;
            flag3 = 0;
        end
        
    end
end

fclose(fileID);

end
    
    
    
    


