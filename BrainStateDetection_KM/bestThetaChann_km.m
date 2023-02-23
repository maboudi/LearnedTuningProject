function thetaChannels = bestThetaChann_km(fileName,sampleRate,nCh,numShanks)

%     fileName=[basics.FileName,'.eeg'];
%     sampleRate = basics.lfpSampleRate;
%     nCh=basics.nChannels;
%     CA1Shanks = basics.Ch.CA1Shanks;
%     Channels=basics.SpkGrps

    fNum= 2^floor(log(600*sampleRate)/log(2)+1);

    fh = fopen(fileName);
    data = [fread(fh,[nCh,fNum],'int16')]';
    fclose(fh);

%     temp=fft(data);

    forder = 256;  % filter order has to be even. The longer the more selective, but the operation
    % will be linearly slower to the filter order. 100-125 for 1.25Hz (500 for 5 KH
    forder = ceil(forder/2)*2; % to make sure filter order is even
    hTheta = 10; % bandpass filter range
    lTheta = 6; %

    hDelta = 5; % bandpass filter range
    lDelta = 1; % 

    fTheta = fir1(forder,[lTheta/sampleRate*2,hTheta/sampleRate*2]); % calculate convolution func
    dTheta=Filter0(fTheta,data);


    avgTheta = mean(dTheta.^2);
    
    nchpershank = nCh/numShanks;
    thetaChannels=[];
    for shank = 1:numShanks
        ch = (shank-1)*nchpershank+1 : shank*nchpershank;
%         [dummy,pos]=max(cv(ch));
        [dummy,pos]=max(avgTheta(ch));
        thetaChannels=[thetaChannels,ch(pos)];
    end
    
    
%     fDelta = fir1(forder,[lDelta/sampleRate*2,hDelta/sampleRate*2]); % calculate convolution func
%     dDelta=Filter0(fDelta,data(:,thetaChannels));
%     avgDelta = mean(dDelta.^2);
%     ratio = avgTheta(thetaChannels)./avgDelta;
%     [dummy,idx]=sort(ratio,'descend');
    
     [dummy,idx]=sort(avgTheta(thetaChannels),'descend');    
    thetaChannels=thetaChannels(idx);
    
end

    


