% reads multi-channel recording file to a matrix
% function [eeg] = function readmulti(fname,numchannel,chselect)
% last argument is optional (if omitted, it will read all the 
% channels

function [eeg] = readmulti(fname,numchannel,chselect,subtract_channel)

if nargin == 2
  datafile = fopen(fname,'r');
  eeg = fread(datafile,[numchannel,inf],'int16');
  fclose(datafile);
  eeg = eeg';
  return
end

if nargin == 3

  % the real buffer will be buffersize * numch * 2 bytes
  % (short = 2bytes)
  
  buffersize = 4096;
  
  % get file size, and calculate the number of samples per channel
  fileinfo = dir(fname);
  numel = ceil(fileinfo(1).bytes / 2 / numchannel);
  
  datafile = fopen(fname,'r');
  
  mmm = sprintf('%d elements',numel);
%  disp(mmm);  
  
  eeg=zeros(length(chselect),numel);
  numel=0;
  numelm=0;
  while ~feof(datafile),
    [data,count] = fread(datafile,[numchannel,buffersize],'int16');
    if count~=0
        numelm = count/numchannel;
        eeg(:,numel+1:numel+numelm) = data(chselect,:);
        numel = numel+numelm;
    end
end
fclose(datafile);
end

if nargin == 4

  % the real buffer will be buffersize * numch * 2 bytes
  % (short = 2bytes)
  
  buffersize = 4096;
  
  % get file size, and calculate the number of samples per channel
  fileinfo = dir(fname);
  numel = ceil(fileinfo.bytes / 2 / numchannel);
  
  datafile = fopen(fname,'r');
  
  mmm = sprintf('%d elements',numel);
%  disp(mmm);  
  
  eeg=zeros(length(chselect),numel);
  numel=0;
  numelm=0;
  while ~feof(datafile),
    [data,count] = fread(datafile,[numchannel,buffersize],'int16');
    if count~=0
        numelm = count/numchannel;
        eeg(:,numel+1:numel+numelm) = data(chselect,:)-repmat(data(subtract_channel,:),length(chselect),1);
        numel = numel+numelm;
    end
  end
fclose(datafile);
end


eeg = eeg';
