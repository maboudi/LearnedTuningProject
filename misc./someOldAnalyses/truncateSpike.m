function truncatedSpike = truncateSpike(spike, whichPart, breakpoint, runBegin, runEnd, Fs, varargin)


if isempty(varargin)
    pix2cm = 1;
else
    pix2cm = varargin{1};
end

if isempty(runBegin)
    runBegin = spike.t(1);
end

if isempty(runEnd)
    runEnd = spike.t(end);
end


%%% determine the track shortening time

%%% The first laps in both long and short track datasets are assigned as
%%% lap number 1

% firstLapSpikes = spike.t(find(spike.lap == 1));
% 
% [~, rmax] = max([0; diff(firstLapSpikes)]);
% 
% breakpoint = firstLapSpikes(rmax); %%% the time of first spike in the first lap after shortening the track

if ~isempty(whichPart)
    figure('Visible','off');

    subplot(2,1,1)
    plot(spike.t(find(spike.t >= runBegin & spike.t < breakpoint))/Fs, spike.x(find(spike.t >= runBegin & spike.t < breakpoint)) * pix2cm, '.', 'color', 'k')
    ylabel('track position', 'fontsize', 20)
    title('Before Track Shortening', 'fontsize', 20)
    set(gca, 'fontsize',16)

    subplot(2,1,2)
    plot(spike.t(find(spike.t >= breakpoint & spike.t <= runEnd))/Fs, spike.x(find(spike.t >= breakpoint & spike.t <= runEnd)) * pix2cm, '.', 'color', 'k')
    xlabel('time', 'fontsize', 20)
    ylabel('track position', 'fontsize', 20)
    title('After Track Shortening', 'fontsize', 20)
    set(gca, 'fontsize',16)
end



dataFields = fieldnames(spike);
noFields = length(dataFields);

truncatedSpike = spike;
for i = 1 : noFields
    
    originalData = getfield(spike, dataFields{i});
    
    if length(originalData) == length(spike.t)

        if isempty(whichPart)
           truncatedData = originalData(find(spike.t >= runBegin & spike.t <= runEnd));
        elseif whichPart == 1
           truncatedData = originalData(find(spike.t >= runBegin & spike.t < breakpoint));
        elseif whichPart == 2
           truncatedData = originalData(find(spike.t >= breakpoint & spike.t <= runEnd));
        end
        
        truncatedSpike = setfield(truncatedSpike, dataFields{i}, truncatedData);
    end 
end


end

