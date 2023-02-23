function   speedprofiles(fileinfo, spike, Fs)


%%% Loading position data

xyt = fileinfo.xyt;
tbegin = fileinfo.tbegin;
tend = fileinfo.tend;

if ~isfield(fileinfo, 'pix2cm')
    pix2cm = 1;
else
    pix2cm = fileinfo.pix2cm;
end


withinRange = find(xyt(:,3) >= tbegin & xyt(:,3) <= tend);
xyt = xyt(withinRange, :);

xpos = xyt(:,1)* pix2cm;
ypos = xyt(:,2)* pix2cm;

tt = (xyt(:,3) - tbegin)/10^6; %%% position sampling times in second

%%% partitioning the the whole session to travel and rest periods (by
%%% travel I mean each time that animal run on the track from one end to
%%% another)

runWin = ones(15,1)/15;

diffx = [0; abs(diff(filtfilt(runWin,1,xpos)))];
diffy = [0; abs(diff(filtfilt(runWin,1,ypos)))];

difftt = [1; diff(tt)];

speed = sqrt(diffx.^2 + diffy.^2)./difftt; %% speed before smoothing


%%%%% smoothing the speed (need do more smoothing of the speed in order to
%%%%% categorize the whole periods to the distinct periods (travel and rest
%%%%% periods)

sigma = 15; %%% the frequency of sampling is about 30 Hz, so the sigma here is about 500 ms
halfwidth = 3 * sigma;

smoothwin = gausswindow(sigma, halfwidth);
speed_sm = conv(speed, smoothwin, 'same');


runposind = find(speed_sm > 15); %%% finding periods with speed (smoothed) of higher than the run threshold


%%% plotting the position and speed profiles together to have a sense of
%%% their covariation

figure%('Visible','off')

subplot(3,1,1)

plot(tt,xpos, 'k')
hold on
plot(tt(runposind), xpos(runposind), '.', 'color','b') %%% highlighting times with speed exceeding the threshold (10 cm/s)
ylabel('Pos (cm)', 'fontsize', 20)
title('Speed > 15 cm/s (blue)', 'fontsize',20)
set(gca,'fontsize',16);

subplot(3,1,2)
spikeInd = find(spike.lap > 0);
plot(tt,xpos, 'k')
hold on
plot(spike.t(spikeInd)/Fs, spike.x(spikeInd)* pix2cm,'.','color','r')
ylabel('Pos (cm)', 'fontsize', 20)
title('Laps (red)', 'fontsize',20)
set(gca,'fontsize',16);

subplot(3,1,3)

plot(tt, speed_sm, 'k')
hold on
plot(tt, speed, 'color', [.8 .8 .8])
ylabel('Speed (cm/s)', 'fontsize', 20)
xlabel('Time (sec)', 'fontsize', 20)
set(gca,'fontsize',16)



%%% finding the start and end times of the periods
%%% example: runposind = [2  3  4  8  9  10 11 12 18 19 20]
%%%%           diffind = [2  1  1  5  1  1  1  1  6  1  1 ]
%%%%                     (b)      (b)            (b)
%%%%% (b)s signifies the break points

diffind = [runposind(1); diff(runposind)];
breaks = find(diffind > 1);

runBegins = runposind(breaks);
runEnds = runposind(breaks(2:end)-1); %% one sample before a break, signifying ...
                                      %%the last sample of the previous high speed period
runEnds = [runEnds; runposind(end)];

runDur = runEnds - runBegins;

figure('Visible','off'); 
hist(runDur* mean(difftt(2:end)),50)
xlabel('Duration(sec)','fontsize',20)
ylabel('Number', 'fontsize', 20)
set(gca,'fontsize',16)

end


