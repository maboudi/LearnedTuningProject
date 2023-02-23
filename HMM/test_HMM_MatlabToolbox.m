ns = 1000;
a = mod(floor((1:ns)/50),2) ==1;
y = zeros(size(a));
y(a==0) = rand(nnz(a==0),1)<.3;
y(a==1) = rand(nnz(a==1),1)<.05;

x(a==0) = rand(nnz(a==0),1)<.3;
x(a==1) = rand(nnz(a==1),1)<.05;

SpkRaster = [y;x];


for unit = 1 : size(SpkRaster, 1)

    timeind = find(SpkRaster(unit, :));
    if ~isempty(timeind)
        plot(timeind, unit, '.', 'color', 'k', 'markersize', 8)
        hold on
    end
    ylim([0 (size(SpkRaster, 1) + 1)])
end


% plot(y);axis([0,1000,0,1.3]);

%%
%The guessed probability of the model switching from the low to the
%high-firing rate state:
switchProb = [.9,.1;.1,.9];
%The guessed probability of emitting the symbols (0 = no spike, 1 = spike)
%in each state:
%Note that state 1 corresponds to low firing rate, and state 2 is
%high firing rate
emissProb = [.9,.1;
             .3,.6];
[T,emissions] = hmmtrain(y,switchProb,emissProb,'symbols',[0,1]);
T
emissions



%%
%Probability of each state
[pstates] = hmmdecode(y,T,emissions,'symbols',[0,1]);
subplot(2,1,1);plot(1:length(y),y,1:length(y),pstates(2,:))
 
%Most probable sequence of states
[ostate] = hmmviterbi(y,T,emissions,'symbols',[0,1]);
subplot(2,1,2);plot(1:length(y),y,1:length(y),1.2*ostate-1.3);
axis([0,1000,-.2,1.3]);