
clear;clc;

for Neuron = 1 : 10
    spike_neuron = [];
    timevec_neuron = [];
    for Timebin = 1 : 10
        ind = rand(1);
        if Timebin == Neuron
            fr = 80;
        elseif Timebin == Neuron + 1 | Timebin == Neuron - 1
            ii = ind > 0.3
            fr = ii*40;
        elseif Timebin == Neuron + 2 | Timebin == Neuron - 2
            ii = ind > 0.5
            fr = ii*20;
        elseif Timebin == Neuron + 3 | Timebin == Neuron - 3
            fr = 0;
        else
            fr = 0;
        end
        
        [spikeMat, tVec] = poissonSpikeGen(fr, 0.25, 1);
        
        spike_neuron = [spike_neuron spikeMat];
        timevec_neuron = [timevec_neuron 0.25*(Timebin-1) + tVec];
    end
    spike(Neuron,:) = spike_neuron;
    timevec(Neuron,:) = timevec_neuron;
end

plotRaster(logical(spike), timevec(1,:)*1000);
xlabel('Time (ms)');
ylabel('Trial Number');

save('spike.mat','spike')