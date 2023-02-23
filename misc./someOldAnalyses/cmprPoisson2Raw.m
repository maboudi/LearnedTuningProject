function cmprPoisson2Raw(poissonEventsBinnedfiring, eventsBinnedfiring, longORshort, fileinfo, binDur)

%%% comapare the poisson data with the raw data (meanFirings and the
%%% interspike interval (ISI) histograms) 

noUnits = size(eventsBinnedfiring{1,1}, 1); 

temp = cell2mat(eventsBinnedfiring(:,2)');
meanFiring = mean(temp, 2) / binDur;



temp = cell2mat(poissonEventsBinnedfiring(:,2)'); %% mean firings
poissonMeanFiring = mean(temp, 2) / binDur;


%%% ISI histograms

concatData = cell2mat(eventsBinnedfiring(:,1)'); %% concatenate all the events

concatPoissonData = cell2mat(poissonEventsBinnedfiring(:,1)');

UnitsMeanISI = zeros(1, noUnits);
UnitsSkewns = zeros(1, noUnits);

UnitsMeanISI_p = zeros(1, noUnits); %% 'p' at the end of variable's name stands for Poisoon
UnitsSkewns_p = zeros(1, noUnits);

for unit = 1 : noUnits
    
    %%% raw data
    
    unitISI = diff(find(concatData(unit, :))); %% interspike interval in milisecond
    UnitsMeanISI(unit) = mean(unitISI); 
    UnitsSkewns(unit) = (sum((unitISI - mean(unitISI)).^3) ./ length(unitISI)) ./ (var(unitISI, 1).^1.5); %% skewness of distribution
    
    rawISI{unit} = unitISI;
    
    %%% posiion simulated data
    
    unitISI_p = diff(find(concatPoissonData(unit, :))); 
    UnitsMeanISI_p(unit) = mean(unitISI_p); 
    UnitsSkewns_p(unit) = (sum((unitISI_p - mean(unitISI_p)).^3) ./ length(unitISI_p)) ./ (var(unitISI_p, 1).^1.5);
    
    poissonISI{unit} = unitISI_p;

end


%%% plotting the mean firing, mean ISI, and skewness scatter plots with two
%%% groups of raw data and Poisson simulated data

figure('Visible','off');

subplot(1,3,1)

plot(meanFiring, poissonMeanFiring, '.', 'color','k')
hold on
plot(meanFiring, meanFiring, 'color', [.7 .7 .7])
xlabel('raw data mean firing(Hz)', 'fontsize', 20)
ylabel('Poisson simulated mean firing(Hz)', 'fontsize', 20 )
set(gca,'fontsize', 16)


subplot(1,3,2)

plot(UnitsMeanISI, UnitsMeanISI_p, '.', 'color','k')
hold on
plot(UnitsMeanISI, UnitsMeanISI, 'color', [.7 .7 .7])
xlabel('raw data mean ISI(ms)', 'fontsize', 20)
ylabel('Poisson simulated mean ISI(ms)', 'fontsize', 20 )
set(gca,'fontsize', 16)


subplot(1,3,3)

plot(UnitsSkewns, UnitsSkewns_p, '.', 'color','k')
hold on
plot(UnitsSkewns, UnitsSkewns, 'color', [.7 .7 .7])
xlabel('raw data skewness(ms)', 'fontsize', 20)
ylabel('Poisson simulated skewness(ms)', 'fontsize', 20 )
set(gca,'fontsize', 16)


currDir = pwd;
FileBase = [currDir '/' fileinfo.name '/data part' num2str(longORshort)  '/Poisson Stimulated Data'];
mkdir(FileBase)

saveas(gcf, [FileBase, '/PoissonVSRaw.fig'])


end
