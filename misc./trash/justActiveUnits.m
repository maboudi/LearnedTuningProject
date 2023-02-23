function newData = justActiveUnits(Data, activeUnits)

noEvents = size(Data, 1);
noBinnigMethods = size(Data, 2);

newData = cell(size(Data));

for evt = 1 : noEvents
    for binnig = 1 : noBinnigMethods
        newData{evt,binnig} = Data{evt,binnig}(activeUnits, :);
    end
    
end

end