function  shuffle_data = shuffle( data, method, shiftrange )

shuffle_data = zeros(size(data));
switch method
    case 1  %%% cell ID shuffle
        for bin = 1 : size(data,2)
            shuffleID = randperm(size(data,1));
            shuffle_data(:,bin) = data(shuffleID , bin);
        end
    case 2  %%% bin shuffle_separately for each cell
        for cell = 1 : size(data, 1)
             shuffle_bin = randperm(size(data,2));
             shuffle_data(cell,:) = data(cell, shuffle_bin);
        end
         
    case 3  %%% bin shuffle (time swapping)
        shuffle_bin = randperm(size(data,2));
        shuffle_data = data(:, shuffle_bin);
    
    case 4 %%% column cycle shuffling
        shifts = - shiftrange : shiftrange;
        shuffle_data = column_cycle_shuffle(data, shifts);
        
    case 5
        shifts = - shiftrange : shiftrange;
        data_trav = data';
        
        shuffle_data_trav = column_cycle_shuffle(data_trav, shifts);
        shuffle_data = shuffle_data_trav';
end

end

