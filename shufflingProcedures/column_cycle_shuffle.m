function  newData  = column_cycle_shuffle(data, randShift)


if iscell(data)
    data = data{:};
end

noRows = size(data , 1);
noColumns = size(data , 2);

newData = zeros(size(data));

for ii = 1: noColumns
   
    newData(:, ii) = circshift(data(:, ii), randShift(ii));
    
end

% 
% %% old way
% newData = zeros(size(data));
% 
% for col = 1 : noColumns
% 
%     shift = randShift(col);
% 
%     if shift < 0 
%        shift = -shift;
% 
%        newData(1 : (noRows - shift), col) = data((shift+1) : noRows, col);
%        newData((noRows - shift + 1) : noRows, col) = data(1 : shift, col);
%        
%     elseif shift > 0
% 
%         newData((shift+1) : noRows, col) = data(1 : (noRows - shift), col);
%         newData(1 : shift, col) = data((noRows - shift + 1) : noRows, col);
%         
%     elseif shift == 0
%         
%         newData(:,col) = data(:,col);
%         
%     end
% end


end

