% function U = unity(A,exclude) 
% A is n x m  matrix, and each column is a signal.
% exclude are pointers for which data will be excluded, and therefore
% should not factor into std or median.
% the function gives a unity signal for each column
% i.e. [signal-mean(signal)] / std(signal)

function [U, stdA] = unity(A,exclude)

ii_ok = ones(size(A,1),1);  % all the indices
if nargin>1;
    fprintf('excluding undesired indices from signal...\n')
    exclude = round(exclude); % make sure these there integers
    for ii = 1:size(exclude,1)
        if size(A,2)>size(A,1)
            error('input should be a column vector! \n')
        end
        ii_ok(exclude(ii,1):exclude(ii,2)) = 0; % zero out unwanted indices
    end
%     A = zeros(length(find(ii_ok)),size(B,2));
%     for jj = 1:size(A,2)
%         A(:,jj) = B(find(ii_ok),jj);
%     end
end
ii_ok = find(ii_ok);
U = zeros(size(A));
U(ii_ok,:) = A(ii_ok,:);  %  U is all zeros, except for the non-excluded elements of A
%     stdA = std(A);
meanA = mean(A(ii_ok,:)); 
stdA = median(abs(A(ii_ok,:))/sqrt(2)/0.6745);
 
U = (U - repmat(meanA,size(A,1),1))./repmat(stdA,size(A,1),1);


