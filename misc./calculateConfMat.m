function [confMat, confMatIndices] = calculateConfMat(Sa, Sb, noEvents, FileBase, dataType)



aOK = find(Sa > 1); %% indices of significant events considering Score a (Sa) 
bOK = find(Sb > 1); %% indices of significant events considereing Score b (Sb)


%%% calculating the confusion matrix

a_b = intersect(aOK, bOK); %%% both the a score and b score are significant

a_nb = aOK(~ismember(aOK, bOK)); %%% events with exclusively a Score significant

na_b = bOK(~ismember(bOK, aOK)); %%% events with exclusively b score significant 

na_nb = 1:noEvents;
na_nb(union(aOK, bOK)) = []; %%%% niether of scores are significant


%%%% confusion Matrix

confMat = [length(a_b) length(a_nb); length(na_b) length(na_nb)];

confMatIndices = zeros(noEvents, 1);

confMatIndices(a_b) = 1;
confMatIndices(a_nb) = 2;
confMatIndices(na_b) = 3;
confMatIndices(na_nb) = 4;

FileBase2 = [FileBase '/' dataType];
mkdir(FileBase2)

save([FileBase2 '/confmat.mat'], 'confMat', 'confMatIndices')

end