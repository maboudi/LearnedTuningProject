% Fet = LoadFet(FileName)
%
% A simple matlab function to load a .fet file

function [Fet, nFeatures] = LoadFet(FileName,Spikes2Load);

if nargin<2;Spikes2Load = inf;end

Fp = fopen(FileName, 'r');

if Fp==-1
    error(['Could not open file ' FileName]);
end

nFeatures = fscanf(Fp, '%d', 1);
Fet = fscanf(Fp, '%f', [nFeatures, Spikes2Load])';

fclose(Fp);