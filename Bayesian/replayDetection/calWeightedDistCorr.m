function weightedDistCorr = calWeightedDistCorr(postPr)


options.dmax     = 1.6;
options.alpha    = 8;
options.option_y = 'exponential';
options.option_x = 'linear_bonded';
options.delta_t  = 5;
options.Nw = 1;
    

[nPosBins, nTimeBins] = size(postPr);

temp = eye(nPosBins);
win = gausswindow(5, 10);

trans_matrix_super = zeros(size(temp));

for ii = 1:nPosBins
    trans_matrix_super(ii, :) = conv(temp(ii, :), win, 'same');
end


dist_states_matrix = get_states_distance_matrix(trans_matrix_super, options);

dist_time = get_times_distance_matrix(nTimeBins, options); % get time-distance matrix

weightedDistCorr = weightedDistanceCorr(postPr', dist_time, dist_states_matrix, options);


end