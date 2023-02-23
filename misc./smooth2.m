function smoothed_tuning = smooth2(tuning, sigma, halfwidth)

win = gausswindow(sigma, halfwidth);

smoothed_tuning = zeros(size(tuning));
for i = 1 : size(tuning, 1)
    
    smoothed_tuning(i,:) = conv(tuning(i,:), win, 'same');
    
end

smoothed_tuning = smoothed_tuning + 0.001;

end