function position2 = replaceNans(position, time, nanIdx)


position_woNans         = position;
position_woNans(nanIdx) = [];


time_woNans         = time;
time_woNans(nanIdx) = [];


intrpltPos = interp1(time_woNans, position_woNans, time(nanIdx));


position2        = position;
position2(nanIdx) = intrpltPos;



end