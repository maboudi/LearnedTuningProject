function [replaydirection, replayOrder] = replaytype(postPrLeft, postPrRight, bestSlope, bestrho)


%%% calculating the sum of probability around the fitting line 

sumLikeLeft = sumalongline(postPrLeft, bestrho, bestSlope); %% RL (leftward) direction
sumLikeRight = sumalongline(postPrRight, bestrho, bestSlope); %% LR (rightward) direction


replaydirection = (sumLikeRight - sumLikeLeft)/ (sumLikeRight + sumLikeLeft);  %% above zero --> replay of LR, less than zero --> replay of RL

replayOrder = replaydirection * sign(bestSlope); %% above zero for forward and less than zero for reverse



end
