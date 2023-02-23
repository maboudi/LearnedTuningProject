function linearPos = linearizePosition(fileinfo, behavior, mazeShape)


time = fileinfo.xyt(:, 3);
xpos = fileinfo.xyt(:, 1);
ypos = fileinfo.xyt(:, 2);

linearPos = nan(length(fileinfo.xyt), 1);

runIdx = time > behavior.time(2,1) & time < behavior.time(2,2); % only the positions on the maze

xpos = xpos(runIdx);
ypos = ypos(runIdx);

nPosBins = length(xpos);

posBinWidth = 0.1; % in cm


runLinearPos = nan(nPosBins, 1);

nonNanBins = find(~isnan(xpos) & ~isnan(ypos));

switch mazeShape
    
    
    case 'linear'
        
%         Ypeaks = localexterma(ypos);
        Ypeaks = [min(ypos) max(ypos)];
        
        trackX = min(xpos): posBinWidth: max(xpos);
        trackX(end) =  max(xpos);
        trackY = mean(Ypeaks) * ones(size(trackX));
        
        % calculate the projections
        
        for ii = 1: length(nonNanBins)
            
            bin = nonNanBins(ii);
            
            projectedPnts = projectPnt([xpos(bin) ypos(bin)], [min(xpos) mean(Ypeaks)], [max(xpos) mean(Ypeaks)], 1);

            runLinearPos(bin) = find(histc(projectedPnts(1), trackX))*posBinWidth;

        end
        

            
    case 'L-shape' % we can add an imaginary arc track where the platform is to avoid omitting the position data on the platform
        
        Xpeaks = localexterma(xpos);
        Ypeaks = localexterma(ypos);
        
        
        trackX_1 = min(xpos): posBinWidth: max(Xpeaks);
        trackX_1(end) = max(Xpeaks);
        
        trackY_1 = min(Ypeaks) * ones(size(trackX_1));
        
        
        
        trackY_2 = min(Ypeaks): posBinWidth: max(ypos);
        trackY_2(end) = max(ypos);
        
        trackX_2 = max(Xpeaks) * ones(size(trackY_2));
        
        trackX = [trackX_1 trackX_2];
        trackY = [trackY_1 trackY_2];
        

        % calcualte the projections

        for ii = 1: length(nonNanBins)
            
             bin = nonNanBins(ii);
            
            [projection1, dist1] = projectPnt([xpos(bin) ypos(bin)], [min(xpos) min(Ypeaks)], [max(Xpeaks) min(Ypeaks)], 1);
            
            [projection2, dist2] = projectPnt([xpos(bin) ypos(bin)], [max(Xpeaks) min(Ypeaks)], [max(Xpeaks) max(ypos)], 2);
            

            [~, idx] = min([dist1 dist2]);
            
            switch idx
                case 1
                    projectedPnts = projection1;
                    projBin = find(histc(projectedPnts(1), trackX_1));
                    
                    runLinearPos(bin) = projBin*posBinWidth;
                    
                case 2
                    projectedPnts = projection2;
                    
                    projBin = find(histc(projectedPnts(2), trackY_2));

                    runLinearPos(bin) = (length(trackX_1)- 1 + projBin) *posBinWidth;

            end
                    
        end

              
        
    case 'U-shape'
        
        Xpeaks = localexterma(xpos);
        Ypeaks = localexterma(ypos);

        
        trackX_1 = min(xpos): posBinWidth: max(Xpeaks);
        trackX_1(end) = max(Xpeaks);
        
        trackY_1 = min(Ypeaks) * ones(size(trackX_1));
        
        trackY_2 = min(Ypeaks): posBinWidth: max(Ypeaks);
        trackY_2(end) = max(Ypeaks);
        
        trackX_2 = max(Xpeaks) * ones(size(trackY_2));
        
        trackX_3 = trackX_1(end:-1:1);
        trackY_3 = max(Ypeaks) * ones(size(trackX_3));
        
        
        trackX = [trackX_1 trackX_2 trackX_3];
        trackY = [trackY_1 trackY_2 trackY_3];
        
        
        
        % calcualte the projections

        for ii = 1: length(nonNanBins)
            
             bin = nonNanBins(ii);
            
            [projection1, dist1] = projectPnt([xpos(bin) ypos(bin)], [min(xpos) min(Ypeaks)], [max(Xpeaks) min(Ypeaks)], 1);
            
            [projection2, dist2] = projectPnt([xpos(bin) ypos(bin)], [max(Xpeaks) min(Ypeaks)], [max(Xpeaks) max(Ypeaks)], 2);
            
            [projection3, dist3] = projectPnt([xpos(bin) ypos(bin)], [min(xpos) max(Ypeaks)], [max(Xpeaks) max(Ypeaks)], 1);
    
            
            [~, idx] = min([dist1 dist2 dist3]);
            
            switch idx
                case 1
                    projectedPnts = projection1;
                    
                    projBin = find(histc(projectedPnts(1), trackX_1));
                    runLinearPos(bin) = projBin*posBinWidth;
                    
                    
                case 2
                    projectedPnts = projection2;
                    
                    projBin = find(histc(projectedPnts(2), trackY_2));
                    runLinearPos(bin) = (length(trackX_1)-1 + projBin)*posBinWidth;
                    
                case 3
                    projectedPnts = projection3;
                    
                    projBin = find(histc(projectedPnts(1), trackX_1));

                    runLinearPos(bin) = (length(trackX_1)-1 + length(trackY_2)-1 + length(trackX_1)-1-projBin)*posBinWidth;
  
        
            end
                

        end
        
    
    case 'circular'
        
        cx = min(xpos) + range(xpos)/2;
        cy = min(ypos) + range(ypos)/2;
        
        radius = nanmean(sqrt((xpos - cx).^2 + (ypos - cy).^2));
        
        xgrid = transpose(linspace(cx-radius, cx+radius, 100));
        
        trackX = [xgrid(end:-1:1); xgrid];
        trackY = [cy+sqrt(radius^2 - (xgrid - cx).^2); cy-sqrt(radius^2 - (xgrid - cx).^2)];
        
        % calcualte the projections

        for ii = 1: length(nonNanBins)
            
             bin = nonNanBins(ii);
            
            distFromCntr = sqrt((xpos(bin) - cx)^2 + (ypos(bin) - cy)^2);
            
            if distFromCntr >= radius
                deltaX = -(distFromCntr - radius)/distFromCntr * (xpos(bin) - cx); 
                deltaY = -(distFromCntr - radius)/distFromCntr * (ypos(bin) - cy);
            else %distFromCntr < radius
                deltaX = (radius - distFromCntr)/radius * (xpos(bin) - cx);
                deltaY = (radius - distFromCntr)/radius * (ypos(bin) - cy);
            end
            
            projectedPnts = [xpos(bin)+deltaX ypos(bin)+deltaY];
            
            px = projectedPnts(1);
            py = projectedPnts(2);
            
            if px >= cx && py >= cy
               theta = atan(abs(py-cy)/abs(px-cx));
            elseif px < cx && py >= cy
               theta = (pi/2-atan(abs(py-cy)/abs(px-cx)))+pi/2;
            elseif px < cx && py < cy
               theta = atan(abs(py-cy)/abs(px-cx))+pi;
            elseif px >= cx && py < cy
               theta = (pi/2-atan(abs(py-cy)/abs(px-cx))) + 3*pi/2;
            end
   
            runLinearPos(bin) = radius* (2*pi-theta); % if the animal is moving clockwise, use (2*pi-theta) and if ccw use theta

        end
        

end

figure;
plot(xpos, ypos, '.', 'markersize', 2, 'color', [.7 .7 .7])

hold on

plot(trackX, trackY, '-r', 'linewidth', 2)

linearPos(runIdx) = runLinearPos;


end

function exterma = localexterma(x)


binWidth = 0.5;
bins = min(x):binWidth:max(x);

middleBin = floor(length(bins)/2);

counts = hist(x, bins);

win = gausswindow(20,40);

counts = conv(counts, win, 'same');
counts = transpose(counts);

counts = [counts (1:length(bins))'];

firstHalf = counts(counts(:, 2) < middleBin, :);
[~, firstMaximum] = max(firstHalf(:, 1));

exterma = bins(firstHalf(firstMaximum, 2));

secondHalf = counts(counts(:, 2) >= middleBin, :);
[~, secondMaximum] = max(secondHalf(:, 1));

exterma = [exterma bins(secondHalf(secondMaximum, 2))];


end

function [Q, dist] = projectPnt(pnt, V1, V2, dim)

v = (V2- V1)/norm(V2-V1);
Q = dot(pnt-V1, v)*v + V1;

dist = norm(pnt-Q);

if dim == 1 % x
    condition = (Q(1) > max(V1(1), V2(1)) || Q(1) < min(V1(1), V2(1)));
else % dim = 2
    condition = (Q(2) > max(V1(2), V2(2)) || Q(2) < min(V1(2), V2(2)));
end

if condition
    
    dist1 = norm(V1-Q);
    dist2 = norm(V2-Q);
    
    if dist1 < dist2
        Q = V1;
    else
        Q = V2;
    end
end


end


