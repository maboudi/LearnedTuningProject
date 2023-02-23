function positions = subplot_pos2(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, spacer, spacec)
    
    subheight = (plotheight - topmargin - bottommargin - spacer * (nor-1.0)) / nor;
    subwidth = (plotwidth - leftmargin - rightmargin - spacec * (noc-1.0)) / noc;
    
    
    positions = cell(4, noc);
    for i = 1:nor
       for j = 1:noc
 
           xfirst = leftmargin + (j-1) * (subwidth + spacec);
           yfirst = bottommargin + ((nor - i + 1) - 1) * (subheight + spacer);
           
           if ismember(i, [4 8 12 16])
              positions{i/4,j} = [xfirst/plotwidth yfirst/plotheight subwidth/plotwidth subheight*4/plotheight];
           elseif i == 17
              positions{5,j} = [xfirst/plotwidth yfirst/plotheight subwidth/plotwidth subheight/plotheight]; 
           end
 
       end
    end
end
