function positions = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, spacer, spacec)
    
    subheight = (plotheight - topmargin - bottommargin - spacer * (nor - 1))/nor;
    subwidth = (plotwidth - leftmargin - rightmargin - spacec * (noc - 1))/noc;
    
    
    positions = cell(nor, noc);
    for i = 1:nor
       for j = 1:noc
 
           xfirst = leftmargin + (j - 1) * (subwidth + spacec);
           yfirst = bottommargin + ((nor - i + 1) - 1) * (subheight + spacer);

          positions{i,j} = [xfirst/plotwidth yfirst/plotheight subwidth/plotwidth subheight/plotheight];

       end
    end
end
