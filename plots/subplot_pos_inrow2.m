function positions = subplot_pos_inrow2(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, spacer, spacec, widthRatios)
    

    totalheight = plotheight - topmargin - bottommargin - spacer * (nor-1.0);
    totalwidth = plotwidth - leftmargin - rightmargin - spacec * (noc-1.0);
    
      
    positions{1,1} = [leftmargin/plotwidth bottommargin/plotheight totalwidth*widthRatios(1)/sum(widthRatios)/plotwidth totalheight/plotheight]; 
    positions{1,2} = [(leftmargin+spacec)/plotwidth+positions{1,1}(3) bottommargin/plotheight totalwidth*widthRatios(2)/sum(widthRatios)/plotwidth totalheight/plotheight]; 
    positions{1,3} = [(leftmargin+2*spacec)/plotwidth+positions{1,1}(3)+positions{1,2}(3) bottommargin/plotheight totalwidth*widthRatios(3)/sum(widthRatios)/plotwidth totalheight/plotheight];
   
end