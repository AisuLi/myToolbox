function timePoint = timeWinToPoint(xmin,timeBin)

		for iTrans = 1:size(timeBin,2)
	        timePoint(1,iTrans) = timeBin(iTrans)-xmin*1000+1;
	    end