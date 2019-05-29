function [myColormap] = makeColormap(num,isGrey,colorNum,colorName)

% size(colormap) = [64,3];

% written by Li Aisu 

if ~isGrey

	myColormap = ones(num,3);
	switch colorNum
		case 1
			step = 1/(num-1);

			switch colorName
				case 'b' 
					for i = 1:num
						myColormap(i,[1,2]) = 1-(i-1)*step;
					end
				case 'r'
					for i = 1:num
						myColormap(i,[2,3]) = 1-(i-1)*step;
					end					
			end

		case 2

			step       = 1/(num/2-1); 

			for i = 1:num/2 
				red(i,[1,2]) = 0+(i-1)*step;
			end

			blue  = flipud(red);

			myColormap([1:num/2],[1,2])     = red;
			myColormap([num/2+1:end],[2,3]) = blue;

			if ~isRedup

				red        = flipud(myColormap([1:num/2],:));
				blue       = flipud(myColormap([num/2+1:num],:));
				myColormap = [blue;red];
					
			end
	end

else
			
	step = 0.4/(num-1);

	for i = 1:num

		myColormap(i,[1,2,3]) = 0.6+(i-1)*step;
	end

			

end




% myColorMap         = (colorMap*transparency+bkColor*(1-transparency))./255;

