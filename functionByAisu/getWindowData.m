

function winData = getWindowData(windowType,rawLengthBin)

	switch windowType
		case 'cosin'
			winData = cos((0:rawLengthBin-1).*((2*pi)/rawLengthBin) - pi)./2+0.5;

			winData = evalc('window(@',windowType,',rawLengthBin);');
		case 'barthannwin'
			winData = evalc('window(@',windowType,',rawLengthBin);');
		case 'bartlett'
			winData = evalc('window(@',windowType,',rawLengthBin);');
		case 'blackman'
			winData = evalc('window(@',windowType,',rawLengthBin);');
		case 'blackmanharris'
			winData = evalc('window(@',windowType,',rawLengthBin);');
		case 'bohmanwin'
			winData = evalc('window(@',windowType,',rawLengthBin);');
		case 'chebwin'
			winData = evalc('window(@',windowType,',rawLengthBin);');
		case 'flattopwin'
			winData = evalc('window(@',windowType,',rawLengthBin);');
		case 'gausswin'
			winData = evalc('window(@',windowType,',rawLengthBin);');
		case 'hamming'
			winData = evalc('window(@',windowType,',rawLengthBin);');
		case 'hann'
			winData = evalc('window(@',windowType,',rawLengthBin);');
		case 'kaiser'
			winData = evalc('window(@',windowType,',rawLengthBin);');
		case 'nuttallwin'
			winData = evalc('window(@',windowType,',rawLengthBin);');
		case 'parzenwin'
			winData = evalc('window(@',windowType,',rawLengthBin);');
		case 'rectwin'
			winData = evalc('window(@',windowType,',rawLengthBin);');
		case 'taylorwin'
			winData = evalc('window(@',windowType,',rawLengthBin);');
		case 'triang'
			winData = evalc('window(@',windowType,',rawLengthBin);');
		case 'tukeywin'
			winData = evalc('window(@',windowType,',rawLengthBin);');
		case 'null'
			winData = ones(1,rawLengthBin);
		otherwise
			error('wrong windowType !');
	end 


end 
