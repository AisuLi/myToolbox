function  [powerData,freBins,chan_names,chanList,rawSignal] = FFTanalysis_bcl_su(AvgFilename,timeInterval,chanList,FFTpoints,windowType,ldrTransType)



if ~exist('chanList','var')||isempty(chanList)
	chanList = 'all';
end 



if ~exist('FFTpoints','var')||isempty(FFTpoints)
	FFTpoints = 1024;
end 


if ~exist('windowType','var')||isempty(windowType)
	windowType = 'cosin';
end 





% Rev. by yang zhang Tue Apr  3 14:53:54 2018
%      Soochow University, China

if isa(AvgFilename,'char')

	[signal,chan_names,variance,pnts,fs,xmin,xmax]=loadavg_bcl(AvgFilename,chanList);
    
elseif isa(AvgFilename,'double')
	signal  	  = AvgFilename;
    fs = 1000;
    chan_names = 'oneChan';
end	



isDetrend = false;

if isDetrend
	signal = detrend(signal);
end



if isa(timeInterval,'double')

		binRange     = binTimeTransf(timeInterval,round(xmin*1000),fs,1); % time to bin
		signal       = signal(binRange(1):binRange(2),:);% each column for each chan

else
	fprintf('%s\n','------You r getting the whole TimeWindow-------');
end
		



rawSignal    = signal;
rawLengthBin = size(signal,1);



[ avgLdr  lpsLdr ] = LdrTransform_bcl(0,'noLDR');


switch ldrTransType

	case 'lpsldr' 

		 signal = (lpsLdr*signal')';

	case 'avgldr'

		 signal = (avgLdr*signal')';
	case 'noldr'

		signal = signal;

end

% add ldrTransType,Li Aisu




%---- added window function -------/
winData = getWindowData(windowType,rawLengthBin);

signal = signal.* repmat(winData(:),1,size(signal,2));
%----------------------------------\


%--- force the data bin to be uneven ---/
if floor(rawLengthBin/2) == (rawLengthBin/2)
	isEven = true;
else
	isEven = false;
end
%---------------------------------------\

Y         = fft(signal,FFTpoints);
powerData = abs(Y)/rawLengthBin;

lengthBin = FFTpoints;
%----- double the power ----/
if isEven
	powerData            = powerData(1:(lengthBin/2),:);
	powerData(2:end-1,:) = powerData(2:end-1,:)*2;
else
	powerData            = powerData(1:(lengthBin/2+1),:);
	powerData(2:end,:)   = powerData(2:end,:)*2;
end
%---------------------------\

freBins   = fs*(0:(lengthBin/2))/lengthBin;


end % function end




function winData = getWindowData(windowType,rawLengthBin)

	switch windowType
		case 'cosin'
			winData = cos((0:rawLengthBin-1).*((2*pi)/rawLengthBin) - pi)./2+0.5;

% 			winData = evalc('window(@',windowType,',rawLengthBin);');
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



