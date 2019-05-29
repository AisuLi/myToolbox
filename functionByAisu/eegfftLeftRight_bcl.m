function  [freBins,ymean,fftType] = eegfftLeftRight_bcl(eegName,freRange,nfft,leftStiType,rightStiType,ldrType,fftType)


if ~exist('freRange','var')||isempty(freRange)
	freRange = [5 25];
end 


if ~exist('nfft','var')||isempty(nfft)
	nfft = 2048;
end 


[signal, accept, typeeeg, rt, response, chan_names, pnts, ...
ntrials, srate, xmin, xmax] = loadeeg(eegName);

signal   = reshape(signal,[size(signal,1),pnts,ntrials]);

% leftStiType=[1:2:40];
% rightStiType=[2:2:40];

if exist('leftStiType','var')
		
		leftTrial=[];
		rightTrial=[];

		for iTrial=1:ntrials
			if ismember(typeeeg(iTrial),leftStiType)
				cleftTrial=iTrial;
				leftTrial=[leftTrial cleftTrial];
			elseif ismember(typeeeg(iTrial),rightStiType)
				crightTrial=iTrial;
				rightTrial=[rightTrial crightTrial];
			end
		end

		leftSignal=signal(:,:,leftTrial);
		rightSignal=signal(:,:,rightTrial);

end 



if exist('ldrType','var')

	if size(ldrType,2)==1
		[transData] = LdrTransform_bcl(signal,ldrType);
	elseif size(ldrType,2)==2

		leftTrans  = LdrTransform_bcl(leftSignal,ldrType{1});
		rightTrans = LdrTransform_bcl(rightSignal,ldrType{2});
		transData  = cat(3,leftTrans,rightTrans);
	end

end 


f        = srate/2*linspace(0,1,nfft/2+1);
idxStart = find(f>freRange(1),1,'first');
idxEnd   = find(f<freRange(2),1,'last');
freBins  = f(idxStart:idxEnd);


switch fftType
	case 'trialFFT'

			for iTrial=1:ntrials

					itrialSignal     = transData([1:64],:,iTrial);
					y                = fft(detrend(itrialSignal'),nfft)/pnts;
					yall(:,:,iTrial) = y(idxStart:idxEnd,:);

			end
			ymean = mean(abs(yall),3);

    case 'avgFFT'

			avgData = mean(transData,3)
			y       = abs(fft(detrend(avgData'),nfft)/pnts);
			ymean   = y(idxStart:idxEnd,:);
end