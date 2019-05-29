function  [avgData,freBins,Power,SNR] = eegfftSNR_bcl(eegName,freRange,timeBin,stimTrigger,fftType,nfft)


if ~exist('freRange','var')||isempty(freRange)
	freRange = [3 45];
end 


if ~exist('nfft','var')||isempty(nfft)
	nfft = [];
end 


snrskipNum = 1;

if isa(eegName,'char')
	[signal, accept, typeeeg, rt, response, chan_names, pnts, ...
	ntrials, srate, xmin, xmax] = loadeeg(eegName);

	timePoint     = timeWinToPoint(xmin,timeBin);
	signal        = reshape(signal,[size(signal,1),pnts,ntrials]);
	signal        = signal(:,[timePoint(1):timePoint(end)],ismember(typeeeg,stimTrigger'));
elseif isa(eegName,'double')
	signal  	  = eegName;
	fftType  	  = 'avgFFT';
	srate         = 1000;
end	

timeBinLength = size(signal,2);

if timeBinLength>2048
	nfft = 2^nextpow2(timeBinLength);
else
	nfft = 1024;
end

f        = srate/2*linspace(0,1,nfft/2+1);

SNRmean  = zeros(size(f,2),size(signal,1));

switch fftType

%%% should be revised later .....
	% case 'trialFFT'
 %        trialSNR = zeros(size(f,2),size(signal,1),ntrials);
	% 	for iTrial=1:ntrials

	% 			itrialSignal         = signal(:,:,iTrial);
	% 			Power                = fft(detrend(itrialSignal'),nfft)/timeBinLength;
	% 			Power                = abs(double(Power(idxStart:idxEnd,:)));
	% 			Powerall(:,:,iTrial) = Power ;

	% 			%%%%%%%%%%%%%%%%%   SNR of each Trial  %%%%%%%%%%%%%
	% 			for ie  = 1:size(Powerall,2)       
	% 				for iFrebin = snrskipNum+2:size(freBins,2)-(snrskipNum+1)
	% 						trialSNR(iFrebin,ie,iTrial) = Power(iFrebin,ie)/mean(Power([iFrebin-(snrskipNum+1),iFrebin+(snrskipNum+1)],ie));
	% 				end
	% 			end	
	% 			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% 			clear Power ;
	% 	end
	% 	avgData   = [];
	% 	SNRmean   = squeeze(mean(trialSNR,3));
	% 	Powermean = squeeze(mean(abs(Powerall),3));

    case 'avgFFT'

		avgData   = mean(signal,3);
		Power     = abs(fft(detrend(avgData'),nfft)/timeBinLength);
		%------------ SNR of avg  ------------------//
		for ie  = 1:size(Power,2)       
			for iFrebin = snrskipNum+2:size(f,2)-(snrskipNum+1)
				SNR(iFrebin,ie) = Power(iFrebin,ie)/mean(Power([iFrebin-(snrskipNum+1),iFrebin+(snrskipNum+1)],ie));
			end
		end	
		%-------------------------------------------\\

end


idxStart = find(f>freRange(1),1,'first');
idxEnd   = find(f<freRange(2),1,'last');
freBins  = f(idxStart:idxEnd);
Power    = Power(idxStart:idxEnd,:);
SNR      = SNR(idxStart:idxEnd,:);