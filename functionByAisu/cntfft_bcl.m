
function [freBins,y,trialSNR] = cntfft_bcl(cntfilename,isEpoched,moreBinNum,stimTrigger,freRange,nfft,fftType,ldrType)

% this funtion is uesd to make FFT to cntFile
% ldrType is  noLDR,m1m2,m1m2Lps,m1m2leftLps,m1m2rightLps, Li Aisu,06-10-2017


trialSNR = [];

cnt = loadcnt_bcl(cntfilename); 


if ~exist('freRange','var')||isempty(freRange)
	freRange = [5 25];
end 


if 	~isEpoched 

	if ~exist('nfft','var')||isempty(nfft)
		nfft = 2^nextpow2(cnt.ldnsamples);
	end 
else
	if ~exist('nfft','var')||isempty(nfft)
		nfft = 2048;
	end 
end


if exist('ldrType','var')
	[transData] = LdrTransform_bcl(cnt.data,ldrType);
else
	transData = cnt.data;
end 


[selecData] = cntSelection_bcl(transData,cnt,isEpoched,moreBinNum,stimTrigger);




f = cnt.header.rate/2*linspace(0,1,nfft/2+1);
idxStart = find(f>freRange(1),1,'first');
idxEnd   = find(f<freRange(2),1,'last');
freBins  = f(idxStart:idxEnd);

if 	~isEpoched

		% y = fft(detrend(selecData'),nfft)/cnt.ldnsamples;

		y = fft(detrend(selecData'),nfft)/size(selecData,2);
		
		% if size(selecData,1)==1 % trigger-Test
		y = abs(double(y(idxStart:idxEnd,:)));
		% else            
				% y = abs(double(y(idxStart:idxEnd,[1:64])));
		% end

elseif isEpoched

	switch fftType
		case 'trialFFT'

				for iTrial = 1:size(selecData,3)

						itrialSignal     = selecData(:,:,iTrial); % c means current
						y                = fft(detrend(itrialSignal'),nfft)/size(selecData,2);
						y 				 = abs(double(y(idxStart:idxEnd,:)))
						yall(:,:,iTrial) = y;

						%%%%%%%%%%%%%%%%%   SNR of each Trial  %%%%%%%%%%%%%
		
						for ie  = 1:size(yall,2)       

							for iFrebin = snrskipNum+2:size(y,1)-(snrskipNum+1)
								trialSNR(iFrebin,ie,iTrial) = y(iFrebin,ie)/mean(y([iFrebin-(snrskipNum+1),iFrebin+(snrskipNum+1)],ie));
							end

						end	
						%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				end
				y = mean(abs(yall),3);


	    case 'avgFFT'

				avgData = mean(selecData,3);
				y   	= abs(fft(detrend(avgData'),nfft)/size(selecData,2));
				y       = y(idxStart:idxEnd,:);

				%------------ SNR of avg  ------------------//
				for ie  = 1:size(y,2)       
					for iFrebin = snrskipNum+2:size(yall,1)-(snrskipNum+1)
						SNRmean(iFrebin,ie) = ymean(iFrebin,ie)/mean(ymean([iFrebin-(snrskipNum+1),iFrebin+(snrskipNum+1)],ie));
					end
				end	
				%-------------------------------------------\\
	end

end

% e        = cnt.electloc;

% for ichan = 1:63
%     subplot(7,9,ichan);
% 	plot(freBins,2*abs(y(:,ichan)));
% 	xlabel('Frequency Hz');
% 	ylabel('|y(f)|');
% 	title(cnt.electloc(ichan).lab);
% end
