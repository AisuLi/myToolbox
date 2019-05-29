
function [trialNum trialBaseR trialRealR avgBaseR avgRealR] = canonical2_eeg(directory,fileFilt,stimTrigger,needTimeBin)

% written by Li Aisu, 2017-12-03

%%%% add harmonic wave (iwave) based on  canonical_eeg, 2017-12-05


%%%%%%%%%%%%%%%%%%%%%fake Data %%%%%%%%%%%%%%%%%
% x = [0:1/1000:1];
% noise = rand(size(x,2),1);

% for iSignal = 1:200

% 	sinY  = iSignal*sin(x)';
% 	cosY  = iSignal*cos(x)';

% 	 for inoise = 1:100
% 			realY=sinY+cosY+inoise*noise;
% 			[A,B,R]=canoncorr(realY, [sinY cosY]);
% 			r(iSignal,inoise)=R;
% 	 end

% 	 clear sinY cosY 
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% needTimeBin = [200:1200];
 


allfre              = [15,12,8,20];
T                   = 1000./allfre; % ms 
T(2,:)              = T./2; 
T(3,:)				= T(1,:)./3; 
omega               = 2*pi./T;


allbaseTimeWinStart = [-400, -500, -500,  -525]; % according to T, caculate the N*T(N=2,3,4,5,6)

% freNum              = find(allfre==fre);
% baseTimeWinStart    = allbaseTimeWinStart(freNum);

eegFile = dir(fullfile(directory,fileFilt));

for iSub = 1:size(eegFile,1)

		eegName = eegFile(iSub).name;
		
		[allsignal, accept, typeeeg, rt, response, chan_names, pnts, ntrials, srate, xmin, xmax] = loadeeg(eegName);


		allsignal  = reshape(allsignal,[size(allsignal,1),pnts,ntrials]);
		
		signal     = allsignal(:,:,ismember(typeeeg,stimTrigger'));

		trialNum(iSub,1)   = size(signal,3);
		
		for iFre = 1:size(allfre,2)

				basePointRange = timeWinToPoint(xmin,[allbaseTimeWinStart(iFre):-1]);
				baseSignal = signal(:,basePointRange,:);
			    realSignal = signal(:,timeWinToPoint(xmin,needTimeBin),:);

				baseX      = [1:size(baseSignal,2)]';
				realX      = [1:size(realSignal,2)]';



				for iWave = 1:size(T,1)


					sinBaseX   = sin(omega(iWave,iFre)*baseX);
					cosBaseX   = cos(omega(iWave,iFre)*baseX);		
					sinRealX   = sin(omega(iWave,iFre)*realX);
					cosRealX   = cos(omega(iWave,iFre)*realX);		


					%%%%%%%%%%%%%%%%%%   trial corr   %%%%%%%%%%%%%%%
						for iEpoch = 1:size(signal,3)

							iEpochBaseSignal = squeeze(baseSignal(:,:,iEpoch))';
							iEpochRealSignal = squeeze(realSignal(:,:,iEpoch))';

							[ibaseA,ibaseB,ibaseR] = canoncorr(iEpochBaseSignal, [sinBaseX cosBaseX]); % ,baseU(:,:,iEpoch),baseV(:,:,iEpoch)
							baseR(1,iEpoch)        = sqrt(ibaseR(1)^2+ibaseR(2)^2);

							[irealA,irealB,irealR] = canoncorr(iEpochRealSignal, [sinRealX cosRealX]);  %,realU(:,:,iEpoch),realV(:,:,iEpoch)
							realR(1,iEpoch) 	   = sqrt(irealR(1)^2+irealR(2)^2);

							clear  iEpochBaseSignal iEpochRealSignal

						
						end

						trialBaseR(iSub,iWave,iFre) = mean(baseR,2);
						% trialBaseU(iSub,iWave,iFre) = mean(baseU,3);
						% trialBaseV(iSub,iWave,iFre) = mean(baseV,3);	
											
						trialRealR(iSub,iWave,iFre) = mean(realR,2);
						% trialRealU(iSub,iWave,iFre) = mean(realU,3);
						% trialRealV(iSub,iWave,iFre) = mean(realV,3);	
					%%%%%%%%%%%%%%%%%%   avg corr   %%%%%%%%%%%%%%%
						avgBaseSignal = squeeze(mean(baseSignal,3))';
						avgRealSignal = squeeze(mean(realSignal,3))';
						
						[iavgBaseA,iavgBaseB,iavgBaseR] = canoncorr(avgBaseSignal, [sinBaseX cosBaseX]);
						avgBaseR(iSub,iWave,iFre)       = sqrt(iavgBaseR(1)^2+iavgBaseR(2)^2);
						
						[iavgRealA,iavgRealB,iavgRealR] = canoncorr(avgRealSignal, [sinRealX cosRealX]);
						avgRealR(iSub,iWave,iFre)       = sqrt(iavgRealR(1)^2+iavgRealR(2)^2);
				end
		end			
end

