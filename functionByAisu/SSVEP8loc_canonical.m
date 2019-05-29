% function  = SSVEP8loc_canonical

% directory   = 'D:\expData\SSVEP_8loc\Exp1_20170328\Data\scanData\20180407\average';
% targetFre   = [8,20];
% T           = 1000./targetFre; % ms 
% T(2,:)      = T./2; 
% T(3,:)      = T(1,:)./3; 
% omega       = 2*pi./T;
% timeBin     = [-500 -1]; % ;201 700; 701 1200
% xmin        = -1;

% for iTimeBin = 1:size(timeBin,1)
% 	TimePointRange(iTimeBin,:) = timeWinToPoint(xmin,timeBin(iTimeBin,:));
% end

% %--------- make sin & cos ----------/
% baseX = 1:3000; 

% for iFre = 1:size(targetFre,2)

% 	for iWave = 1:size(T,1)

% 		sinX(:,iWave,iFre)   = sin(omega(iWave,iFre)*baseX);
% 		cosX(:,iWave,iFre)   = cos(omega(iWave,iFre)*baseX);		
	
% 	end
% end
% %----------------------------------\



% subNum        = [1:21];
% condition     = {'allcue13_2~45Hz','allcue21_2~45Hz','allcue34_2~45Hz','allcue42_2~45Hz'};
% [excludeChan,standard_eName]=getElectrodeNo_bcl('haveM1',{'M1','M2','VEO','HEO'});


% for iCondition = 1:size(condition,2)

% 		cd(fullfile(directory,condition{iCondition}));

% 		for iSub = subNum

% 				iSubFileName = [num2str(iSub),'-',condition{iCondition},'.avg'];

% 			    [signal,chan_names,variance, pnts, rate, xmin, xmax,nsweeps] = loadavg_bcl(iSubFileName);

% 			    signal(:,excludeChan) = [];
% 			    for iTimeBin = 1:size(timeBin,1)

% 						cTimePoint = TimePointRange(iTimeBin,:);
% 						cSignal    = signal([cTimePoint(1):cTimePoint(2)],:);

% 						for iFre = 1:size(targetFre,2)

% 								for iWave = 1:size(sinX,2)
							
% 									[A1,B1,R1,U1,V1] = canoncorr(cSignal, sinX([1:size(cSignal,1)],iWave,iFre));

% 									[A2,B2,R2,U2,V2] = canoncorr(cSignal, cosX([1:size(cSignal,1)],iWave,iFre));

% 									allR(iSub,iWave,iFre,iCondition,iTimeBin) = sqrt(R1^2+R2^2);
% 									allA1(:,iSub,iWave,iFre,iCondition,iTimeBin)= A1;
% 									allU1(:,iSub,iWave,iFre,iCondition,iTimeBin)= U1;
% 									allA2(:,iSub,iWave,iFre,iCondition,iTimeBin)= A2;
% 									allU2(:,iSub,iWave,iFre,iCondition,iTimeBin)= U2;
% 								end
% 						end

% 						clear  cTimePoint cSignal;
% 				end
% 		end

% 		cd('..');
% end 

% allR_minusBase(:,:,:,:,1) = allR(:,:,:,:,2)-allR(:,:,:,:,1);
% allR_minusBase(:,:,:,:,2) = allR(:,:,:,:,3)-allR(:,:,:,:,1);

% %--------shift-------------/
% allR_minusBase(:,:,1,:,:) =  allR_minusBase(:,:,1,[1 3 2 4],:);
% allR_minusBase(:,:,2,:,:) =  allR_minusBase(:,:,2,[4 2 3 1],:);
% %--------------------------\


% allR_minusBase(:,:,3,:,:) = mean(allR_minusBase,3); % mean of 8hz & 20hz
% allR_minusBase(:,:,:,:,3) = mean(allR_minusBase,5); % mean of 201To700 & 701To1200


% allR_minusBase(size(allR_minusBase,1)+1,:,:,:,:) = mean(allR_minusBase,1);


% for iTimeRange = 1:size(allR_minusBase,5)
% 	for iFre = 1:size(allR_minusBase,3)
% 		for iWave = 1:size(allR_minusBase,2)

% 			cR = squeeze(allR_minusBase(:,iWave,iFre,:,iTimeRange));
% 	 		for iCmp = 1:3

% 	 			[H,pvalue(iCmp,iWave,iFre,iTimeRange),CI] = ttest(cR(:,1),cR(:,iCmp+1));

% 	 		end
% 	 	end
% 	end
% end






% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%  compute postBaseline[201 700; 701 1200] based on baseline's weight to get time course ----------

directory   = 'D:\expData\SSVEP_8loc\Exp1_20170328\Data\scanData\20180407\average';
targetFre   = [8,20];
T           = 1000./targetFre; % ms 
T(2,:)      = T./2; 
T(3,:)      = T(1,:)./3; 
omega       = 2*pi./T;
timeBin     = [-500 -1; 201 700 ;701 1200]; % ;201 700; 701 1200
xmin        = -1;

for iTimeBin = 1:size(timeBin,1)
	TimePointRange(iTimeBin,:) = timeWinToPoint(xmin,timeBin(iTimeBin,:));
end

%--------- make sin & cos ----------/
baseX = 1:3000; 

for iFre = 1:size(targetFre,2)

	for iWave = 1:size(T,1)

		sinX(:,iWave,iFre)   = sin(omega(iWave,iFre)*baseX);
		cosX(:,iWave,iFre)   = cos(omega(iWave,iFre)*baseX);		
	
	end
end
%----------------------------------\



subNum        = [1:21];
condition     = {'allcue13_2~45Hz','allcue21_2~45Hz','allcue34_2~45Hz','allcue42_2~45Hz'};
[excludeChan,standard_eName]=getElectrodeNo_bcl('haveM1',{'M1','M2','VEO','HEO'});


for iCondition = 1:size(condition,2)

		cd(fullfile(directory,condition{iCondition}));

		for iSub = subNum

				iSubFileName = [num2str(iSub),'-',condition{iCondition},'.avg'];

			    [signal,chan_names,variance, pnts, rate, xmin, xmax,nsweeps] = loadavg_bcl(iSubFileName);

			    signal(:,excludeChan) = [];


				baseTimePoint = TimePointRange(iTimeBin,:);
				baseSignal    = signal([baseTimePoint(1):baseTimePoint(2)],:);


				for iFre = 1:size(targetFre,2)

						for iWave = 1:size(sinX,2)
					
					        %------ baseline canoncorr ------/
							[A1,B1,R1,U1,V1] = canoncorr(baseSignal, sinX([1:size(baseSignal,1)],iWave,iFre)); % 1:sin

							[A2,B2,R2,U2,V2] = canoncorr(baseSignal, cosX([1:size(baseSignal,1)],iWave,iFre)); % 2:cos

							allR(iSub,iWave,iFre,iCondition,iTimeBin) = sqrt(R1^2+R2^2);
							allA1(:,iSub,iWave,iFre,iCondition,iTimeBin)= A1;
							allA2(:,iSub,iWave,iFre,iCondition,iTimeBin)= A2;
 							%------ baseline canoncorr ------\

							for iTimeBin = 2:size(timeBin,1) %%%%%% attention %%%%%% 

									timePoint_postBase = TimePointRange(iTimeBin,:);
									signal_postBase    = signal([timePoint_postBase(1):timePoint_postBase(2)],:);

								       % U = (X - repmat(mean(X),N,1))*A
								       % V = (Y - repmat(mean(Y),N,1))*B
								    ts_sin = (signal_postBase - repmat(mean(signal_postBase),size(signal_postBase,1),1))*A1;
								    ts_cos = (signal_postBase - repmat(mean(signal_postBase),size(signal_postBase,1),1))*A2;
									% allU1_postBase(:,iSub,iWave,iFre,iCondition,iTimeBin-1) = ts_sin;
									% allU2_postBase(:,iSub,iWave,iFre,iCondition,iTimeBin-1) = ts_cos;

									for i = 1:size(ts_sin,1)
										ts_sincos(i) = sqrt(ts_sin(i)^2+ts_cos(i)^2);
									end 

									allTs_sincos(:,iSub,iWave,iFre,iCondition,iTimeBin-1) = ts_sincos;
							end
						end
				end

		end

		cd('..');
end 



%------- read rawFreData (not harmonic), AND THEN adjusted conditions -----------/
freSignal = squeeze(allTs_sincos(:,:,1,:,:,:));  
signal    = cat(3,freSignal(:,:,1,[1 3 2 4],:),freSignal(:,:,2,[4 2 3 1],:)); 
% size: 1-timeCourse,2-sub,3-fre,4-condition,5-timeRange
%------------------- have adjusted conditions, the end --------------------------\





%%---do fft of post-baseline for 4 conditions of 2 fre,separately---/
freRange     = [5 25];
targetFre    = [8;20]; 
timeInterval = 'all';
windowType   = 'cosin';
ldrTransType = 'noldr';
snrskipNum   = 1;

for iTimeRange = 1:size(signal,5)

	for iCondition = 1:size(signal,4)

		for iFre = 1:size(signal,3)

			for iSub = 1:size(signal,2)

				iSubFileName =  squeeze(signal(:,iSub,iFre,iCondition,iTimeRange));

				[powerData,f,chan_names,chanList,rawSignal] = FFTanalysis_bcl_su(iSubFileName,timeInterval,[],[],windowType,ldrTransType);

				allPower_4con_allTimeRange(:,iSub,iCondition,iFre,iTimeRange)  = powerData; % size: 1-freBins;  2-sub; 3-condition; 4-freStim; 5-timeRange(base,postBase)

				[SNR] = computeSNR(powerData,snrskipNum);

				allSNR_4con_allTimeRange(:,iSub,iCondition,iFre,iTimeRange)    = SNR;

				clear SNR powerData;

			end
		end
	end
end 

%---------do fft for 4 conditions of 2 fre,separately, the end--------\

[freIdx]                   = cutFreBins(f,freRange); 
allSNR_4con_allTimeRange   = allSNR_4con_allTimeRange(freIdx,:,:,:,:);
allPower_4con_allTimeRange = allPower_4con_allTimeRange(freIdx,:,:,:,:);
freBins                    = f(freIdx);


allSNR_4con_allTimeRange(:,size(allSNR_4con_allTimeRange,2)+1,:,:,:)     = mean(allSNR_4con_allTimeRange,2); 
allPower_4con_allTimeRange(:,size(allPower_4con_allTimeRange,2)+1,:,:,:) = mean(allPower_4con_allTimeRange,2); 

%-------find closest frePoint-------/
for iFre = 1:size(targetFre,1)
	[X, freIndex]         = sort(abs(freBins-targetFre(iFre,1)),'ascend');
	freTimePoint(iFre,1)  = freIndex(1);
	realFre(iFre,1) 	  = freBins(freIndex(1));
end
%-------have found, the end---------\


allSNR_4con_cano   = allSNR_4con_allTimeRange(:,:,:,:,2);
allPower_4con_cano = allPower_4con_allTimeRange(:,:,:,:,2);

subNumPerFig = 5;
figNum       = ceil(size(allSNR_4con_cano,2)./subNumPerFig);

%-------------------------plot  SNR & Power----------------------
for iPage = 1:figNum
	figure;
	set(gcf,'Position',get(0,'ScreenSize'),'color','w');
	suptitle(['canoncorr     left:    ',num2str(realFre(1,:)),'Hz     con1		con2   con3   con4     ',    '          right:  '   num2str(realFre(2,:))  'Hz     con1		con2   con3   con4     ']); 

	cSub = [1:subNumPerFig]+(iPage-1)*subNumPerFig;

	for iSub = cSub

		for iFre = 1:size(allSNR_4con_cano,4)

			for iCondition = 1:size(allSNR_4con_cano,3)
				
				if iSub<=size(allSNR_4con_cano,2)
					%-------plot topo-------/
					if iFre==1
						figureLoc = 9*(iSub-1-(iPage-1)*subNumPerFig)+iCondition;
					elseif iFre==2
						figureLoc = 9*(iSub-1-(iPage-1)*subNumPerFig)+iCondition+5;
					end	
					subplot(subNumPerFig,9,figureLoc);

					%%%%%%size: 1-freBins;  2-sub; 3-condition; 4-freStim; 
					[AX,H1,H2] = plotyy(freBins,allPower_4con_cano(:,iSub,iCondition,iFre),freBins,allSNR_4con_cano(:,iSub,iCondition,iFre),'plot');

					box off;

					set(AX(1),'XColor','k','YColor','b','Ylim',[0,0.4],'ytick',[0:0.1:0.4],'Xlim',freRange);
					set(AX(2),'XColor','k','YColor','r','Ylim',[0,4],'ytick',[0:1:4],'Xlim',freRange);
					set(H1,'LineStyle','-','color','b');
					set(H2,'LineStyle','-','color','r');
					line(repmat(realFre(iFre,:),1,2),[0,.2],'linestyle',':','linewidth',2, 'color',[0 0 0 ]); %


					% line(repmat(realFre(1,:),1,2),[0,.2],'linestyle',':','linewidth',2, 'color',[0 0 0 ]); %
					% line(repmat(realFre(2,:),1,2),[0,.2],'linestyle',':','linewidth',2, 'color',[0 0 0 ]); %
					% line(repmat(realFre(3,:),1,2),[0,.2],'linestyle',':','linewidth',2, 'color',[0 0 0 ]); %
					% line(repmat(realFre(4,:),1,2),[0,.2],'linestyle',':','linewidth',2, 'color',[0 0 0 ]); %
					% line(repmat(realFre(5,:),1,2),[0,.2],'linestyle',':','linewidth',2, 'color',[0 0 0 ]); %
					% line(repmat(realFre(6,:),1,2),[0,.2],'linestyle',':','linewidth',2, 'color',[0 0 0 ]); %
					end
			end
		end

		subplot(subNumPerFig,9,9*(iSub-1-(iPage-1)*subNumPerFig)+5);
		if iSub<size(allPower_4con_cano,2)
				title(['sub:' num2str(iSub)]);
		elseif iSub==size(allPower_4con_cano,2)
				title('allSubMean');
				break;
		end

	end
	l=legend([H1,H2],{'power';'SNR'});
	legend('boxoff');
	set(l,'Position',[0.85199652777777 0.539659128154703 0.0479166666666667 0.0412979351032448]);
	set(gcf,'color','w');
	% saveas(gcf,['S',num2str(cSub(1)),'~',num2str(cSub(end)),'_',fftType,'_SNRtopo.pdf']);
end
%------------------------------\\

