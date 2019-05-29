
%% avgFFT 
clear all;
directory  = 'D:\expData\SSVEP_8loc\Exp1_20170328\Data\scanData\20180407\average';
chanList   = 'all';
subNum     = [1:21];
FFTpoints  = 2048;
windowType = 'cosin';
snrskipNum = 1;
freRange   = [5 25];
targetFre  = [8;20]; 
ldrTransType = 'noldr';


%%%% step1
%-------------- do fft of [-500 1200] for all conditions-------------/
timeBin    = [-1000,2000];
maxChanNum = 5;
condition  = {'all_2~45Hz'};



for iSub = subNum

	iSubFileName = [num2str(iSub),'-',condition{:},'.avg'];
	
	[powerData,f,chan_names,chanList,rawSignal] = FFTanalysis_bcl_su(iSubFileName,timeBin,chanList,FFTpoints,windowType,ldrTransType);

	allPower(:,:,iSub) = powerData;

	[SNR] = computeSNR(powerData,snrskipNum);

	allSNR(:,:,iSub)  = SNR;  
	fprintf('%d\n', iSub);
	clear powerData SNR

end

[freIdx]   = cutFreBins(f,freRange); 
allPower   = allPower(freIdx,:,:);
allSNR     = allSNR(freIdx,:,:);
freBins    = f(freIdx);

%-------find closest frePoint-------/
for iFre = 1:size(targetFre,1)
	[X, freIndex]         = sort(abs(freBins-targetFre(iFre,1)),'ascend');
	freTimePoint(iFre,1)  = freIndex(1);
	realFre(iFre,1) 	  = freBins(freIndex(1));
end
%-----------------------------------\

frePointSNR   = allSNR(freTimePoint(:),:,:);
frePointPower = allPower(freTimePoint(:),:,:);
%---------- do fft of [-500 1200] for all conditions,the end---------\ 



[none,standard_eName] = getElectrodeNo_bcl('haveM1',[]);


%%%%% step2
%%%%%%%%%%% method-1: choose eachSub's max chans and sort the results %%%%%%%%%%%%%%%%%%%
%%--------plot topo and wave of SNR  to check chosen chans ---------//

for iSub = 1:size(frePointSNR,3)
	for iFre = 1:size(targetFre,1)
				
		[X, chanIndex]                  = sort(frePointSNR(iFre,[1:size(standard_eName,2)-2],iSub),'descend');	

		% chanIndex                       = chanIndex+23;					
		chooseChanSNR(iSub,:,iFre)      = chanIndex(1:maxChanNum);
		frePointMaxChanSNR(iFre,:,iSub) = frePointSNR(iFre,chanIndex(1:maxChanNum),iSub);	
		MaxChanSNR_SNR(:,:,iSub,iFre)   = allSNR(:,chanIndex(1:maxChanNum),iSub);
		MaxChanSNR_Power(:,:,iSub,iFre) = allPower(:,chanIndex(1:maxChanNum),iSub);
		clear X chanIndex

		[X, chanIndex]                    = sort(frePointPower(iFre,[1:size(standard_eName,2)-2],iSub),'descend');		

		% chanIndex                         = chanIndex+23;				
		chooseChanPower(iSub,:,iFre)      = chanIndex(1:maxChanNum);	
		frePointMaxChanPower(iFre,:,iSub) = frePointPower(iFre,chanIndex(1:maxChanNum),iSub);		
		MaxChanPower_SNR(:,:,iSub,iFre)   = allSNR(:,chanIndex(1:maxChanNum),iSub);
		MaxChanPower_Power(:,:,iSub,iFre) = allPower(:,chanIndex(1:maxChanNum),iSub);
		clear X chanIndex	
	end
end


%%%%%%%%%%%  method-2: choose chans according to subavg's max chans

subavgPower = squeeze(mean(frePointPower(:,:,:),3));
subavgSNR   = squeeze(mean(frePointSNR(:,:,:),3));

for iFre = 1:2
		[X, chanIndex]                   = sort(subavgSNR(iFre,[1:size(standard_eName,2)-2]),'descend');
		chooseChanSubavgSNR(iFre,:)      = chanIndex;

		[X, chanIndex]                   = sort(subavgPower(iFre,[1:size(standard_eName,2)-2]),'descend');
		chooseChanSubavgPower(iFre,:)    = chanIndex;
end     

save(['avgFFT_mergeAllconFindChan_-1000To2000','noldr_skip' num2str(snrskipNum)]); 

[myColormap] = makeColormap(20,0,1,'r');
figure;
for iFre = 1:2
	subplot(1,2,iFre);
	topoplot_bcl(subavgSNR(iFre,:),'EEGChans',standard_eName,'maplimitsDouble',[0 2.5],'colormap',myColormap);
	hold on;
end

[myColormap] = makeColormap(20,0,1,'b');
figure;
for iFre = 1:2
	subplot(1,2,iFre);
	topoplot_bcl(subavgPower(iFre,:),'EEGChans',standard_eName,'maplimitsDouble',[0 0.1],'colormap',myColormap);
	hold on;
end








%%%%%  
%%%%% step2:sort the SNR  
%%---------sort the SNR of allSub according to eachSub's eAvg-----------/
%%%%-----compute eAvg--------
frePointMaxChanSNR(:,maxChanNum+1,:) = mean(frePointMaxChanSNR,2);

for iFre = 1:size(frePointMaxChanSNR,1)
		[X, subIndex]                         = sort(frePointMaxChanSNR(iFre,maxChanNum+1,:),'descend');
		sortSubSNR(iFre,:)                    = subIndex;
		frePointSNR_sorted(iFre,:,:)          = frePointSNR(iFre,:,subIndex);
		MaxChanSNR_SNR_sorted(:,:,:,iFre)     = MaxChanSNR_SNR(:,:,subIndex,iFre);
		MaxChanSNR_Power_sorted(:,:,:,iFre)   = MaxChanSNR_Power(:,:,subIndex,iFre);
		chooseChanSNR_sorted(:,:,iFre)   	  = chooseChanSNR(subIndex,:,iFre);	
		clear X subIndex;
end

frePointSNR_sorted(:,:,size(sortSubSNR,2)+1)        = mean(frePointSNR_sorted,3);
MaxChanSNR_SNR_sorted(:,:,size(sortSubSNR,2)+1,:)   = mean(MaxChanSNR_SNR_sorted,3);
MaxChanSNR_SNR_sorted(:,maxChanNum+1,:,:)           = mean(MaxChanSNR_SNR_sorted,2);
MaxChanSNR_Power_sorted(:,:,size(sortSubSNR,2)+1,:) = mean(MaxChanSNR_Power_sorted,3);
MaxChanSNR_Power_sorted(:,maxChanNum+1,:,:)         = mean(MaxChanSNR_Power_sorted,2);

waveNum = maxChanNum+1;
[myColormap] = makeColormap(20,0,1,'r');
subNumPerFig = 5;
figNum       = ceil(22./subNumPerFig);
% [P1eNo,standard_eName]=getElectrodeNo_bcl('haveM1',[]);

for iPage = 1:figNum
	figure;
	set(gcf,'Position',get(0,'ScreenSize'),'color','w');

	% suptitle([fftType,'SNR     -500To1200    chooseChan   left:    ', num2str(realFre(1)) 'Hz          right:    ',num2str(realFre(2))  'Hz']);
	
	cSub = [1:subNumPerFig]+(iPage-1)*subNumPerFig;

	for iSub = cSub
				
		for iFre = 1:size(targetFre,1)
				
			if iSub<=size(chooseChanSNR,1)+1

				
				%-------plot topo-------/
				if iFre==1
					figureLoc = 2*(waveNum+1)*(iSub-1-(iPage-1)*subNumPerFig)+1;
				elseif iFre==2
					figureLoc = 2*(waveNum+1)*(iSub-1-(iPage-1)*subNumPerFig)+waveNum+2;
				end	

				subplot(subNumPerFig,2*(waveNum+1),figureLoc);
				if iSub<= size(chooseChanSNR,1)
					topoplot_bcl(frePointSNR_sorted(iFre,:,iSub),'EEGChans',standard_eName,'maplimitsDouble',[0 4],'colormap',myColormap);
				end
				%-----------------------\	

				for iWave = 1:waveNum
					figureLoc  = figureLoc+1;
					subplot(subNumPerFig,2*(waveNum+1),figureLoc);

					[AX,H1,H2] = plotyy(freBins,MaxChanSNR_Power_sorted(:,iWave,iSub,iFre),freBins,MaxChanSNR_SNR_sorted(:,iWave,iSub,iFre),'plot');
					
					if iSub<=size(chooseChanSNR,1)
						if iWave<=maxChanNum
							title([num2str(sortSubSNR(iFre,iSub)),standard_eName(squeeze(chooseChanSNR_sorted(iSub,iWave,iFre)))]);
						else
							title(['eAvg']);					
						end
					else
						title(['allSubAvg']);
					end
					
					box off;

					set(AX(1),'XColor','k','YColor','b','Ylim',[0,0.2],'ytick',[0:0.05:0.2],'Xlim',freRange);		 				
					set(AX(2),'XColor','k','YColor','r','Ylim',[0,4],'ytick',[0:1:4],'Xlim',freRange);
					set(H1,'LineStyle','-','color','b');
					set(H2,'LineStyle','-','color','r');
					line(repmat(realFre(iFre,:),1,2),[0,.2],'linestyle',':','linewidth',2, 'color',[0 0 0 ]); %
					line([5 25],[0.05,0.05],'linestyle',':','linewidth',1, 'color',[0 0 0 ]); %
				end
			end
		end
	end
	legend([H1,H2],{'power';'SNR'},'Location','North');
	legend('boxoff');
	set(gcf,'color','w');
	% saveas(gcf,['S',num2str(cSub(1)),'~',num2str(cSub(end)),'_',fftType,'_SNRtopo.pdf']);
end

%%%%%%%%%%%%%%  have ploted sorted SNR, the end %%%%%%%%%%%%%%%%%%%%%%%%%%%







%%%%  test_820=[squeeze(mean(MaxChanSNR_SNR(freTimePoint(1),:,:,1),2)) squeeze(mean(MaxChanSNR_SNR(freTimePoint(2),:,:,2),2))];

%%%%% bar(test_820,'grouped','EdgeColor','none');









% %%%% step3
% %---- do fft of baseline and post-baseline([201 700;701 1200; 201:1200]) for 4 conditions,separately ------------------/
% clear all;
% directory     = 'D:\expData\SSVEP_8loc\Exp1_20170328\Data\scanData\20180407\average';
% chanList      = 'all';
% subNum        = [1:21];
% FFTpoints     = 1024;
% windowType    = 'cosin';
% snrskipNum    = 1;
% freRange      = [5 25];
% targetFre     = [8;20]; 
% condition     = {'allcue13_2~45Hz','allcue21_2~45Hz','allcue34_2~45Hz','allcue42_2~45Hz'};
% timeBin       = [-500,-1;201,700;701 1200;201 1200;401 1200];
% allPower_4con = [];
% allSNR_4con   = [];
% ldrTransType  = 'noldr';


% for iCondition = 1:size(condition,2)

% 		cd(fullfile(directory,condition{iCondition}));

% 		for iSub = subNum

% 				iSubFileName = [num2str(iSub),'-',condition{iCondition},'.avg'];

% 				for iTimeRange = 1:size(timeBin,1)  

% 					[powerData,f,chan_names,chanList,rawSignal] = FFTanalysis_bcl(iSubFileName,timeBin(iTimeRange,:),chanList,FFTpoints,windowType,ldrTransType);

% 					allPower_4con_allTimeRange(:,:,iSub,iCondition,iTimeRange)  = powerData; % size: 1-freBins; 2-chan; 3-sub; 4-condition; 5-timeRange(base)

% 					[SNR] = computeSNR(powerData,snrskipNum);

% 					allSNR_4con_allTimeRange(:,:,iSub,iCondition,iTimeRange)    = SNR;

% 					clear SNR powerData;

% 				end
% 		end

% 		cd('..');
% end 


% [freIdx]                   = cutFreBins(f,freRange); 
% allSNR_4con_allTimeRange   = allSNR_4con_allTimeRange(freIdx,:,:,:,:);
% allPower_4con_allTimeRange = allPower_4con_allTimeRange(freIdx,:,:,:,:);
% freBins                    = f(freIdx);

% %-------find closest frePoint-------/
% for iFre = 1:size(targetFre,1)
% 	[X, freIndex]         = sort(abs(freBins-targetFre(iFre,1)),'ascend');
% 	freTimePoint(iFre,1)  = freIndex(1);
% 	realFre(iFre,1) 	  = freBins(freIndex(1));
% end
% %-------have found, the end---------\


% allSNRMinusBase_4con_allTimeRange(:,:,:,:,1) = allSNR_4con_allTimeRange(:,:,:,:,2)-allSNR_4con_allTimeRange(:,:,:,:,1); % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-timeRange(post1,2)
% allSNRMinusBase_4con_allTimeRange(:,:,:,:,2) = allSNR_4con_allTimeRange(:,:,:,:,3)-allSNR_4con_allTimeRange(:,:,:,:,1); % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-timeRange(post1,2)
% allSNRMinusBase_4con_allTimeRange(:,:,:,:,3) = allSNR_4con_allTimeRange(:,:,:,:,4)-allSNR_4con_allTimeRange(:,:,:,:,1); % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-timeRange(post1,2)
% allSNRMinusBase_4con_allTimeRange(:,:,:,:,4) = allSNR_4con_allTimeRange(:,:,:,:,5)-allSNR_4con_allTimeRange(:,:,:,:,1); % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-timeRange(post1,2)

% allFrePointSNRMinusBase_4con_allTimeRange = cat(1,allSNRMinusBase_4con_allTimeRange(freTimePoint(1),:,:,[2 3 4 1],:), allSNRMinusBase_4con_allTimeRange(freTimePoint(2),:,:,[3 2 1 4],:)); % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-timeRange(post1,2)
% allFrePointSNRMinusBase_4con_3Range       = allFrePointSNRMinusBase_4con_allTimeRange(:,:,:,:,[1 2]);
% allFrePointSNRMinusBase_4con_2Range       = allFrePointSNRMinusBase_4con_allTimeRange(:,:,:,:,[3]);
% allFrePointSNRMinusBase_4con_2Range800    = allFrePointSNRMinusBase_4con_allTimeRange(:,:,:,:,[4]);



% allPowerMinusBase_4con_allTimeRange(:,:,:,:,1) = allPower_4con_allTimeRange(:,:,:,:,2)-allPower_4con_allTimeRange(:,:,:,:,1); % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-timeRange(post1,2)
% allPowerMinusBase_4con_allTimeRange(:,:,:,:,2) = allPower_4con_allTimeRange(:,:,:,:,3)-allPower_4con_allTimeRange(:,:,:,:,1); % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-timeRange(post1,2)
% allPowerMinusBase_4con_allTimeRange(:,:,:,:,3) = allPower_4con_allTimeRange(:,:,:,:,4)-allPower_4con_allTimeRange(:,:,:,:,1); % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-timeRange(post1,2)
% allPowerMinusBase_4con_allTimeRange(:,:,:,:,4) = allPower_4con_allTimeRange(:,:,:,:,5)-allPower_4con_allTimeRange(:,:,:,:,1); % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-timeRange(post1,2)

% allFrePointPowerMinusBase_4con_allTimeRange = cat(1,allPowerMinusBase_4con_allTimeRange(freTimePoint(1),:,:,[2 3 4 1],:), allPowerMinusBase_4con_allTimeRange(freTimePoint(2),:,:,[3 2 1 4],:)); % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-timeRange(post1,2)
% allFrePointPowerMinusBase_4con_3Range       = allFrePointPowerMinusBase_4con_allTimeRange(:,:,:,:,[1 2]);
% allFrePointPowerMinusBase_4con_2Range       = allFrePointPowerMinusBase_4con_allTimeRange(:,:,:,:,[3]);
% allFrePointPowerMinusBase_4con_2Range800    = allFrePointPowerMinusBase_4con_allTimeRange(:,:,:,:,[4]);



% save(['fft_4con_allTimeRange_skip' num2str(snrskipNum) ]);
% %---------- It's the end..  Do fft  for 4 conditions,separately. It's the end---------------\





%%%-------------- read eachSub's Max 3 chans' Data ----------------



load('D:\expData\SSVEP_8loc\Exp1_20170328\Data\scanData\20180407\average\spssData_final_20180422\avgFFT_mergeAllconFindChan_-500To1200noldr_skip1.mat')
% load('D:\expData\SSVEP_8loc\Exp1_20170328\Data\scanData\20180407\average\spssData_final_20180422\avgFFT_mergeAllconFindChan_-1000To2000noldr_skip1.mat')
% clearvars -except chooseChanSNR  chooseChanPower  chooseChanSubavgSNR   chooseChanSubavgPower;
load('D:\expData\SSVEP_8loc\Exp1_20170328\Data\scanData\20180407\average\spssData_final_20180422\fft_4con_allTimeRange_skip1.mat')


method =  'eachSubAvg';


%%%%%%%%% 3Range  %%%%%%%%%%
SNR_8hz20hz_3Range        = [];
Power_8hz20hz_3Range      = [];
PowerOfSNR_8hz20hz_3Range = [];

validSubNum =  [1:21]; % [1:21] 1:3,5:7,9:11,13:18,20,21
switch method 

	case 'eachSubAvg'     %%%%%%%%%%   method-1 :according to eachSub's max SNR  %%%%%%%%
		
		for iFre = 1:size(allFrePointSNRMinusBase_4con_3Range,1)  % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-TimeRange

			for iSub = 1:size(allFrePointSNRMinusBase_4con_3Range,3) 
					cSNR = allFrePointSNRMinusBase_4con_3Range(iFre,chooseChanSNR(iSub,:,iFre),iSub,:,:); 
					SNR_8hz20hz_3Range(iFre,:,iSub,:,:) = cSNR;

					cPower = allFrePointPowerMinusBase_4con_3Range(iFre,chooseChanPower(iSub,:,iFre),iSub,:,:); 
					Power_8hz20hz_3Range(iFre,:,iSub,:,:) = cPower;


					ccPower = allFrePointPowerMinusBase_4con_3Range(iFre,chooseChanSNR(iSub,:,iFre),iSub,:,:);

					PowerOfSNR_8hz20hz_3Range(iFre,:,iSub,:,:) = ccPower;

			end
		end


	case 'allsubAvg'   %%%%%%%  method-2 :according to subAvg's max SNR  %%%%%%%
	



		subavgPower = squeeze(mean(frePointPower(:,:,validSubNum),3));
		subavgSNR   = squeeze(mean(frePointSNR(:,:,validSubNum),3));


 

		for iFre = 1:2

			[X, chanIndex]                   = sort(subavgSNR(iFre,[1:size(standard_eName,2)-2]),'descend');
			chooseChanSubavgSNR(iFre,:)      = chanIndex;
			[X, chanIndex]                   = sort(subavgPower(iFre,[1:size(standard_eName,2)-2]),'descend');
			chooseChanSubavgPower(iFre,:)    = chanIndex;

			SNR_8hz20hz_3Range(iFre,:,:,:,:)         = allFrePointSNRMinusBase_4con_3Range(iFre,chooseChanSubavgSNR(iFre,[1:5]),:,:,:);
			PowerOfSNR_8hz20hz_3Range(iFre,:,:,:,:)  = allFrePointPowerMinusBase_4con_3Range(iFre,chooseChanSubavgSNR(iFre,[1:5]),:,:,:);
			Power_8hz20hz_3Range(iFre,:,:,:,:)       = allFrePointPowerMinusBase_4con_3Range(iFre,chooseChanSubavgPower(iFre,[1:5]),:,:,:);

		end
end


SNR_8hz20hz_3Range(3,:,:,:,:)                           = mean(SNR_8hz20hz_3Range,1);
SNR_8hz20hz_3Range(:,6,:,:,:)                           = mean(SNR_8hz20hz_3Range,2);
SNR_8hz20hz_3Range(:,:,numel(validSubNum)+1,:,:)        = mean(SNR_8hz20hz_3Range,3);


Power_8hz20hz_3Range(3,:,:,:,:)                         = mean(Power_8hz20hz_3Range,1);
Power_8hz20hz_3Range(:,6,:,:,:)                         = mean(Power_8hz20hz_3Range,2);
Power_8hz20hz_3Range(:,:,numel(validSubNum)+1,:,:)      = mean(Power_8hz20hz_3Range,3);


PowerOfSNR_8hz20hz_3Range(3,:,:,:,:)                    = mean(PowerOfSNR_8hz20hz_3Range,1);
PowerOfSNR_8hz20hz_3Range(:,6,:,:,:)                    = mean(PowerOfSNR_8hz20hz_3Range,2);
PowerOfSNR_8hz20hz_3Range(:,:,numel(validSubNum)+1,:,:) = mean(PowerOfSNR_8hz20hz_3Range,3);



SNR_8hz20hz_3Range(:,:,:,5,:) = mean(SNR_8hz20hz_3Range(:,:,:,[1 3],:),4);




%%%%%%%%%%%%   spss  %%%%%%%%%%
clear SNR_merge820_3Range pvalue_SNR_merge820_3Range
clear SNR_merge820_2Range pvalue_SNR_merge820_2Range
clear SNR_merge820_2Range800 pvalue_SNR_merge820_2Range800



SNR_merge820_3Range = SNR_8hz20hz_3Range(:,:,validSubNum,:,:);

for iFre = 1:size(SNR_merge820_3Range,1)

	for iTimeRange = 1:size(SNR_merge820_3Range,5)
	 	
	 	for  iChan = 1:size(SNR_merge820_3Range,2)

	 			[H,pvalue_SNR_merge820_3Range(iCmp,iChan,iTimeRange,iFre),CI] = ttest(squeeze(SNR_merge820_3Range(iFre,iChan,:,4,iTimeRange)), squeeze(SNR_merge820_3Range(iFre,iChan,:,5,iTimeRange)));
	 	end
	end
end











SNR_merge820_3Range = SNR_8hz20hz_3Range(:,:,validSubNum,:,:);

for iFre = 1:size(SNR_merge820_3Range,1)

	for iTimeRange = 1:size(SNR_merge820_3Range,5)
	 	
	 	for  iChan = 1:size(SNR_merge820_3Range,2)

	 		for iCmp = 1:3

	 			[H,pvalue_SNR_merge820_3Range(iCmp,iChan,iTimeRange,iFre),CI] = ttest(squeeze(SNR_merge820_3Range(iFre,iChan,:,4,iTimeRange)), squeeze(SNR_merge820_3Range(iFre,iChan,:,iCmp,iTimeRange)));

	 		end
	 	end
	end
end






Power_merge820_3Range = Power_8hz20hz_3Range(:,:,validSubNum,:,:);

for iFre = 1:size(Power_merge820_3Range,1)

	for iTimeRange = 1:size(Power_merge820_3Range,5)
	 	
	 	for  iChan = 1:size(Power_merge820_3Range,2)

	 		for iCmp = 1:3

	 			[H,pvalue_Power_merge820_3Range(iCmp,iChan,iTimeRange,iFre),CI] = ttest(squeeze(Power_merge820_3Range(iFre,iChan,:,4,iTimeRange)), squeeze(Power_merge820_3Range(iFre,iChan,:,iCmp,iTimeRange)));

	 		end
	 	end
	end
end






%%%%%%%%%%%%   spss  %%%%%%%%%%
clear PowerOfSNR_merge820_3Range pvalue_PowerOfSNR_merge820_3Range
clear PowerOfSNR_merge820_2Range pvalue_PowerOfSNR_merge820_2Range
clear PowerOfSNR_merge820_2Range800 pvalue_PowerOfSNR_merge820_2Range800


PowerOfSNR_merge820_3Range = PowerOfSNR_8hz20hz_3Range(:,:,validSubNum,:,:);

for iFre = 1:size(PowerOfSNR_merge820_3Range,1)

	for iTimeRange = 1:size(PowerOfSNR_merge820_3Range,5)
	 	
	 	for  iChan = 1:size(PowerOfSNR_merge820_3Range,2)

	 		for iCmp = 1:3

	 			[H,pvalue_PowerOfSNR_merge820_3Range(iCmp,iChan,iTimeRange,iFre),CI] = ttest(squeeze(PowerOfSNR_merge820_3Range(iFre,iChan,:,4,iTimeRange)), squeeze(PowerOfSNR_merge820_3Range(iFre,iChan,:,iCmp,iTimeRange)));

	 		end
	 	end
	end
end






powerData = [];



for iTimeRange = 1:2

	for iChan = 1:size(Power_merge820_3Range,2)

		cPowerData = squeeze(Power_merge820_3Range(3,iChan,:,:,iTimeRange));
		powerData  = [powerData cPowerData];
	end 
end




%%  for figure of dissertation




load('D:\expData\SSVEP_8loc\Exp1_20170328\Data\scanData\20180407\average\spssData_final_20180422\avgFFT_mergeAllconFindChan_-1000To2000noldr_skip1.mat')
% topo of 8hz SNR [-1000 2000] 21subjects
[myColormap] = makeColormap(20,0,1,'r');
figure;
set(gcf,'Position',[1 1 1600 1080*2/3],'color','w');
for iSub = 1:size(allSNR,3)
	subplot(4,6,iSub);
	topoplot_bcl(allSNR(freTimePoint(1),:,iSub),'EEGChans',standard_eName,'maplimitsDouble',[0 2.5],'colormap',myColormap);
end
set(gcf,'color','w');
print('-depsc','-painters',['D:\expData\SSVEP_8loc\Exp1_20170328\Figure\rawFigure\-1000To2000\SNR_8hz_S21.eps']);


% %% topo of 20hz SNR [-1000 2000] 21subjects
figure;
set(gcf,'Position',[1 1 1600 1080*2/3],'color','w');
for iSub = 1:size(allSNR,3)
	subplot(4,6,iSub);
	topoplot_bcl(allSNR(freTimePoint(2),:,iSub),'EEGChans',standard_eName,'maplimitsDouble',[0 2.5],'colormap',myColormap);
end
set(gcf,'color','w');
print('-depsc','-painters',['D:\expData\SSVEP_8loc\Exp1_20170328\Figure\rawFigure\-1000To2000\SNR_20hz_S21.eps']);


% topo of subAvgSNR 21sajCue_201To700_820_e1
ajUncue_201To700_820_e1
sepCue_201To700_820_e1
sepUncue_201To700_820_e1
ajCue_201To700_820_e2
ajUncue_201To700_820_e2
sepCue_201To700_820_e2
sepUncue_201To700_820_e2ujects [-1000 2000]
figure;
validSubNum  = [1:21];
[myColormap] = makeColormap(20,0,1,'r');
subplot(1,2,1);
topoplot_bcl(squeeze(mean(allSNR(freTimePoint(1),:,validSubNum),3)),'EEGChans',standard_eName,'maplimitsDouble',[0 2],'colormap',myColormap);
subplot(1,2,2);
topoplot_bcl(squeeze(mean(allSNR(freTimePoint(2),:,validSubNum),3)),'EEGChans',standard_eName,'maplimitsDouble',[0 2],'colormap',myColormap);
set(gcf,'color','w');
print('-depsc','-painters',['D:\expData\SSVEP_8loc\Exp1_20170328\Figure\rawFigure\-1000To2000\SNRavg_820_S21.eps']);





% % topo of subAvgPower
% [myColormap] = makeColormap(20,0,1,'b');
% figure;
% topoplot_bcl(squeeze(mean(allPower(freTimePoint(1),:,validSubNum),3)),'EEGChans',standard_eName,'maplimitsDouble',[0 0.15],'colormap',myColormap);
% set(gcf,'color','w');
% print('-depsc','-painters',['D:\expData\SSVEP_8loc\Exp1_20170328\Figure\rawFigure\Poweravg_8hz.eps']);

% figure;
% topoplot_bcl(squeeze(mean(allPower(freTimePoint(2),:,validSubNum),3)),'EEGChans',standard_eName,'maplimitsDouble',[0 0.05],'colormap',myColormap);
% set(gcf,'color','w');
% print('-depsc','-painters',['D:\expData\SSVEP_8loc\Exp1_20170328\Figure\rawFigure\Poweravg_20hz.eps']);

 print('-depsc','-painters',['D:\expData\SSVEP_8loc\Exp1_20170328\Figure\rawFigure\triggerDelay.eps']);





% illustration of computeSNR
[AX,H1,H2] = plotyy([1:13], power./3,[3:11],SNR,'plot');
set(AX(1),'XColor','k','YColor','b','Ylim',[0,0.2],'ytick',[0:0.05:0.2],'Xlim',[0.5 13.5],'Xtick',[0:1:13]);		 				
set(AX(2),'XColor','k','YColor','r','Ylim',[0,4],'ytick',[0:1:4],'Xlim',[0.5 13.5]);
set(H1,'LineStyle','-','color','b');
set(H2,'LineStyle','-','color','r');
print('-depsc','-painters',['D:\expData\SSVEP_8loc\Exp1_20170328\Figure\rawFigure\computeSNR.eps']);




















%%%%%%%%%%%%%  compute errorBar %%%%%%%%%

for iFre = 1:size(SNR_8hz20hz_3Range,1)
	for iChan = 1:size(SNR_8hz20hz_3Range,2)
		for iRange = 1:size(SNR_8hz20hz_3Range,5)
			[gm,gsem] = grpstats(squeeze(SNR_8hz20hz_3Range(iFre,iChan,[1:end-1],:,iRange)),{},{'mean','sem'});
			errorBar_SNR_8hz20hz_3Range(iFre,iChan,:,iRange) = gsem;
			clear gsem;
			[gm,gsem] = grpstats(squeeze(Power_8hz20hz_3Range(iFre,iChan,[1:end-1],:,iRange)),{},{'mean','sem'});
			errorBar_Power_8hz20hz_3Range(iFre,iChan,:,iRange) = gsem;
			clear gsem;
			[gm,gsem] = grpstats(squeeze(PowerOfSNR_8hz20hz_3Range(iFre,iChan,[1:end-1],:,iRange)),{},{'mean','sem'});
			errorBar_PowerOfSNR_8hz20hz_3Range(iFre,iChan,:,iRange) = gsem;
			clear gsem;

		end
	end
end

maxChanNum = 3;
%%-------------- start to plot -------------------/
ngroups    = size(SNR_8hz20hz_3Range,2);
nbars      = size(SNR_8hz20hz_3Range,4);
groupwidth = min(0.8, nbars/(nbars+1.5));	
barColor   = [235 95 52; 254 158 15;63 96 169; 0 0 200]/255;

% % %%%%%%   bar of SNR %%%%%%%%
% pptTitleName = ['SNR_max', num2str(maxChanNum), 'chan_3Range_skip' num2str(snrskipNum) ];
% mkdir(pptTitleName);
% [opt,ppt] = figPPT(pptTitleName);

% plotSNRRange = [-3 3; -1  1];
% freStr = {'8hz','20hz','mean of 8+20'};
% timeBinStr = {'201To700','701To1200','mean of 2Range'};

% for iSub = 1:size(allFrePointSNRMinusBase_4con_3Range,3)+1
				
% 	figure;
%     set (gcf,'Position',[50,50,1800,900], 'color','w');

% 	for iRange = 1:size(SNR_8hz20hz_3Range,5)		

% 		for iFre = 1:size(SNR_8hz20hz_3Range,1)  % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-TimeRange
		
% 			subplot(3,3,3*(iRange-1)+iFre);

% 			cSNR = squeeze(SNR_8hz20hz_3Range(iFre,:,iSub,:,iRange)); 
% 			% cError = SNR_8hz20hz_3Range(iFre,:,iSub+1,:,iRange);			
% 			H=bar(cSNR,'grouped','EdgeColor','none');
% 		    ch = get(H,'children');
% 	        set(ch{1},'Facecolor',barColor(1,:));
% 	        set(ch{2},'Facecolor',barColor(2,:));
% 	        set(ch{3},'Facecolor',barColor(3,:));
% 	        set(ch{4},'Facecolor',barColor(4,:));
% 			hold on;box off;
% 			ylim(plotSNRRange(1,:));

% 			if iSub==size(allFrePointSNRMinusBase_4con_3Range,3)+1
% 				subStr = 'allSubMean';
% 				cError = squeeze(errorBar_SNR_8hz20hz_3Range(iFre,:,:,iRange));
% 				for i = 1:nbars
% 				    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
% 				    errorbar(x,cSNR(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
% 				    hold on;
% 					ylim(plotSNRRange(2,:));
% 				end

% 			elseif iSub<=size(allFrePointSNRMinusBase_4con_3Range,3)
% 					subStr = ['Sub:',num2str(iSub)];
% 			end

% 			title([subStr,freStr(iFre),timeBinStr(iRange)]);
% 			clear cSNR cError;
% 		end
% 	end
% 	l = legend([ch{1} ch{2} ch{3} ch{4} ],'sepCue--Uncue','ajCue--Uncue','ajCue--Cue','SepCue--Cue','Location',[0.908333333333333 0.248888888888889 0.0816666666666667 0.0881481481481481]);
% 	legend('boxoff');
% 	cfigFileName = fullfile(fullfile(pwd,pptTitleName),[pptTitleName,'_',subStr([1:3,5:end]),'.bmp']);
% 	saveas(gcf,cfigFileName);
% 	ppt = ppt.addImageSlide(['SNR_',subStr(1:3),'_',subStr(5:end)],cfigFileName,opt);
% 	close all;
% end

% ppt.saveAs(fullfile(pwd,pptTitleName));

% ppt.close;




% % %%%%%%   bar of Power %%%%%%%%
% pptTitleName = ['Power_max', num2str(maxChanNum), 'chan_3Range_skip' num2str(snrskipNum) ];
% mkdir(pptTitleName);
% [opt,ppt] = figPPT(pptTitleName);

% plotPowerRange = [-0.3 0.3; -0.1  0.1];
% freStr = {'8hz','20hz','mean of 8+20'};
% timeBinStr = {'201To700','701To1200','mean of 2Range'};

% for iSub = 1:size(allFrePointPowerMinusBase_4con_3Range,3)+1
				
% 	figure;
%     set (gcf,'Position',[50,50,1800,900], 'color','w');

% 	for iRange = 1:size(Power_8hz20hz_3Range,5)		

% 		for iFre = 1:size(Power_8hz20hz_3Range,1)  % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-TimeRange
		
% 			subplot(3,3,3*(iRange-1)+iFre);

% 			cPower = squeeze(Power_8hz20hz_3Range(iFre,:,iSub,:,iRange)); 
% 			% cError = Power_8hz20hz_3Range(iFre,:,iSub+1,:,iRange);			
% 			H=bar(cPower,'grouped','EdgeColor','none');
% 		    ch = get(H,'children');
% 	        set(ch{1},'Facecolor',barColor(1,:));
% 	        set(ch{2},'Facecolor',barColor(2,:));
% 	        set(ch{3},'Facecolor',barColor(3,:));
% 	        set(ch{4},'Facecolor',barColor(4,:));
% 			hold on;box off;
% 			ylim(plotPowerRange(1,:));

% 			if iSub==size(allFrePointPowerMinusBase_4con_3Range,3)+1
% 				subStr = 'allSubMean';
% 				cError = squeeze(errorBar_Power_8hz20hz_3Range(iFre,:,:,iRange));
% 				for i = 1:nbars
% 				    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
% 				    errorbar(x,cPower(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
% 				    hold on;
% 					ylim(plotPowerRange(2,:));
% 				end

% 			elseif iSub<=size(allFrePointPowerMinusBase_4con_3Range,3)
% 					subStr = ['Sub:',num2str(iSub)];
% 			end

% 			title([subStr,freStr(iFre),timeBinStr(iRange)]);
% 			clear cPower cError;
% 		end
% 	end
% 	l = legend([ch{1} ch{2} ch{3} ch{4} ],'sepCue--Uncue','ajCue--Uncue','ajCue--Cue','SepCue--Cue','Location',[0.908333333333333 0.248888888888889 0.0816666666666667 0.0881481481481481]);
% 	legend('boxoff');
% 	cfigFileName = fullfile(fullfile(pwd,pptTitleName),[pptTitleName,'_',subStr([1:3,5:end]),'.bmp']);
% 	saveas(gcf,cfigFileName);
% 	ppt = ppt.addImageSlide(['Power_',subStr(1:3),'_',subStr(5:end)],cfigFileName,opt);
% 	close all;
% end

% ppt.saveAs(fullfile(pwd,pptTitleName));

% ppt.close;






% % %%%%%%   bar of PowerOfSNR %%%%%%%%
% pptTitleName = ['PowerOfSNR_max', num2str(maxChanNum), 'chan_3Range_skip' num2str(snrskipNum) ];
% mkdir(pptTitleName);
% [opt,ppt] = figPPT(pptTitleName);

% plotPowerOfSNRRange = [-0.3 0.3; -0.1  0.1];
% freStr = {'8hz','20hz','mean of 8+20'};
% timeBinStr = {'201To700','701To1200','mean of 2Range'};

% for iSub = 1:size(allFrePointPowerMinusBase_4con_3Range,3)+1
				
% 	figure;
%     set (gcf,'Position',[50,50,1800,900], 'color','w');

% 	for iRange = 1:size(PowerOfSNR_8hz20hz_3Range,5)		

% 		for iFre = 1:size(PowerOfSNR_8hz20hz_3Range,1)  % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-TimeRange
		
% 			subplot(3,3,3*(iRange-1)+iFre);

% 			cPowerOfSNR = squeeze(PowerOfSNR_8hz20hz_3Range(iFre,:,iSub,:,iRange)); 
% 			% cError = PowerOfSNR_8hz20hz_3Range(iFre,:,iSub+1,:,iRange);			
% 			H=bar(cPowerOfSNR,'grouped','EdgeColor','none');
% 		    ch = get(H,'children');
% 	        set(ch{1},'Facecolor',barColor(1,:));
% 	        set(ch{2},'Facecolor',barColor(2,:));
% 	        set(ch{3},'Facecolor',barColor(3,:));
% 	        set(ch{4},'Facecolor',barColor(4,:));
% 			hold on;box off;
% 			ylim(plotPowerOfSNRRange(1,:));

% 			if iSub==size(allFrePointPowerMinusBase_4con_3Range,3)+1
% 				subStr = 'allSubMean';
% 				cError = squeeze(errorBar_PowerOfSNR_8hz20hz_3Range(iFre,:,:,iRange));
% 				for i = 1:nbars
% 				    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
% 				    errorbar(x,cPowerOfSNR(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
% 				    hold on;
% 					ylim(plotPowerOfSNRRange(2,:));
% 				end

% 			elseif iSub<=size(allFrePointPowerMinusBase_4con_3Range,3)
% 					subStr = ['Sub:',num2str(iSub)];
% 			end

% 			title([subStr,freStr(iFre),timeBinStr(iRange)]);
% 			clear cPowerOfSNR cError;
% 		end
% 	end
% 	l = legend([ch{1} ch{2} ch{3} ch{4} ],'sepCue--Uncue','ajCue--Uncue','ajCue--Cue','SepCue--Cue','Location',[0.908333333333333 0.248888888888889 0.0816666666666667 0.0881481481481481]);
% 	legend('boxoff');
% 	cfigFileName = fullfile(fullfile(pwd,pptTitleName),[pptTitleName,'_',subStr([1:3,5:end]),'.bmp']);
% 	saveas(gcf,cfigFileName);
% 	ppt = ppt.addImageSlide(['PowerOfSNR_',subStr(1:3),'_',subStr(5:end)],cfigFileName,opt);
% 	close all;
% end

% ppt.saveAs(fullfile(pwd,pptTitleName));

% ppt.close;





%% method-2 :according to subAvg's max SNR




%%%%% 2Range

%%%%%%%%% 2Range  %%%%%%%%%%
SNR_8hz20hz_2Range        = [];
Power_8hz20hz_2Range      = [];
PowerOfSNR_8hz20hz_2Range = [];

for iFre = 1:size(allFrePointSNRMinusBase_4con_2Range,1)  % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-TimeRange

	for iSub = 1:size(allFrePointSNRMinusBase_4con_2Range,3) 
			cSNR = allFrePointSNRMinusBase_4con_2Range(iFre,chooseChanSNR(iSub,:,iFre),iSub,:,:); 
			SNR_8hz20hz_2Range(iFre,:,iSub,:,:) = cSNR;

			cPower = allFrePointPowerMinusBase_4con_2Range(iFre,chooseChanPower(iSub,:,iFre),iSub,:,:); 
			Power_8hz20hz_2Range(iFre,:,iSub,:,:) = cPower;

			ccPower = allFrePointPowerMinusBase_4con_2Range(iFre,chooseChanSNR(iSub,:,iFre),iSub,:,:); 
			PowerOfSNR_8hz20hz_2Range(iFre,:,iSub,:,:) = ccPower;

	end
end

SNR_8hz20hz_2Range(3,:,:,:,:)                                = mean(SNR_8hz20hz_2Range,1);
SNR_8hz20hz_2Range(:,size(chooseChanSNR,2)+1,:,:,:)          = mean(SNR_8hz20hz_2Range,2);
SNR_8hz20hz_2Range(:,:,size(chooseChanSNR,1)+1,:,:)          = mean(SNR_8hz20hz_2Range,3);

Power_8hz20hz_2Range(3,:,:,:,:)                              = mean(Power_8hz20hz_2Range,1);
Power_8hz20hz_2Range(:,size(chooseChanPower,2)+1,:,:,:)      = mean(Power_8hz20hz_2Range,2);
Power_8hz20hz_2Range(:,:,size(chooseChanPower,1)+1,:,:)      = mean(Power_8hz20hz_2Range,3);

PowerOfSNR_8hz20hz_2Range(3,:,:,:,:)                         = mean(PowerOfSNR_8hz20hz_2Range,1);
PowerOfSNR_8hz20hz_2Range(:,size(chooseChanPower,2)+1,:,:,:) = mean(PowerOfSNR_8hz20hz_2Range,2);
PowerOfSNR_8hz20hz_2Range(:,:,size(chooseChanPower,1)+1,:,:) = mean(PowerOfSNR_8hz20hz_2Range,3);



%%%%%%%%%%%%%  compute errorBar %%%%%%%%%

for iFre = 1:size(SNR_8hz20hz_2Range,1)
	for iChan = 1:size(SNR_8hz20hz_2Range,2)
		for iRange = 1:size(SNR_8hz20hz_2Range,5)
			[gm,gsem] = grpstats(squeeze(SNR_8hz20hz_2Range(iFre,iChan,[1:end-1],:,iRange)),{},{'mean','sem'});
			errorBar_SNR_8hz20hz_2Range(iFre,iChan,:,iRange) = gsem;
			clear gsem;
			[gm,gsem] = grpstats(squeeze(Power_8hz20hz_2Range(iFre,iChan,[1:end-1],:,iRange)),{},{'mean','sem'});
			errorBar_Power_8hz20hz_2Range(iFre,iChan,:,iRange) = gsem;
			clear gsem;
			[gm,gsem] = grpstats(squeeze(PowerOfSNR_8hz20hz_2Range(iFre,iChan,[1:end-1],:,iRange)),{},{'mean','sem'});
			errorBar_PowerOfSNR_8hz20hz_2Range(iFre,iChan,:,iRange) = gsem;
			clear gsem;
		end
	end
end

maxChanNum = 5;
%%-------------- start to plot -------------------/
ngroups    = size(SNR_8hz20hz_2Range,2);
nbars      = size(SNR_8hz20hz_2Range,4);
groupwidth = min(0.8, nbars/(nbars+1.5));	
barColor   = [235 95 52; 254 158 15;63 96 169; 0 0 200]/255;

% %%%%%%   bar of SNR %%%%%%%%
% pptTitleName = ['SNR_max', num2str(maxChanNum), 'chan_2Range_skip' num2str(snrskipNum)];
% mkdir(pptTitleName);
% [opt,ppt] = figPPT(pptTitleName);

% plotSNRRange = [-3 3; -1  1];
% freStr = {'8hz','20hz','mean of 8+20'};
% timeBinStr = {'201To1200'};

% for iSub = 1:size(allFrePointSNRMinusBase_4con_2Range,3)+1
				
% 	figure;
%     set (gcf,'Position',[50,50,1800,900], 'color','w');

% 	for iRange = 1:size(SNR_8hz20hz_2Range,5)		

% 		for iFre = 1:size(SNR_8hz20hz_2Range,1)  % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-TimeRange
		
% 			subplot(3,3,3*(iRange-1)+iFre);

% 			cSNR = squeeze(SNR_8hz20hz_2Range(iFre,:,iSub,:,iRange)); 
% 			% cError = SNR_8hz20hz_2Range(iFre,:,iSub+1,:,iRange);			
% 			H=bar(cSNR,'grouped','EdgeColor','none');
% 		    ch = get(H,'children');
% 	        set(ch{1},'Facecolor',barColor(1,:));
% 	        set(ch{2},'Facecolor',barColor(2,:));
% 	        set(ch{3},'Facecolor',barColor(3,:));
% 	        set(ch{4},'Facecolor',barColor(4,:));
% 			hold on;box off;
% 			ylim(plotSNRRange(1,:));

% 			if iSub==size(allFrePointSNRMinusBase_4con_2Range,3)+1
% 				subStr = 'allSubMean';
% 				cError = squeeze(errorBar_SNR_8hz20hz_2Range(iFre,:,:,iRange));
% 				for i = 1:nbars
% 				    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
% 				    errorbar(x,cSNR(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
% 				    hold on;
% 					ylim(plotSNRRange(2,:));
% 				end

% 			elseif iSub<=size(allFrePointSNRMinusBase_4con_2Range,3)
% 					subStr = ['Sub:',num2str(iSub)];
% 			end

% 			title([subStr,freStr(iFre),timeBinStr(iRange)]);
% 			clear cSNR cError;
% 		end
% 	end
% 	l = legend([ch{1} ch{2} ch{3} ch{4} ],'sepCue--Uncue','ajCue--Uncue','ajCue--Cue','SepCue--Cue','Location',[0.908333333333333 0.248888888888889 0.0816666666666667 0.0881481481481481]);
% 	legend('boxoff');
% 	cfigFileName = fullfile(fullfile(pwd,pptTitleName),[pptTitleName,'_',subStr([1:3,5:end]),'.bmp']);
% 	saveas(gcf,cfigFileName);
% 	ppt = ppt.addImageSlide(['SNR_',subStr(1:3),'_',subStr(5:end)],cfigFileName,opt);
% 	close all;
% end

% ppt.saveAs(fullfile(pwd,pptTitleName));

% ppt.close;




% %%%%%%   bar of Power %%%%%%%%
% pptTitleName = ['Power_max', num2str(maxChanNum), 'chan_2Range_skip' num2str(snrskipNum)];
% mkdir(pptTitleName);
% [opt,ppt] = figPPT(pptTitleName);

% plotPowerRange = [-0.3 0.3; -0.1  0.1];
% freStr = {'8hz','20hz','mean of 8+20'};
% timeBinStr = {'201To1200'};

% for iSub = 1: size(allFrePointPowerMinusBase_4con_2Range,3)+1
				
% 	figure;
%     set (gcf,'Position',[50,50,1800,900], 'color','w');

% 	for iRange = 1:size(Power_8hz20hz_2Range,5)		

% 		for iFre = 1:size(Power_8hz20hz_2Range,1)  % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-TimeRange
		
% 			subplot(3,3,3*(iRange-1)+iFre);

% 			cPower = squeeze(Power_8hz20hz_2Range(iFre,:,iSub,:,iRange)); 
% 			% cError = Power_8hz20hz_2Range(iFre,:,iSub+1,:,iRange);			
% 			H=bar(cPower,'grouped','EdgeColor','none');
% 		    ch = get(H,'children');
% 	        set(ch{1},'Facecolor',barColor(1,:));
% 	        set(ch{2},'Facecolor',barColor(2,:));
% 	        set(ch{3},'Facecolor',barColor(3,:));
% 	        set(ch{4},'Facecolor',barColor(4,:));
% 			hold on;box off;
% 			ylim(plotPowerRange(1,:));

% 			if iSub==size(allFrePointPowerMinusBase_4con_2Range,3)+1
% 				subStr = 'allSubMean';
% 				cError = squeeze(errorBar_Power_8hz20hz_2Range(iFre,:,:,iRange));
% 				for i = 1:nbars
% 				    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
% 				    errorbar(x,cPower(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
% 				    hold on;
% 					ylim(plotPowerRange(2,:));
% 				end

% 			elseif iSub<=size(allFrePointPowerMinusBase_4con_2Range,3)
% 					subStr = ['Sub:',num2str(iSub)];
% 			end

% 			title([subStr,freStr(iFre),timeBinStr(iRange)]);
% 			clear cPower cError;
% 		end
% 	end
% 	l = legend([ch{1} ch{2} ch{3} ch{4} ],'sepCue--Uncue','ajCue--Uncue','ajCue--Cue','SepCue--Cue','Location',[0.908333333333333 0.248888888888889 0.0816666666666667 0.0881481481481481]);
% 	legend('boxoff');
% 	cfigFileName = fullfile(fullfile(pwd,pptTitleName),[pptTitleName,'_',subStr([1:3,5:end]),'.bmp']);
% 	saveas(gcf,cfigFileName);
% 	ppt = ppt.addImageSlide(['Power_',subStr(1:3),'_',subStr(5:end)],cfigFileName,opt);
% 	close all;
% end

% ppt.saveAs(fullfile(pwd,pptTitleName));

% ppt.close;





% %%%%%%   bar of PowerOfSNR %%%%%%%%
pptTitleName = ['PowerOfSNR_max', num2str(maxChanNum), 'chan_2Range_skip' num2str(snrskipNum) ];
mkdir(pptTitleName);
[opt,ppt] = figPPT(pptTitleName);

plotPowerOfSNRRange = [-0.3 0.3; -0.1  0.1];
freStr = {'8hz','20hz','mean of 8+20'};
timeBinStr = {'201To1200'};

for iSub = 1:size(allFrePointPowerMinusBase_4con_2Range,3)+1
				
	figure;
    set (gcf,'Position',[50,50,1800,900], 'color','w');

	for iRange = 1:size(PowerOfSNR_8hz20hz_2Range,5)		

		for iFre = 1:size(PowerOfSNR_8hz20hz_2Range,1)  % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-TimeRange
		
			subplot(3,3,3*(iRange-1)+iFre);

			cPowerOfSNR = squeeze(PowerOfSNR_8hz20hz_2Range(iFre,:,iSub,:,iRange)); 
			% cError = PowerOfSNR_8hz20hz_3Range(iFre,:,iSub+1,:,iRange);			
			H=bar(cPowerOfSNR,'grouped','EdgeColor','none');
		    ch = get(H,'children');
	        set(ch{1},'Facecolor',barColor(1,:));
	        set(ch{2},'Facecolor',barColor(2,:));
	        set(ch{3},'Facecolor',barColor(3,:));
	        set(ch{4},'Facecolor',barColor(4,:));
			hold on;box off;
			ylim(plotPowerOfSNRRange(1,:));

			if iSub==size(allFrePointPowerMinusBase_4con_2Range,3)+1
				subStr = 'allSubMean';
				cError = squeeze(errorBar_PowerOfSNR_8hz20hz_2Range(iFre,:,:,iRange));
				for i = 1:nbars
				    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
				    errorbar(x,cPowerOfSNR(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
				    hold on;
					ylim(plotPowerOfSNRRange(2,:));
				end

			elseif iSub<=size(allFrePointPowerMinusBase_4con_2Range,3)
					subStr = ['Sub:',num2str(iSub)];
			end

			title([subStr,freStr(iFre),timeBinStr(iRange)]);
			clear cPowerOfSNR cError;
		end
	end
	l = legend([ch{1} ch{2} ch{3} ch{4} ],'sepCue--Uncue','ajCue--Uncue','ajCue--Cue','SepCue--Cue','Location',[0.908333333333333 0.248888888888889 0.0816666666666667 0.0881481481481481]);
	legend('boxoff');
	cfigFileName = fullfile(fullfile(pwd,pptTitleName),[pptTitleName,'_',subStr([1:3,5:end]),'.bmp']);
	saveas(gcf,cfigFileName);
	ppt = ppt.addImageSlide(['PowerOfSNR_',subStr(1:3),'_',subStr(5:end)],cfigFileName,opt);
	close all;
end

ppt.saveAs(fullfile(pwd,pptTitleName));

ppt.close;
















%%%%% 2Range800

%%%%%%%%% 2Range800  %%%%%%%%%%
SNR_8hz20hz_2Range800        = [];
Power_8hz20hz_2Range800      = [];
PowerOfSNR_8hz20hz_2Range800 = [];

for iFre = 1:size(allFrePointSNRMinusBase_4con_2Range800,1)  % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-TimeRange

	for iSub = 1:size(allFrePointSNRMinusBase_4con_2Range800,3) 
			cSNR = allFrePointSNRMinusBase_4con_2Range800(iFre,chooseChanSNR(iSub,:,iFre),iSub,:,:); 
			SNR_8hz20hz_2Range800(iFre,:,iSub,:,:) = cSNR;

			cPower = allFrePointPowerMinusBase_4con_2Range800(iFre,chooseChanPower(iSub,:,iFre),iSub,:,:); 
			Power_8hz20hz_2Range800(iFre,:,iSub,:,:) = cPower;

			ccPower = allFrePointPowerMinusBase_4con_2Range800(iFre,chooseChanSNR(iSub,:,iFre),iSub,:,:); 
			PowerOfSNR_8hz20hz_2Range800(iFre,:,iSub,:,:) = ccPower;

	end
end

SNR_8hz20hz_2Range800(3,:,:,:,:)                                = mean(SNR_8hz20hz_2Range800,1);
SNR_8hz20hz_2Range800(:,size(chooseChanSNR,2)+1,:,:,:)          = mean(SNR_8hz20hz_2Range800,2);
SNR_8hz20hz_2Range800(:,:,size(chooseChanSNR,1)+1,:,:)          = mean(SNR_8hz20hz_2Range800,3);

Power_8hz20hz_2Range800(3,:,:,:,:)                              = mean(Power_8hz20hz_2Range800,1);
Power_8hz20hz_2Range800(:,size(chooseChanPower,2)+1,:,:,:)      = mean(Power_8hz20hz_2Range800,2);
Power_8hz20hz_2Range800(:,:,size(chooseChanPower,1)+1,:,:)      = mean(Power_8hz20hz_2Range800,3);

PowerOfSNR_8hz20hz_2Range800(3,:,:,:,:)                         = mean(PowerOfSNR_8hz20hz_2Range800,1);
PowerOfSNR_8hz20hz_2Range800(:,size(chooseChanPower,2)+1,:,:,:) = mean(PowerOfSNR_8hz20hz_2Range800,2);
PowerOfSNR_8hz20hz_2Range800(:,:,size(chooseChanPower,1)+1,:,:) = mean(PowerOfSNR_8hz20hz_2Range800,3);



%%%%%%%%%%%%%  compute errorBar %%%%%%%%%

for iFre = 1:size(SNR_8hz20hz_2Range800,1)
	for iChan = 1:size(SNR_8hz20hz_2Range800,2)
		for iRange = 1:size(SNR_8hz20hz_2Range800,5)
			[gm,gsem] = grpstats(squeeze(SNR_8hz20hz_2Range800(iFre,iChan,[1:end-1],:,iRange)),{},{'mean','sem'});
			errorBar_SNR_8hz20hz_2Range800(iFre,iChan,:,iRange) = gsem;
			clear gsem;
			[gm,gsem] = grpstats(squeeze(Power_8hz20hz_2Range800(iFre,iChan,[1:end-1],:,iRange)),{},{'mean','sem'});
			errorBar_Power_8hz20hz_2Range800(iFre,iChan,:,iRange) = gsem;
			clear gsem;
			[gm,gsem] = grpstats(squeeze(PowerOfSNR_8hz20hz_2Range800(iFre,iChan,[1:end-1],:,iRange)),{},{'mean','sem'});
			errorBar_PowerOfSNR_8hz20hz_2Range800(iFre,iChan,:,iRange) = gsem;
			clear gsem;
		end
	end
end

maxChanNum = 5;
%%-------------- start to plot -------------------/
ngroups    = size(SNR_8hz20hz_2Range800,2);
nbars      = size(SNR_8hz20hz_2Range800,4);
groupwidth = min(0.8, nbars/(nbars+1.5));	
barColor   = [235 95 52; 254 158 15;63 96 169; 0 0 200]/255;

% %%%%%%   bar of SNR %%%%%%%%
% pptTitleName = ['SNR_max', num2str(maxChanNum), 'chan_2Range800_skip' num2str(snrskipNum)];
% mkdir(pptTitleName);
% [opt,ppt] = figPPT(pptTitleName);

% plotSNRRange = [-3 3; -1  1];
% freStr = {'8hz','20hz','mean of 8+20'};
% timeBinStr = {'401To1200'};

% for iSub = 1:size(allFrePointSNRMinusBase_4con_2Range800,3)+1
				
% 	figure;
%     set (gcf,'Position',[50,50,1800,900], 'color','w');

% 	for iRange = 1:size(SNR_8hz20hz_2Range800,5)		

% 		for iFre = 1:size(SNR_8hz20hz_2Range800,1)  % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-TimeRange
		
% 			subplot(3,3,3*(iRange-1)+iFre);

% 			cSNR = squeeze(SNR_8hz20hz_2Range800(iFre,:,iSub,:,iRange)); 
% 			% cError = SNR_8hz20hz_2Range800(iFre,:,iSub+1,:,iRange);			
% 			H=bar(cSNR,'grouped','EdgeColor','none');
% 		    ch = get(H,'children');
% 	        set(ch{1},'Facecolor',barColor(1,:));
% 	        set(ch{2},'Facecolor',barColor(2,:));
% 	        set(ch{3},'Facecolor',barColor(3,:));
% 	        set(ch{4},'Facecolor',barColor(4,:));
% 			hold on;box off;
% 			ylim(plotSNRRange(1,:));

% 			if iSub==size(allFrePointSNRMinusBase_4con_2Range800,3)+1
% 				subStr = 'allSubMean';
% 				cError = squeeze(errorBar_SNR_8hz20hz_2Range800(iFre,:,:,iRange));
% 				for i = 1:nbars
% 				    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
% 				    errorbar(x,cSNR(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
% 				    hold on;
% 					ylim(plotSNRRange(2,:));
% 				end

% 			elseif iSub<=size(allFrePointSNRMinusBase_4con_2Range800,3)
% 					subStr = ['Sub:',num2str(iSub)];
% 			end

% 			title([subStr,freStr(iFre),timeBinStr(iRange)]);
% 			clear cSNR cError;
% 		end
% 	end
% 	l = legend([ch{1} ch{2} ch{3} ch{4} ],'sepCue--Uncue','ajCue--Uncue','ajCue--Cue','SepCue--Cue','Location',[0.908333333333333 0.248888888888889 0.0816666666666667 0.0881481481481481]);
% 	legend('boxoff');
% 	cfigFileName = fullfile(fullfile(pwd,pptTitleName),[pptTitleName,'_',subStr([1:3,5:end]),'.bmp']);
% 	saveas(gcf,cfigFileName);
% 	ppt = ppt.addImageSlide(['SNR_',subStr(1:3),'_',subStr(5:end)],cfigFileName,opt);
% 	close all;
% end

% ppt.saveAs(fullfile(pwd,pptTitleName));

% ppt.close;




% %%%%%%   bar of Power %%%%%%%%
% pptTitleName = ['Power_max', num2str(maxChanNum), 'chan_2Range800_skip' num2str(snrskipNum)];
% mkdir(pptTitleName);
% [opt,ppt] = figPPT(pptTitleName);

% plotPowerRange = [-0.3 0.3; -0.1  0.1];
% freStr = {'8hz','20hz','mean of 8+20'};
% timeBinStr = {'401To1200'};

% for iSub = 1: size(allFrePointPowerMinusBase_4con_2Range800,3)+1
				
% 	figure;
%     set (gcf,'Position',[50,50,1800,900], 'color','w');

% 	for iRange = 1:size(Power_8hz20hz_2Range800,5)		

% 		for iFre = 1:size(Power_8hz20hz_2Range800,1)  % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-TimeRange
		
% 			subplot(3,3,3*(iRange-1)+iFre);

% 			cPower = squeeze(Power_8hz20hz_2Range800(iFre,:,iSub,:,iRange)); 
% 			% cError = Power_8hz20hz_2Range800(iFre,:,iSub+1,:,iRange);			
% 			H=bar(cPower,'grouped','EdgeColor','none');
% 		    ch = get(H,'children');
% 	        set(ch{1},'Facecolor',barColor(1,:));
% 	        set(ch{2},'Facecolor',barColor(2,:));
% 	        set(ch{3},'Facecolor',barColor(3,:));
% 	        set(ch{4},'Facecolor',barColor(4,:));
% 			hold on;box off;
% 			ylim(plotPowerRange(1,:));

% 			if iSub==size(allFrePointPowerMinusBase_4con_2Range800,3)+1
% 				subStr = 'allSubMean';
% 				cError = squeeze(errorBar_Power_8hz20hz_2Range800(iFre,:,:,iRange));
% 				for i = 1:nbars
% 				    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
% 				    errorbar(x,cPower(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
% 				    hold on;
% 					ylim(plotPowerRange(2,:));
% 				end

% 			elseif iSub<=size(allFrePointPowerMinusBase_4con_2Range800,3)
% 					subStr = ['Sub:',num2str(iSub)];
% 			end

% 			title([subStr,freStr(iFre),timeBinStr(iRange)]);
% 			clear cPower cError;
% 		end
% 	end
% 	l = legend([ch{1} ch{2} ch{3} ch{4} ],'sepCue--Uncue','ajCue--Uncue','ajCue--Cue','SepCue--Cue','Location',[0.908333333333333 0.248888888888889 0.0816666666666667 0.0881481481481481]);
% 	legend('boxoff');
% 	cfigFileName = fullfile(fullfile(pwd,pptTitleName),[pptTitleName,'_',subStr([1:3,5:end]),'.bmp']);
% 	saveas(gcf,cfigFileName);
% 	ppt = ppt.addImageSlide(['Power_',subStr(1:3),'_',subStr(5:end)],cfigFileName,opt);
% 	close all;
% end

% ppt.saveAs(fullfile(pwd,pptTitleName));

% ppt.close;





% %%%%%%   bar of PowerOfSNR %%%%%%%%
pptTitleName = ['PowerOfSNR_max', num2str(maxChanNum), 'chan_2Range800_skip' num2str(snrskipNum) ];
mkdir(pptTitleName);
[opt,ppt] = figPPT(pptTitleName);

plotPowerOfSNRRange = [-0.3 0.3; -0.1  0.1];
freStr = {'8hz','20hz','mean of 8+20'};
timeBinStr = {'401To1200'};

for iSub = 1:size(allFrePointPowerMinusBase_4con_2Range800,3)+1
				
	figure;
    set (gcf,'Position',[50,50,1800,900], 'color','w');

	for iRange = 1:size(PowerOfSNR_8hz20hz_2Range800,5)		

		for iFre = 1:size(PowerOfSNR_8hz20hz_2Range800,1)  % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-TimeRange
		
			subplot(3,3,3*(iRange-1)+iFre);

			cPowerOfSNR = squeeze(PowerOfSNR_8hz20hz_2Range800(iFre,:,iSub,:,iRange)); 
			% cError = PowerOfSNR_8hz20hz_3Range(iFre,:,iSub+1,:,iRange);			
			H=bar(cPowerOfSNR,'grouped','EdgeColor','none');
		    ch = get(H,'children');
	        set(ch{1},'Facecolor',barColor(1,:));
	        set(ch{2},'Facecolor',barColor(2,:));
	        set(ch{3},'Facecolor',barColor(3,:));
	        set(ch{4},'Facecolor',barColor(4,:));
			hold on;box off;
			ylim(plotPowerOfSNRRange(1,:));

			if iSub==size(allFrePointPowerMinusBase_4con_2Range800,3)+1
				subStr = 'allSubMean';
				cError = squeeze(errorBar_PowerOfSNR_8hz20hz_2Range800(iFre,:,:,iRange));
				for i = 1:nbars
				    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
				    errorbar(x,cPowerOfSNR(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
				    hold on;
					ylim(plotPowerOfSNRRange(2,:));
				end

			elseif iSub<=size(allFrePointPowerMinusBase_4con_2Range800,3)
					subStr = ['Sub:',num2str(iSub)];
			end

			title([subStr,freStr(iFre),timeBinStr(iRange)]);
			clear cPowerOfSNR cError;
		end
	end
	l = legend([ch{1} ch{2} ch{3} ch{4} ],'sepCue--Uncue','ajCue--Uncue','ajCue--Cue','SepCue--Cue','Location',[0.908333333333333 0.248888888888889 0.0816666666666667 0.0881481481481481]);
	legend('boxoff');
	cfigFileName = fullfile(fullfile(pwd,pptTitleName),[pptTitleName,'_',subStr([1:3,5:end]),'.bmp']);
	saveas(gcf,cfigFileName);
	ppt = ppt.addImageSlide(['PowerOfSNR_',subStr(1:3),'_',subStr(5:end)],cfigFileName,opt);
	close all;
end

ppt.saveAs(fullfile(pwd,pptTitleName));

ppt.close;
















%%%%%%%%%%%%   spss  %%%%%%%%%%
% PowerForSpss_eachSubMax_3Range = [];

% for iFre = 1:size(Power_8hz20hz_3Range,1)
% 	for iChan = 1:size(Power_8hz20hz_3Range,2)
% 		cData = squeeze(Power_8hz20hz_3Range(iFre,iChan,[1:size(Power_8hz20hz_3Range,3)-1],:));
% 		PowerForSpss_eachSubMax_3Range = [PowerForSpss_eachSubMax_3Range cData];
% 	end
% end



% clear SNR_merge820_3Range pvalue_SNR_merge820_3Range
% clear SNR_merge820_2Range pvalue_SNR_merge820_2Range
% clear SNR_merge820_2Range800 pvalue_SNR_merge820_2Range800


% SNR_merge820_3Range = SNR_8hz20hz_3Range(3,:,[1:21],:,:);

% for iTimeRange = 1:size(SNR_merge820_3Range,5)
 	
%  	for  iChan = 1:size(SNR_merge820_3Range,2)

%  		for iCmp = 1:3

%  			[H,pvalue_SNR_merge820_3Range(iCmp,iChan,iTimeRange),CI] = ttest(squeeze(SNR_merge820_3Range(1,iChan,:,1,iTimeRange)), squeeze(SNR_merge820_3Range(1,iChan,:,iCmp+1,iTimeRange)));

%  		end
%  	end
% end




% SNR_merge820_2Range = SNR_8hz20hz_2Range(3,:,[1:21],:,:);

% for iTimeRange = 1:size(SNR_merge820_2Range,5)
 	
%  	for  iChan = 1:size(SNR_merge820_2Range,2)

%  		for iCmp = 1:3

%  			[H,pvalue_SNR_merge820_2Range(iCmp,iChan,iTimeRange),CI] = ttest(squeeze(SNR_merge820_2Range(1,iChan,:,1,iTimeRange)), squeeze(SNR_merge820_2Range(1,iChan,:,iCmp+1,iTimeRange)));

%  		end
%  	end
% end


% SNR_merge820_2Range800 = SNR_8hz20hz_2Range800(3,:,[1:21],:,:);

% for iTimeRange = 1:size(SNR_merge820_2Range800,5)
 	
%  	for  iChan = 1:size(SNR_merge820_2Range800,2)

%  		for iCmp = 1:3

%  			[H,pvalue_SNR_merge820_2Range800(iCmp,iChan,iTimeRange),CI] = ttest(squeeze(SNR_merge820_2Range800(1,iChan,:,1,iTimeRange)), squeeze(SNR_merge820_2Range800(1,iChan,:,iCmp+1,iTimeRange)));

%  		end
%  	end
% end





% clear Power_merge820_3Range pvalue_Power_merge820_3Range
% clear Power_merge820_2Range pvalue_Power_merge820_2Range
% clear Power_merge820_2Range800 pvalue_Power_merge820_2Range800


% Power_merge820_3Range = Power_8hz20hz_3Range(3,:,[1:21],:,:);

% for iTimeRange = 1:size(Power_merge820_3Range,5)
 	
%  	for  iChan = 1:size(Power_merge820_3Range,2)

%  		for iCmp = 1:3

%  			[H,pvalue_Power_merge820_3Range(iCmp,iChan,iTimeRange),CI] = ttest(squeeze(Power_merge820_3Range(1,iChan,:,1,iTimeRange)), squeeze(Power_merge820_3Range(1,iChan,:,iCmp+1,iTimeRange)));

%  		end
%  	end
% end




% Power_merge820_2Range = Power_8hz20hz_2Range(3,:,[1:21],:,:);

% for iTimeRange = 1:size(Power_merge820_2Range,5)
 	
%  	for  iChan = 1:size(Power_merge820_2Range,2)

%  		for iCmp = 1:3

%  			[H,pvalue_Power_merge820_2Range(iCmp,iChan,iTimeRange),CI] = ttest(squeeze(Power_merge820_2Range(1,iChan,:,1,iTimeRange)), squeeze(Power_merge820_2Range(1,iChan,:,iCmp+1,iTimeRange)));

%  		end
%  	end
% end


% Power_merge820_2Range800 = Power_8hz20hz_2Range800(3,:,[1:21],:,:);

% for iTimeRange = 1:size(Power_merge820_2Range800,5)
 	
%  	for  iChan = 1:size(Power_merge820_2Range800,2)

%  		for iCmp = 1:3

%  			[H,pvalue_Power_merge820_2Range800(iCmp,iChan,iTimeRange),CI] = ttest(squeeze(Power_merge820_2Range800(1,iChan,:,1,iTimeRange)), squeeze(Power_merge820_2Range800(1,iChan,:,iCmp+1,iTimeRange)));

%  		end
%  	end
% end




clear PowerOfSNR_merge820_3Range pvalue_PowerOfSNR_merge820_3Range
clear PowerOfSNR_merge820_2Range pvalue_PowerOfSNR_merge820_2Range
clear PowerOfSNR_merge820_2Range800 pvalue_PowerOfSNR_merge820_2Range800


PowerOfSNR_merge820_3Range = PowerOfSNR_8hz20hz_3Range(3,:,[1:21],:,:);

for iTimeRange = 1:size(PowerOfSNR_merge820_3Range,5)
 	
 	for  iChan = 1:size(PowerOfSNR_merge820_3Range,2)

 		for iCmp = 1:3

 			[H,pvalue_PowerOfSNR_merge820_3Range(iCmp,iChan,iTimeRange),CI] = ttest(squeeze(PowerOfSNR_merge820_3Range(1,iChan,:,1,iTimeRange)), squeeze(PowerOfSNR_merge820_3Range(1,iChan,:,iCmp+1,iTimeRange)));

 		end
 	end
end




PowerOfSNR_merge820_2Range = PowerOfSNR_8hz20hz_2Range(3,:,[1:21],:,:);

for iTimeRange = 1:size(PowerOfSNR_merge820_2Range,5)
 	
 	for  iChan = 1:size(PowerOfSNR_merge820_2Range,2)

 		for iCmp = 1:3

 			[H,pvalue_PowerOfSNR_merge820_2Range(iCmp,iChan,iTimeRange),CI] = ttest(squeeze(PowerOfSNR_merge820_2Range(1,iChan,:,1,iTimeRange)), squeeze(PowerOfSNR_merge820_2Range(1,iChan,:,iCmp+1,iTimeRange)));

 		end
 	end
end


PowerOfSNR_merge820_2Range800 = PowerOfSNR_8hz20hz_2Range800(3,:,[1:21],:,:);

for iTimeRange = 1:size(PowerOfSNR_merge820_2Range800,5)
 	
 	for  iChan = 1:size(PowerOfSNR_merge820_2Range800,2)

 		for iCmp = 1:3

 			[H,pvalue_PowerOfSNR_merge820_2Range800(iCmp,iChan,iTimeRange),CI] = ttest(squeeze(PowerOfSNR_merge820_2Range800(1,iChan,:,1,iTimeRange)), squeeze(PowerOfSNR_merge820_2Range800(1,iChan,:,iCmp+1,iTimeRange)));

 		end
 	end
end


