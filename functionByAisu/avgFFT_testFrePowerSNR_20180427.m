
%% avgFFT 

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
%-------------- do fft of [-1000 1200] for all conditions-------------/
timeBin    = [-1000,2000]; % 
maxChanNum = 3;
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

save(['avgFFT_mergeAllconFindChan_-1000To1200','noldr_skip' num2str(snrskipNum)]); 




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










%%%% step3
%---- do fft of baseline and post-baseline([201 700;701 1200; 201:1200]) for 4 conditions,separately ------------------/
clear all;
directory     = 'D:\expData\SSVEP_8loc\Exp1_20170328\Data\scanData\20180407\average';
chanList      = 'all';
subNum        = [1:21];
FFTpoints     = 1024;
windowType    = 'cosin';
snrskipNum    = 1;
freRange      = [5 25];
targetFre     = [8;20]; 
condition     = {'allcue13_2~45Hz','allcue21_2~45Hz','allcue34_2~45Hz','allcue42_2~45Hz'};
timeBin       = [-500,-1;201,700;701 1200;201 1200;401 1200];
allPower_4con = [];
allSNR_4con   = [];
ldrTransType  = 'noldr';


for iCondition = 1:size(condition,2)

		cd(fullfile(directory,condition{iCondition}));

		for iSub = subNum

				iSubFileName = [num2str(iSub),'-',condition{iCondition},'.avg'];

				for iTimeRange = 1:size(timeBin,1)  

					[powerData,f,chan_names,chanList,rawSignal] = FFTanalysis_bcl(iSubFileName,timeBin(iTimeRange,:),chanList,FFTpoints,windowType,ldrTransType);

					allPower_4con_allTimeRange(:,:,iSub,iCondition,iTimeRange)  = powerData; % size: 1-freBins; 2-chan; 3-sub; 4-condition; 5-timeRange(base)

					[SNR] = computeSNR(powerData,snrskipNum);

					allSNR_4con_allTimeRange(:,:,iSub,iCondition,iTimeRange)    = SNR;

					clear SNR powerData;

				end
		end

		cd('..');
end 


[freIdx]                   = cutFreBins(f,freRange); 
allSNR_4con_allTimeRange   = allSNR_4con_allTimeRange(freIdx,:,:,:,:);
allPower_4con_allTimeRange = allPower_4con_allTimeRange(freIdx,:,:,:,:);
freBins                    = f(freIdx);
%-------find closest frePoint-------/
for iFre = 1:size(targetFre,1)
	[X, freIndex]         = sort(abs(freBins-targetFre(iFre,1)),'ascend');
	freTimePoint(iFre,1)  = freIndex(1);
	realFre(iFre,1) 	  = freBins(freIndex(1));
end
%-------have found, the end---------\



allSNRMinusBase_4con_allTimeRange(:,:,:,:,1) = allSNR_4con_allTimeRange(:,:,:,:,2)-allSNR_4con_allTimeRange(:,:,:,:,1); % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-timeRange(post1,2)
allSNRMinusBase_4con_allTimeRange(:,:,:,:,2) = allSNR_4con_allTimeRange(:,:,:,:,3)-allSNR_4con_allTimeRange(:,:,:,:,1); % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-timeRange(post1,2)
allSNRMinusBase_4con_allTimeRange(:,:,:,:,3) = allSNR_4con_allTimeRange(:,:,:,:,4)-allSNR_4con_allTimeRange(:,:,:,:,1); % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-timeRange(post1,2)
allSNRMinusBase_4con_allTimeRange(:,:,:,:,4) = allSNR_4con_allTimeRange(:,:,:,:,5)-allSNR_4con_allTimeRange(:,:,:,:,1); % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-timeRange(post1,2)

allFrePointSNRMinusBase_4con_allTimeRange = cat(1,allSNRMinusBase_4con_allTimeRange(freTimePoint(1),:,:,[2 3 4 1],:), allSNRMinusBase_4con_allTimeRange(freTimePoint(2),:,:,[3 2 1 4],:)); % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-timeRange(post1,2)
allFrePointSNRMinusBase_4con_3Range       = allFrePointSNRMinusBase_4con_allTimeRange(:,:,:,:,[1 2]);
allFrePointSNRMinusBase_4con_2Range       = allFrePointSNRMinusBase_4con_allTimeRange(:,:,:,:,[3]);
allFrePointSNRMinusBase_4con_2Range800    = allFrePointSNRMinusBase_4con_allTimeRange(:,:,:,:,[4]);



allPowerMinusBase_4con_allTimeRange(:,:,:,:,1) = allPower_4con_allTimeRange(:,:,:,:,2)-allPower_4con_allTimeRange(:,:,:,:,1); % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-timeRange(post1,2)
allPowerMinusBase_4con_allTimeRange(:,:,:,:,2) = allPower_4con_allTimeRange(:,:,:,:,3)-allPower_4con_allTimeRange(:,:,:,:,1); % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-timeRange(post1,2)
allPowerMinusBase_4con_allTimeRange(:,:,:,:,3) = allPower_4con_allTimeRange(:,:,:,:,4)-allPower_4con_allTimeRange(:,:,:,:,1); % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-timeRange(post1,2)
allPowerMinusBase_4con_allTimeRange(:,:,:,:,4) = allPower_4con_allTimeRange(:,:,:,:,5)-allPower_4con_allTimeRange(:,:,:,:,1); % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-timeRange(post1,2)

allFrePointPowerMinusBase_4con_allTimeRange = cat(1,allPowerMinusBase_4con_allTimeRange(freTimePoint(1),:,:,[2 3 4 1],:), allPowerMinusBase_4con_allTimeRange(freTimePoint(2),:,:,[3 2 1 4],:)); % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-timeRange(post1,2)
allFrePointPowerMinusBase_4con_3Range       = allFrePointPowerMinusBase_4con_allTimeRange(:,:,:,:,[1 2]);
allFrePointPowerMinusBase_4con_2Range       = allFrePointPowerMinusBase_4con_allTimeRange(:,:,:,:,[3]);
allFrePointPowerMinusBase_4con_2Range800    = allFrePointPowerMinusBase_4con_allTimeRange(:,:,:,:,[4]);



save(['fft_4con_allTimeRange_skip' num2str(snrskipNum) ]);
%---------- It's the end..  Do fft  for 4 conditions,separately. It's the end---------------\


%%%-------------- read eachSub's Max 3 chans' Data ----------------


clear all;
load('D:\expData\SSVEP_8loc\Exp1_20170328\Data\scanData\20180407\average\spssData_final_20180422\avgFFT_mergeAllconFindChan_-500To1200noldr_skip1.mat')
clearvars -except chooseChanSNR  chooseChanPower  ;
load('D:\expData\SSVEP_8loc\Exp1_20170328\Data\scanData\20180407\average\spssData_final_20180422\fft_4con_allTimeRange_skip1.mat')



Power_8hz20hz_3Range = [];
SNR_8hz20hz_3Range = [];



for iFre = 1:size(allFrePointPowerMinusBase_4con_3Range,1)  % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-TimeRange

	for iSub = 1:size(allFrePointPowerMinusBase_4con_3Range,3) 
			cPower = allFrePointPowerMinusBase_4con_3Range(iFre,chooseChanPower(iSub,[1:5],iFre),iSub,:,:); 
			Power_8hz20hz_3Range(iFre,:,iSub,:,:) = cPower;

			cSNR = allFrePointSNRMinusBase_4con_3Range(iFre,chooseChanSNR(iSub,[1:5],iFre),iSub,:,:); 
			SNR_8hz20hz_3Range(iFre,:,iSub,:,:) = cPower;

	end
end



Power_8hz20hz_3Range(3,:,:,:,:)                              = mean(Power_8hz20hz_3Range,1);
Power_8hz20hz_3Range(:,size(Power_8hz20hz_3Range,2)+1,:,:,:) = mean(Power_8hz20hz_3Range,2);


SNR_8hz20hz_3Range(3,:,:,:,:)                                = mean(SNR_8hz20hz_3Range,1);
SNR_8hz20hz_3Range(:,size(SNR_8hz20hz_3Range,2)+1,:,:,:)     = mean(SNR_8hz20hz_3Range,2);

%%%%%%%%%%%%   spss  %%%%%%%%%%
clear Power_merge820_3Range pvalue_Power_merge820_3Range
clear Power_merge820_2Range pvalue_Power_merge820_2Range
clear Power_merge820_2Range800 pvalue_Power_merge820_2Range800
clear SNR_merge820_3Range pvalue_SNR_merge820_3Range
clear SNR_merge820_2Range pvalue_SNR_merge820_2Range
clear SNR_merge820_2Range800 pvalue_SNR_merge820_2Range800

validSubNum =   [1:3,5:7,9:11,13:18,20,21];   %  [1:21]  %

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





%%%  for figure of dissertation

%% topo of 8hz SNR [-500 1200] 21subjects
% [myColormap] = makeColormap(20,0,1,'r');
% figure;
% set(gcf,'Position',[1 1 1600 1080*2/3],'color','w');
% for iSub = 1:size(allSNR,3)
% 	subplot(4,6,iSub);
% 	topoplot_bcl(allSNR(freTimePoint(1),:,iSub),'EEGChans',standard_eName,'maplimitsDouble',[0 3],'colormap',myColormap);
% end
% set(gcf,'color','w');
% print('-depsc','-painters',['D:\expData\SSVEP_8loc\Exp1_20170328\Figure\rawFigure\SNR_8hz.eps']);


% %% topo of 20hz SNR [-500 1200] 21subjects
% set(gcf,'Position',[1 1 1600 1080*2/3],'color','w');
% for iSub = 1:size(allSNR,3)
% 	subplot(4,6,iSub);
% 	topoplot_bcl(allSNR(freTimePoint(2),:,iSub),'EEGChans',standard_eName,'maplimitsDouble',[0 3],'colormap',myColormap);
% end
% set(gcf,'color','w');
% print('-depsc','-painters',['D:\expData\SSVEP_8loc\Exp1_20170328\Figure\rawFigure\SNR_20hz.eps']);


% %% topo of 8hz Power [-500 1200] 21subjects
% [myColormap] = makeColormap(20,0,1,'b');
% figure;
% set(gcf,'Position',[1 1 1600 1080*2/3],'color','w');
% for iSub = 1:size(allPower,3)
% 	subplot(4,6,iSub);
% 	topoplot_bcl(allPower(freTimePoint(1),:,iSub),'EEGChans',standard_eName,'maplimitsDouble',[0 0.1],'colormap',myColormap);
% end
% set(gcf,'color','w');
% print('-depsc','-painters',['D:\expData\SSVEP_8loc\Exp1_20170328\Figure\rawFigure\Power_8hz.eps']);

% %% topo of 20hz Power [-500 1200] 21subjects
% figure;
% set(gcf,'Position',[1 1 1600 1080*2/3],'color','w');
% for iSub = 1:size(allPower,3)
% 	subplot(4,6,iSub);
% 	topoplot_bcl(allPower(freTimePoint(2),:,iSub),'EEGChans',standard_eName,'maplimitsDouble',[0 0.05],'colormap',myColormap);
% end
% set(gcf,'color','w');
% print('-depsc','-painters',['D:\expData\SSVEP_8loc\Exp1_20170328\Figure\rawFigure\Power_20hz.eps']);




% topo of subAvgSNR 21sujects

validSubNum  = [1:21];
[myColormap] = makeColormap(20,0,1,'r');
subplot(1,2,1);
topoplot_bcl(squeeze(mean(allSNR(freTimePoint(1),:,validSubNum),3)),'EEGChans',standard_eName,'maplimitsDouble',[0 3],'colormap',myColormap);
subplot(1,2,2);
topoplot_bcl(squeeze(mean(allSNR(freTimePoint(2),:,validSubNum),3)),'EEGChans',standard_eName,'maplimitsDouble',[0 3],'colormap',myColormap);
set(gcf,'color','w');
print('-depsc','-painters',['D:\expData\SSVEP_8loc\Exp1_20170328\Figure\rawFigure\SNRavg_820.eps']);


% topo of subAvgPower
[myColormap] = makeColormap(20,0,1,'b');
figure;
topoplot_bcl(squeeze(mean(allPower(freTimePoint(1),:,validSubNum),3)),'EEGChans',standard_eName,'maplimitsDouble',[0 0.15],'colormap',myColormap);
set(gcf,'color','w');
print('-depsc','-painters',['D:\expData\SSVEP_8loc\Exp1_20170328\Figure\rawFigure\Poweravg_8hz.eps']);

figure;
topoplot_bcl(squeeze(mean(allPower(freTimePoint(2),:,validSubNum),3)),'EEGChans',standard_eName,'maplimitsDouble',[0 0.05],'colormap',myColormap);
set(gcf,'color','w');
print('-depsc','-painters',['D:\expData\SSVEP_8loc\Exp1_20170328\Figure\rawFigure\Poweravg_20hz.eps']);









% average Power of max 5chans for 8Hz & 20 Hz 
SNR_820=[squeeze(mean(MaxChanSNR_SNR(freTimePoint(1),:,:,1),2)) squeeze(mean(MaxChanPower_SNR(freTimePoint(2),:,:,2),2))];
[X, freIndex] = sort(SNR_820(:,1),'descend');
bar(SNR_820(freIndex,:),'grouped','EdgeColor','none');
xlim([0.3 21.5]);
set(gca,'xtick',[1:1:21]);  
set(gca,'xTickLabel',freIndex);
box off;
print('-depsc','-painters',['D:\expData\SSVEP_8loc\Exp1_20170328\Figure\rawFigure\SNR_820.eps']);



% % illustration of computeSNR
% [AX,H1,H2] = plotyy([1:11], MaxChanSNR_Power([10:20],3,10,2),[3:9],MaxChanSNR_SNR([12:18],3,10,2),'plot');
% set(AX(1),'XColor','k','YColor','b','Ylim',[0,0.2],'ytick',[0:0.05:0.2],'Xlim',[0.5 11.5]);		 				
% set(AX(2),'XColor','k','YColor','r','Ylim',[0,5],'ytick',[0:1:5],'Xlim',[0.5 11.5]);
% set(H1,'LineStyle','-','color','b');
% set(H2,'LineStyle','-','color','r');
% print('-depsc','-painters',['D:\expData\SSVEP_8loc\Exp1_20170328\Figure\rawFigure\computeSNR.eps']);






