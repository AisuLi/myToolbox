
%% avgFFT 

directory  = 'D:\expData\SSVEP_8loc\Exp1_20170328\Data\scanData\20180407\average';
chanList   = 'all';
subNum     = [1:9];
FFTpoints  = 1024;
windowType = 'cosin';
snrskipNum = 1;
freRange   = [5 25];
targetFre  = [8;20]; 
ldrTransType = 'noldr';


%%%% step1
%-------------- do fft of [-500 1200] for all conditions-------------/
timeBin    = [-500,1200];
maxChanNum = 5;
condition  = {'all_2~45Hz'};


for iSub = subNum

	iSubFileName = [num2str(iSub),'-',condition{:},'.avg'];
	
	[powerData,f,chan_names,chanList,rawSignal] = FFTanalysis_bcl_su(iSubFileName,timeBin,chanList,FFTpoints,windowType,ldrTransType);

	allPower(:,:,iSub) = powerData;

	[SNR] = computeSNR(powerData,snrskipNum);

	allSNR(:,:,iSub)  = SNR;  

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

save(['avgFFT_mergeAllconFindChan_-500To1200','noldr_allChan']); 




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
[P1eNo,standard_eName]=getElectrodeNo_bcl('haveM1',[]);

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




%%%%% step2:sort the Power
%%---------sort the SNR of allSub according to eachSub's eAvg-----------/
%%%%-----compute eAvg--------
frePointMaxChanPower(:,maxChanNum+1,:) = mean(frePointMaxChanPower,2);
for iFre = 1:size(frePointMaxChanPower,1)
		[X, subIndex]                           = sort(frePointMaxChanPower(iFre,maxChanNum+1,:),'descend');
		sortSubPower(iFre,:)                    = subIndex;
		frePointPower_sorted(iFre,:,:)          = frePointPower(iFre,:,subIndex);
		MaxChanPower_SNR_sorted(:,:,:,iFre)     = MaxChanPower_SNR(:,:,subIndex,iFre);
		MaxChanPower_Power_sorted(:,:,:,iFre)   = MaxChanPower_Power(:,:,subIndex,iFre);
		chooseChanPower_sorted(:,:,iFre)        = chooseChanPower(subIndex,:,iFre);
		clear X subIndex;
end

frePointPower_sorted(:,:,size(sortSubPower,2)+1)        = mean(frePointPower_sorted,3);

MaxChanPower_SNR_sorted(:,:,size(sortSubPower,2)+1,:)   = mean(MaxChanPower_SNR_sorted,3);
MaxChanPower_SNR_sorted(:,maxChanNum+1,:,:)             = mean(MaxChanPower_SNR_sorted,2);
MaxChanPower_Power_sorted(:,:,size(sortSubPower,2)+1,:) = mean(MaxChanPower_Power_sorted,3);
MaxChanPower_Power_sorted(:,maxChanNum+1,:,:)           = mean(MaxChanPower_Power_sorted,2);

waveNum                = maxChanNum+1;
[myColormap]           = makeColormap(20,0,1,'b');
subNumPerFig           = 5;
figNum                 = ceil(22./subNumPerFig);
[P1eNo,standard_eName] =getElectrodeNo_bcl('haveM1',[]);



for iPage = 1:figNum
	figure;
	set(gcf,'Position',get(0,'ScreenSize'),'color','w');

	% suptitle([fftType,'Power     -500To1200    chooseChan   left:    ', num2str(realFre(1)) 'Hz          right:    ',num2str(realFre(2))  'Hz']);
	
	cSub = [1:subNumPerFig]+(iPage-1)*subNumPerFig;

	for iSub = cSub
				

		for iFre = 1:size(targetFre,1)
				
			if iSub<=size(chooseChanPower,1)+1

				
				%-------plot topo-------/
				if iFre==1
					figureLoc = 2*(waveNum+1)*(iSub-1-(iPage-1)*subNumPerFig)+1;
				elseif iFre==2
					figureLoc = 2*(waveNum+1)*(iSub-1-(iPage-1)*subNumPerFig)+waveNum+2;
				end	

				subplot(subNumPerFig,2*(waveNum+1),figureLoc);
				if iSub<= size(chooseChanPower,1)
					topoplot_bcl(frePointPower_sorted(iFre,:,iSub),'EEGChans',standard_eName,'maplimitsDouble',[0 0.15],'colormap',myColormap);
				end
				%-----------------------\	

				for iWave = 1:waveNum
					figureLoc  = figureLoc+1;
					subplot(subNumPerFig,2*(waveNum+1),figureLoc);

					[AX,H1,H2] = plotyy(freBins,MaxChanPower_Power_sorted(:,iWave,iSub,iFre),freBins,MaxChanPower_SNR_sorted(:,iWave,iSub,iFre),'plot');
					
					if iSub<=size(chooseChanPower,1)
						if iWave<=maxChanNum
							title([num2str(sortSubPower(iFre,iSub)),standard_eName(squeeze(chooseChanPower_sorted(iSub,iWave,iFre)))]);
						else
							title(['eAvg']);					
						end
					else
						title(['allSubAvg']);
					end
					
					box off;

					set(AX(1),'XColor','k','YColor','b','Ylim',[0,0.25],'ytick',[0:0.05:0.25],'Xlim',freRange);		 				
					set(AX(2),'XColor','k','YColor','r','Ylim',[0,5],'ytick',[0:1:5],'Xlim',freRange);
					set(H1,'LineStyle','-','color','b');
					set(H2,'LineStyle','-','color','r');
					line(repmat(realFre(iFre,:),1,2),[0,.2],'linestyle',':','linewidth',2, 'color',[0 0 0 ]); 
					line([5 25],[0.05,0.05],'linestyle',':','linewidth',1, 'color',[0 0 0 ]); %
				end
			end
		end
	end
	legend([H1,H2],{'power';'Power'},'Location','North');
	legend('boxoff');
	set(gcf,'color','w');
	% saveas(gcf,['S',num2str(cSub(1)),'~',num2str(cSub(end)),'_',fftType,'_Powertopo.pdf']);
end
%%%%%%%%%%%%%%  have ploted sorted Power, the end %%%%%%%%%%%%%%%%%%%%%%%%%%%






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
ldrTransType = 'noldr';


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

allFrePointSNRMinusBase_4con_allTimeRange = cat(1,allSNRMinusBase_4con_allTimeRange(freTimePoint(1),:,:,[1 3 2 4],:), allSNRMinusBase_4con_allTimeRange(freTimePoint(2),:,:,[4 2 3 1],:)); % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-timeRange(post1,2)
allFrePointSNRMinusBase_4con_3Range       = allFrePointSNRMinusBase_4con_allTimeRange(:,:,:,:,[1 2]);
allFrePointSNRMinusBase_4con_2Range       = allFrePointSNRMinusBase_4con_allTimeRange(:,:,:,:,[3]);
allFrePointSNRMinusBase_4con_2Range800    = allFrePointSNRMinusBase_4con_allTimeRange(:,:,:,:,[4]);



allPowerMinusBase_4con_allTimeRange(:,:,:,:,1) = allPower_4con_allTimeRange(:,:,:,:,2)-allPower_4con_allTimeRange(:,:,:,:,1); % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-timeRange(post1,2)
allPowerMinusBase_4con_allTimeRange(:,:,:,:,2) = allPower_4con_allTimeRange(:,:,:,:,3)-allPower_4con_allTimeRange(:,:,:,:,1); % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-timeRange(post1,2)
allPowerMinusBase_4con_allTimeRange(:,:,:,:,3) = allPower_4con_allTimeRange(:,:,:,:,4)-allPower_4con_allTimeRange(:,:,:,:,1); % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-timeRange(post1,2)
allPowerMinusBase_4con_allTimeRange(:,:,:,:,4) = allPower_4con_allTimeRange(:,:,:,:,5)-allPower_4con_allTimeRange(:,:,:,:,1); % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-timeRange(post1,2)

allFrePointPowerMinusBase_4con_allTimeRange = cat(1,allPowerMinusBase_4con_allTimeRange(freTimePoint(1),:,:,[1 3 2 4],:), allPowerMinusBase_4con_allTimeRange(freTimePoint(2),:,:,[4 2 3 1],:)); % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-timeRange(post1,2)
allFrePointPowerMinusBase_4con_3Range       = allFrePointPowerMinusBase_4con_allTimeRange(:,:,:,:,[1 2]);
allFrePointPowerMinusBase_4con_2Range       = allFrePointPowerMinusBase_4con_allTimeRange(:,:,:,:,[3]);
allFrePointPowerMinusBase_4con_2Range800    = allFrePointPowerMinusBase_4con_allTimeRange(:,:,:,:,[4]);



save(['fft_4con_allTimeRange' ldrTransType]);
%---------- It's the end..  Do fft  for 4 conditions,separately. It's the end---------------\





%%%-------------- read eachSub's Max 3 chans' Data ----------------



clearvars -except chooseChanSNR chooseChanPower



%%%%%%%%% 3Range  %%%%%%%%%%
SNR_8hz20hz_3Range   = [];
Power_8hz20hz_3Range = [];

for iFre = 1:size(allFrePointSNRMinusBase_4con_3Range,1)  % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-TimeRange

	for iSub = 1:size(allFrePointSNRMinusBase_4con_3Range,3) 
			cSNR = allFrePointSNRMinusBase_4con_3Range(iFre,chooseChanSNR(iSub,:,iFre),iSub,:,:); 
			SNR_8hz20hz_3Range(iFre,:,iSub,:,:) = cSNR;

			cPower = allFrePointPowerMinusBase_4con_3Range(iFre,chooseChanPower(iSub,:,iFre),iSub,:,:); 
			Power_8hz20hz_3Range(iFre,:,iSub,:,:) = cPower;


	end
end

SNR_8hz20hz_3Range(3,:,:,:,:)                                                 = mean(SNR_8hz20hz_3Range,1);
SNR_8hz20hz_3Range(:,size(chooseChanSNR,2)+1,:,:,:)                           = mean(SNR_8hz20hz_3Range,2);
SNR_8hz20hz_3Range(:,:,size(chooseChanSNR,1)+1,:,:)                           = mean(SNR_8hz20hz_3Range,3);
SNR_8hz20hz_3Range(:,:,:,:,size(allFrePointSNRMinusBase_4con_3Range,5)+1)     = mean(SNR_8hz20hz_3Range,5);
Power_8hz20hz_3Range(3,:,:,:,:)                                               = mean(Power_8hz20hz_3Range,1);
Power_8hz20hz_3Range(:,size(chooseChanPower,2)+1,:,:,:)                       = mean(Power_8hz20hz_3Range,2);
Power_8hz20hz_3Range(:,:,size(chooseChanPower,1)+1,:,:)                       = mean(Power_8hz20hz_3Range,3);
Power_8hz20hz_3Range(:,:,:,:,size(allFrePointPowerMinusBase_4con_3Range,5)+1) = mean(Power_8hz20hz_3Range,5);


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
		end
	end
end

maxChanNum = 5;
%%-------------- start to plot -------------------/
ngroups    = size(SNR_8hz20hz_3Range,2);
nbars      = size(SNR_8hz20hz_3Range,4);
groupwidth = min(0.8, nbars/(nbars+1.5));	
barColor   = [235 95 52; 254 158 15;63 96 169; 0 0 200]/255;

% %%%%%%   bar of SNR %%%%%%%%
pptTitleName = ['SNR_max', num2str(maxChanNum), 'chan_3Range'];
mkdir(pptTitleName);
[opt,ppt] = figPPT(pptTitleName);

plotSNRRange = [-3 3; -1  1];
freStr = {'8hz','20hz','mean of 8+20'};
timeBinStr = {'201To700','701To1200','mean of 2Range'};

for iSub = 1:size(allFrePointSNRMinusBase_4con_3Range,3)+1
				
	figure;
    set (gcf,'Position',[50,50,1800,900], 'color','w');

	for iRange = 1:size(SNR_8hz20hz_3Range,5)		

		for iFre = 1:size(SNR_8hz20hz_3Range,1)  % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-TimeRange
		
			subplot(3,3,3*(iRange-1)+iFre);

			cSNR = squeeze(SNR_8hz20hz_3Range(iFre,:,iSub,:,iRange)); 
			% cError = SNR_8hz20hz_3Range(iFre,:,iSub+1,:,iRange);			
			H=bar(cSNR,'grouped','EdgeColor','none');
		    ch = get(H,'children');
	        set(ch{1},'Facecolor',barColor(1,:));
	        set(ch{2},'Facecolor',barColor(2,:));
	        set(ch{3},'Facecolor',barColor(3,:));
	        set(ch{4},'Facecolor',barColor(4,:));
			hold on;box off;
			ylim(plotSNRRange(1,:));

			if iSub==size(allFrePointSNRMinusBase_4con_3Range,3)+1
				subStr = 'allSubMean';
				cError = squeeze(errorBar_SNR_8hz20hz_3Range(iFre,:,:,iRange));
				for i = 1:nbars
				    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
				    errorbar(x,cSNR(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
				    hold on;
					ylim(plotSNRRange(2,:));
				end

			elseif iSub<=size(allFrePointSNRMinusBase_4con_3Range,3)
					subStr = ['Sub:',num2str(iSub)];
			end

			title([subStr,freStr(iFre),timeBinStr(iRange)]);
			clear cSNR cError;
		end
	end
	l = legend([ch{1} ch{2} ch{3} ch{4} ],'sepCue--Uncue','ajCue--Uncue','ajCue--Cue','SepCue--Cue','Location',[0.908333333333333 0.248888888888889 0.0816666666666667 0.0881481481481481]);
	legend('boxoff');
	cfigFileName = fullfile(fullfile(pwd,pptTitleName),[pptTitleName,'_',subStr([1:3,5:end]),'.bmp']);
	saveas(gcf,cfigFileName);
	ppt = ppt.addImageSlide(['SNR_',subStr(1:3),'_',subStr(5:end)],cfigFileName,opt);
	close all;
end

ppt.saveAs(fullfile(pwd,pptTitleName));

ppt.close;




% %%%%%%   bar of Power %%%%%%%%
pptTitleName = ['Power_max', num2str(maxChanNum), 'chan_3Range'];
mkdir(pptTitleName);
[opt,ppt] = figPPT(pptTitleName);

plotPowerRange = [-0.3 0.3; -0.1  0.1];
freStr = {'8hz','20hz','mean of 8+20'};
timeBinStr = {'201To700','701To1200','mean of 2Range'};

for iSub = 1:size(allFrePointPowerMinusBase_4con_3Range,3)+1
				
	figure;
    set (gcf,'Position',[50,50,1800,900], 'color','w');

	for iRange = 1:size(Power_8hz20hz_3Range,5)		

		for iFre = 1:size(Power_8hz20hz_3Range,1)  % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-TimeRange
		
			subplot(3,3,3*(iRange-1)+iFre);

			cPower = squeeze(Power_8hz20hz_3Range(iFre,:,iSub,:,iRange)); 
			% cError = Power_8hz20hz_3Range(iFre,:,iSub+1,:,iRange);			
			H=bar(cPower,'grouped','EdgeColor','none');
		    ch = get(H,'children');
	        set(ch{1},'Facecolor',barColor(1,:));
	        set(ch{2},'Facecolor',barColor(2,:));
	        set(ch{3},'Facecolor',barColor(3,:));
	        set(ch{4},'Facecolor',barColor(4,:));
			hold on;box off;
			ylim(plotPowerRange(1,:));

			if iSub==size(allFrePointPowerMinusBase_4con_3Range,3)+1
				subStr = 'allSubMean';
				cError = squeeze(errorBar_Power_8hz20hz_3Range(iFre,:,:,iRange));
				for i = 1:nbars
				    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
				    errorbar(x,cPower(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
				    hold on;
					ylim(plotPowerRange(2,:));
				end

			elseif iSub<=size(allFrePointPowerMinusBase_4con_3Range,3)
					subStr = ['Sub:',num2str(iSub)];
			end

			title([subStr,freStr(iFre),timeBinStr(iRange)]);
			clear cPower cError;
		end
	end
	l = legend([ch{1} ch{2} ch{3} ch{4} ],'sepCue--Uncue','ajCue--Uncue','ajCue--Cue','SepCue--Cue','Location',[0.908333333333333 0.248888888888889 0.0816666666666667 0.0881481481481481]);
	legend('boxoff');
	cfigFileName = fullfile(fullfile(pwd,pptTitleName),[pptTitleName,'_',subStr([1:3,5:end]),'.bmp']);
	saveas(gcf,cfigFileName);
	ppt = ppt.addImageSlide(['Power_',subStr(1:3),'_',subStr(5:end)],cfigFileName,opt);
	close all;
end

ppt.saveAs(fullfile(pwd,pptTitleName));

ppt.close;







%%%%% 2Range

%%%%%%%%% 2Range  %%%%%%%%%%
SNR_8hz20hz_2Range   = [];
Power_8hz20hz_2Range = [];

for iFre = 1:size(allFrePointSNRMinusBase_4con_2Range,1)  % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-TimeRange

	for iSub = 1:size(allFrePointSNRMinusBase_4con_2Range,3) 
			cSNR = allFrePointSNRMinusBase_4con_2Range(iFre,chooseChanSNR(iSub,:,iFre),iSub,:,:); 
			SNR_8hz20hz_2Range(iFre,:,iSub,:,:) = cSNR;

			cPower = allFrePointPowerMinusBase_4con_2Range(iFre,chooseChanPower(iSub,:,iFre),iSub,:,:); 
			Power_8hz20hz_2Range(iFre,:,iSub,:,:) = cPower;


	end
end

SNR_8hz20hz_2Range(3,:,:,:,:)                                                 = mean(SNR_8hz20hz_2Range,1);
SNR_8hz20hz_2Range(:,size(chooseChanSNR,2)+1,:,:,:)                           = mean(SNR_8hz20hz_2Range,2);
SNR_8hz20hz_2Range(:,:,size(chooseChanSNR,1)+1,:,:)                           = mean(SNR_8hz20hz_2Range,3);

Power_8hz20hz_2Range(3,:,:,:,:)                                               = mean(Power_8hz20hz_2Range,1);
Power_8hz20hz_2Range(:,size(chooseChanPower,2)+1,:,:,:)                       = mean(Power_8hz20hz_2Range,2);
Power_8hz20hz_2Range(:,:,size(chooseChanPower,1)+1,:,:)                       = mean(Power_8hz20hz_2Range,3);



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
		end
	end
end

maxChanNum = 5;
%%-------------- start to plot -------------------/
ngroups    = size(SNR_8hz20hz_2Range,2);
nbars      = size(SNR_8hz20hz_2Range,4);
groupwidth = min(0.8, nbars/(nbars+1.5));	
barColor   = [235 95 52; 254 158 15;63 96 169; 0 0 200]/255;

%%%%%%   bar of SNR %%%%%%%%
pptTitleName = ['SNR_max', num2str(maxChanNum), 'chan_2Range'];
mkdir(pptTitleName);
[opt,ppt] = figPPT(pptTitleName);

plotSNRRange = [-3 3; -1  1];
freStr = {'8hz','20hz','mean of 8+20'};
timeBinStr = {'201To1200'};

for iSub = 1:size(allFrePointSNRMinusBase_4con_2Range,3)+1
				
	figure;
    set (gcf,'Position',[50,50,1800,900], 'color','w');

	for iRange = 1:size(SNR_8hz20hz_2Range,5)		

		for iFre = 1:size(SNR_8hz20hz_2Range,1)  % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-TimeRange
		
			subplot(3,3,3*(iRange-1)+iFre);

			cSNR = squeeze(SNR_8hz20hz_2Range(iFre,:,iSub,:,iRange)); 
			% cError = SNR_8hz20hz_2Range(iFre,:,iSub+1,:,iRange);			
			H=bar(cSNR,'grouped','EdgeColor','none');
		    ch = get(H,'children');
	        set(ch{1},'Facecolor',barColor(1,:));
	        set(ch{2},'Facecolor',barColor(2,:));
	        set(ch{3},'Facecolor',barColor(3,:));
	        set(ch{4},'Facecolor',barColor(4,:));
			hold on;box off;
			ylim(plotSNRRange(1,:));

			if iSub==size(allFrePointSNRMinusBase_4con_2Range,3)+1
				subStr = 'allSubMean';
				cError = squeeze(errorBar_SNR_8hz20hz_2Range(iFre,:,:,iRange));
				for i = 1:nbars
				    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
				    errorbar(x,cSNR(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
				    hold on;
					ylim(plotSNRRange(2,:));
				end

			elseif iSub<=size(allFrePointSNRMinusBase_4con_2Range,3)
					subStr = ['Sub:',num2str(iSub)];
			end

			title([subStr,freStr(iFre),timeBinStr(iRange)]);
			clear cSNR cError;
		end
	end
	l = legend([ch{1} ch{2} ch{3} ch{4} ],'sepCue--Uncue','ajCue--Uncue','ajCue--Cue','SepCue--Cue','Location',[0.908333333333333 0.248888888888889 0.0816666666666667 0.0881481481481481]);
	legend('boxoff');
	cfigFileName = fullfile(fullfile(pwd,pptTitleName),[pptTitleName,'_',subStr([1:3,5:end]),'.bmp']);
	saveas(gcf,cfigFileName);
	ppt = ppt.addImageSlide(['SNR_',subStr(1:3),'_',subStr(5:end)],cfigFileName,opt);
	close all;
end

ppt.saveAs(fullfile(pwd,pptTitleName));

ppt.close;




%%%%%%   bar of Power %%%%%%%%
pptTitleName = ['Power_max', num2str(maxChanNum), 'chan_2Range'];
mkdir(pptTitleName);
[opt,ppt] = figPPT(pptTitleName);

plotPowerRange = [-0.3 0.3; -0.1  0.1];
freStr = {'8hz','20hz','mean of 8+20'};
timeBinStr = {'201To1200'};

for iSub = 1: size(allFrePointPowerMinusBase_4con_2Range,3)+1
				
	figure;
    set (gcf,'Position',[50,50,1800,900], 'color','w');

	for iRange = 1:size(Power_8hz20hz_2Range,5)		

		for iFre = 1:size(Power_8hz20hz_2Range,1)  % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-TimeRange
		
			subplot(3,3,3*(iRange-1)+iFre);

			cPower = squeeze(Power_8hz20hz_2Range(iFre,:,iSub,:,iRange)); 
			% cError = Power_8hz20hz_2Range(iFre,:,iSub+1,:,iRange);			
			H=bar(cPower,'grouped','EdgeColor','none');
		    ch = get(H,'children');
	        set(ch{1},'Facecolor',barColor(1,:));
	        set(ch{2},'Facecolor',barColor(2,:));
	        set(ch{3},'Facecolor',barColor(3,:));
	        set(ch{4},'Facecolor',barColor(4,:));
			hold on;box off;
			ylim(plotPowerRange(1,:));

			if iSub==size(allFrePointPowerMinusBase_4con_2Range,3)+1
				subStr = 'allSubMean';
				cError = squeeze(errorBar_Power_8hz20hz_2Range(iFre,:,:,iRange));
				for i = 1:nbars
				    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
				    errorbar(x,cPower(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
				    hold on;
					ylim(plotPowerRange(2,:));
				end

			elseif iSub<=size(allFrePointPowerMinusBase_4con_2Range,3)
					subStr = ['Sub:',num2str(iSub)];
			end

			title([subStr,freStr(iFre),timeBinStr(iRange)]);
			clear cPower cError;
		end
	end
	l = legend([ch{1} ch{2} ch{3} ch{4} ],'sepCue--Uncue','ajCue--Uncue','ajCue--Cue','SepCue--Cue','Location',[0.908333333333333 0.248888888888889 0.0816666666666667 0.0881481481481481]);
	legend('boxoff');
	cfigFileName = fullfile(fullfile(pwd,pptTitleName),[pptTitleName,'_',subStr([1:3,5:end]),'.bmp']);
	saveas(gcf,cfigFileName);
	ppt = ppt.addImageSlide(['Power_',subStr(1:3),'_',subStr(5:end)],cfigFileName,opt);
	close all;
end

ppt.saveAs(fullfile(pwd,pptTitleName));

ppt.close;












%%%%% 2Range800

%%%%%%%%% 2Range800  %%%%%%%%%%
SNR_8hz20hz_2Range800   = [];
Power_8hz20hz_2Range800 = [];

for iFre = 1:size(allFrePointSNRMinusBase_4con_2Range800,1)  % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-TimeRange

	for iSub = 1:size(allFrePointSNRMinusBase_4con_2Range800,3) 
			cSNR = allFrePointSNRMinusBase_4con_2Range800(iFre,chooseChanSNR(iSub,:,iFre),iSub,:,:); 
			SNR_8hz20hz_2Range800(iFre,:,iSub,:,:) = cSNR;

			cPower = allFrePointPowerMinusBase_4con_2Range800(iFre,chooseChanPower(iSub,:,iFre),iSub,:,:); 
			Power_8hz20hz_2Range800(iFre,:,iSub,:,:) = cPower;


	end
end

SNR_8hz20hz_2Range800(3,:,:,:,:)                                                 = mean(SNR_8hz20hz_2Range800,1);
SNR_8hz20hz_2Range800(:,size(chooseChanSNR,2)+1,:,:,:)                           = mean(SNR_8hz20hz_2Range800,2);
SNR_8hz20hz_2Range800(:,:,size(chooseChanSNR,1)+1,:,:)                           = mean(SNR_8hz20hz_2Range800,3);

Power_8hz20hz_2Range800(3,:,:,:,:)                                               = mean(Power_8hz20hz_2Range800,1);
Power_8hz20hz_2Range800(:,size(chooseChanPower,2)+1,:,:,:)                       = mean(Power_8hz20hz_2Range800,2);
Power_8hz20hz_2Range800(:,:,size(chooseChanPower,1)+1,:,:)                       = mean(Power_8hz20hz_2Range800,3);



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
		end
	end
end

maxChanNum = 5;
%%-------------- start to plot -------------------/
ngroups    = size(SNR_8hz20hz_2Range800,2);
nbars      = size(SNR_8hz20hz_2Range800,4);
groupwidth = min(0.8, nbars/(nbars+1.5));	
barColor   = [235 95 52; 254 158 15;63 96 169; 0 0 200]/255;

%%%%%%   bar of SNR %%%%%%%%
pptTitleName = ['SNR_max', num2str(maxChanNum), 'chan_2Range800'];
mkdir(pptTitleName);
[opt,ppt] = figPPT(pptTitleName);

plotSNRRange = [-3 3; -1  1];
freStr = {'8hz','20hz','mean of 8+20'};
timeBinStr = {'401To1200'};

for iSub = 1:size(allFrePointSNRMinusBase_4con_2Range800,3)+1
				
	figure;
    set (gcf,'Position',[50,50,1800,900], 'color','w');

	for iRange = 1:size(SNR_8hz20hz_2Range800,5)		

		for iFre = 1:size(SNR_8hz20hz_2Range800,1)  % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-TimeRange
		
			subplot(3,3,3*(iRange-1)+iFre);

			cSNR = squeeze(SNR_8hz20hz_2Range800(iFre,:,iSub,:,iRange)); 
			% cError = SNR_8hz20hz_2Range800(iFre,:,iSub+1,:,iRange);			
			H=bar(cSNR,'grouped','EdgeColor','none');
		    ch = get(H,'children');
	        set(ch{1},'Facecolor',barColor(1,:));
	        set(ch{2},'Facecolor',barColor(2,:));
	        set(ch{3},'Facecolor',barColor(3,:));
	        set(ch{4},'Facecolor',barColor(4,:));
			hold on;box off;
			ylim(plotSNRRange(1,:));

			if iSub==size(allFrePointSNRMinusBase_4con_2Range800,3)+1
				subStr = 'allSubMean';
				cError = squeeze(errorBar_SNR_8hz20hz_2Range800(iFre,:,:,iRange));
				for i = 1:nbars
				    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
				    errorbar(x,cSNR(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
				    hold on;
					ylim(plotSNRRange(2,:));
				end

			elseif iSub<=size(allFrePointSNRMinusBase_4con_2Range800,3)
					subStr = ['Sub:',num2str(iSub)];
			end

			title([subStr,freStr(iFre),timeBinStr(iRange)]);
			clear cSNR cError;
		end
	end
	l = legend([ch{1} ch{2} ch{3} ch{4} ],'sepCue--Uncue','ajCue--Uncue','ajCue--Cue','SepCue--Cue','Location',[0.908333333333333 0.248888888888889 0.0816666666666667 0.0881481481481481]);
	legend('boxoff');
	cfigFileName = fullfile(fullfile(pwd,pptTitleName),[pptTitleName,'_',subStr([1:3,5:end]),'.bmp']);
	saveas(gcf,cfigFileName);
	ppt = ppt.addImageSlide(['SNR_',subStr(1:3),'_',subStr(5:end)],cfigFileName,opt);
	close all;
end

ppt.saveAs(fullfile(pwd,pptTitleName));

ppt.close;




%%%%%%   bar of Power %%%%%%%%
pptTitleName = ['Power_max', num2str(maxChanNum), 'chan_2Range800'];
mkdir(pptTitleName);
[opt,ppt] = figPPT(pptTitleName);

plotPowerRange = [-0.3 0.3; -0.1  0.1];
freStr = {'8hz','20hz','mean of 8+20'};
timeBinStr = {'401To1200'};

for iSub = 1: size(allFrePointPowerMinusBase_4con_2Range800,3)+1
				
	figure;
    set (gcf,'Position',[50,50,1800,900], 'color','w');

	for iRange = 1:size(Power_8hz20hz_2Range800,5)		

		for iFre = 1:size(Power_8hz20hz_2Range800,1)  % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-TimeRange
		
			subplot(3,3,3*(iRange-1)+iFre);

			cPower = squeeze(Power_8hz20hz_2Range800(iFre,:,iSub,:,iRange)); 
			% cError = Power_8hz20hz_2Range800(iFre,:,iSub+1,:,iRange);			
			H=bar(cPower,'grouped','EdgeColor','none');
		    ch = get(H,'children');
	        set(ch{1},'Facecolor',barColor(1,:));
	        set(ch{2},'Facecolor',barColor(2,:));
	        set(ch{3},'Facecolor',barColor(3,:));
	        set(ch{4},'Facecolor',barColor(4,:));
			hold on;box off;
			ylim(plotPowerRange(1,:));

			if iSub==size(allFrePointPowerMinusBase_4con_2Range800,3)+1
				subStr = 'allSubMean';
				cError = squeeze(errorBar_Power_8hz20hz_2Range800(iFre,:,:,iRange));
				for i = 1:nbars
				    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
				    errorbar(x,cPower(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
				    hold on;
					ylim(plotPowerRange(2,:));
				end

			elseif iSub<=size(allFrePointPowerMinusBase_4con_2Range800,3)
					subStr = ['Sub:',num2str(iSub)];
			end

			title([subStr,freStr(iFre),timeBinStr(iRange)]);
			clear cPower cError;
		end
	end
	l = legend([ch{1} ch{2} ch{3} ch{4} ],'sepCue--Uncue','ajCue--Uncue','ajCue--Cue','SepCue--Cue','Location',[0.908333333333333 0.248888888888889 0.0816666666666667 0.0881481481481481]);
	legend('boxoff');
	cfigFileName = fullfile(fullfile(pwd,pptTitleName),[pptTitleName,'_',subStr([1:3,5:end]),'.bmp']);
	saveas(gcf,cfigFileName);
	ppt = ppt.addImageSlide(['Power_',subStr(1:3),'_',subStr(5:end)],cfigFileName,opt);
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



clear SNR_merge820_3Range pvalue_SNR_merge820_3Range
clear SNR_merge820_2Range pvalue_SNR_merge820_2Range
clear SNR_merge820_2Range800 pvalue_SNR_merge820_2Range800


SNR_merge820_3Range = SNR_8hz20hz_3Range(3,:,[1:21],:,:);

for iTimeRange = 1:size(SNR_merge820_3Range,5)
 	
 	for  iChan = 1:size(SNR_merge820_3Range,2)

 		for iCmp = 1:3

 			[H,pvalue_SNR_merge820_3Range(iCmp,iChan,iTimeRange),CI] = ttest(squeeze(SNR_merge820_3Range(1,iChan,:,1,iTimeRange)), squeeze(SNR_merge820_3Range(1,iChan,:,iCmp+1,iTimeRange)));

 		end
 	end
end




SNR_merge820_2Range = SNR_8hz20hz_2Range(3,:,[1:21],:,:);

for iTimeRange = 1:size(SNR_merge820_2Range,5)
 	
 	for  iChan = 1:size(SNR_merge820_2Range,2)

 		for iCmp = 1:3

 			[H,pvalue_SNR_merge820_2Range(iCmp,iChan,iTimeRange),CI] = ttest(squeeze(SNR_merge820_2Range(1,iChan,:,1,iTimeRange)), squeeze(SNR_merge820_2Range(1,iChan,:,iCmp+1,iTimeRange)));

 		end
 	end
end


SNR_merge820_2Range800 = SNR_8hz20hz_2Range800(3,:,[1:21],:,:);

for iTimeRange = 1:size(SNR_merge820_2Range800,5)
 	
 	for  iChan = 1:size(SNR_merge820_2Range800,2)

 		for iCmp = 1:3

 			[H,pvalue_SNR_merge820_2Range800(iCmp,iChan,iTimeRange),CI] = ttest(squeeze(SNR_merge820_2Range800(1,iChan,:,1,iTimeRange)), squeeze(SNR_merge820_2Range800(1,iChan,:,iCmp+1,iTimeRange)));

 		end
 	end
end





clear Power_merge820_3Range pvalue_Power_merge820_3Range
clear Power_merge820_2Range pvalue_Power_merge820_2Range
clear Power_merge820_2Range800 pvalue_Power_merge820_2Range800


Power_merge820_3Range = Power_8hz20hz_3Range(3,:,[1:21],:,:);

for iTimeRange = 1:size(Power_merge820_3Range,5)
 	
 	for  iChan = 1:size(Power_merge820_3Range,2)

 		for iCmp = 1:3

 			[H,pvalue_Power_merge820_3Range(iCmp,iChan,iTimeRange),CI] = ttest(squeeze(Power_merge820_3Range(1,iChan,:,1,iTimeRange)), squeeze(Power_merge820_3Range(1,iChan,:,iCmp+1,iTimeRange)));

 		end
 	end
end




Power_merge820_2Range = Power_8hz20hz_2Range(3,:,[1:21],:,:);

for iTimeRange = 1:size(Power_merge820_2Range,5)
 	
 	for  iChan = 1:size(Power_merge820_2Range,2)

 		for iCmp = 1:3

 			[H,pvalue_Power_merge820_2Range(iCmp,iChan,iTimeRange),CI] = ttest(squeeze(Power_merge820_2Range(1,iChan,:,1,iTimeRange)), squeeze(Power_merge820_2Range(1,iChan,:,iCmp+1,iTimeRange)));

 		end
 	end
end


Power_merge820_2Range800 = Power_8hz20hz_2Range800(3,:,[1:21],:,:);

for iTimeRange = 1:size(Power_merge820_2Range800,5)
 	
 	for  iChan = 1:size(Power_merge820_2Range800,2)

 		for iCmp = 1:3

 			[H,pvalue_Power_merge820_2Range800(iCmp,iChan,iTimeRange),CI] = ttest(squeeze(Power_merge820_2Range800(1,iChan,:,1,iTimeRange)), squeeze(Power_merge820_2Range800(1,iChan,:,iCmp+1,iTimeRange)));

 		end
 	end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  first avg, then choose chan, so all the subjects choose the same chans  %%%%%%%%%

% load('D:\expData\SSVEP_8loc\Exp1_20170328\Data\scanData\20180407\average\all_2~45Hz\avgFFT_mergeAllconFindChan_-500To1200.mat')

clearvars -except frePointSNR frePointPower allPower allSNR freBins freRange realFre
[P1eNo,standard_eName]=getElectrodeNo_bcl('haveM1',[]);
maxChanNum = 5;

frePointSNR(:,:,size(frePointSNR,3)+1)     = mean(frePointSNR,3);
frePointPower(:,:,size(frePointPower,3)+1) = mean(frePointPower,3);

allSNR(:,:,size(allSNR,3)+1)			   = mean(allSNR,3);
allPower(:,:,size(allPower,3)+1)		   = mean(allPower,3);
 

for iFre = 1:size(frePointSNR,1)

	[X, chanIndex]                = sort(frePointSNR(iFre,[1:size(standard_eName,2)-2,end]),'descend');	

	% chanIndex                     = chanIndex+23;
	chooseChanSNR_subSame(iFre,:) = chanIndex(1:maxChanNum);
	frePointMaxChanSNR(iFre,:,:)  = frePointSNR(iFre,chanIndex(1:maxChanNum),:);	
	MaxChanSNR_SNR(:,:,:,iFre)    = allSNR(:,chanIndex(1:maxChanNum),:);
	MaxChanSNR_Power(:,:,:,iFre)  = allPower(:,chanIndex(1:maxChanNum),:);
	clear X chanIndex

	[X, chanIndex]                  = sort(frePointPower(iFre,[1:size(standard_eName,2)-2,end]),'descend');	

	% chanIndex                       = chanIndex+23;
	chooseChanPower_subSame(iFre,:) = chanIndex(1:maxChanNum);
	frePointMaxChanPower(iFre,:,:)  = frePointPower(iFre,chanIndex(1:maxChanNum),:);	
	MaxChanPower_SNR(:,:,:,iFre)    = allSNR(:,chanIndex(1:maxChanNum),:);
	MaxChanPower_Power(:,:,:,iFre)  = allPower(:,chanIndex(1:maxChanNum),:);
end




% %%%%% step2:sort the SNR
% %%---------sort the SNR of allSub according to eachSub's eAvg-----------/
% %%%%-----compute eAvg--------
% frePointMaxChanSNR(:,maxChanNum+1,:) = mean(frePointMaxChanSNR,2);


% for iFre = 1:size(frePointMaxChanSNR,1)
% 		[X, subIndex]                       = sort(frePointMaxChanSNR(iFre,maxChanNum+1,[1:21]),'descend');
% 		sortSubSNR(iFre,:)                  = subIndex;
% 		frePointSNR_sorted(iFre,:,:)        = frePointSNR(iFre,:,subIndex);
% 		MaxChanSNR_SNR_sorted(:,:,:,iFre)   = MaxChanSNR_SNR(:,:,subIndex,iFre);
% 		MaxChanSNR_Power_sorted(:,:,:,iFre) = MaxChanSNR_Power(:,:,subIndex,iFre);

% end


% MaxChanSNR_Power_sorted(:,maxChanNum+1,:,:) = mean(MaxChanSNR_Power_sorted,2);
% MaxChanSNR_SNR_sorted(:,maxChanNum+1,:,:)   = mean(MaxChanSNR_SNR_sorted,2);

% waveNum                = maxChanNum+1;
% [myColormap]           = makeColormap(20,0,1,'r');
% subNumPerFig           = 5;
% figNum                 = ceil(22./subNumPerFig);
% [P1eNo,standard_eName] =getElectrodeNo_bcl('haveM1',[]);

% for iPage = 1:figNum
% 	figure;
% 	set(gcf,'Position',get(0,'ScreenSize'),'color','w');

% 	% suptitle([fftType,'SNR     -500To1200    chooseChan   left:    ', num2str(realFre(1)) 'Hz          right:    ',num2str(realFre(2))  'Hz']);
	
% 	cSub = [1:subNumPerFig]+(iPage-1)*subNumPerFig;

% 	for iSub = cSub
				
% 		for iFre = 1:2
				
% 			if iSub<=21

				
% 				%-------plot topo-------/
% 				if iFre==1
% 					figureLoc = 2*(waveNum+1)*(iSub-1-(iPage-1)*subNumPerFig)+1;
% 				elseif iFre==2
% 					figureLoc = 2*(waveNum+1)*(iSub-1-(iPage-1)*subNumPerFig)+waveNum+2;
% 				end	

% 				subplot(subNumPerFig,2*(waveNum+1),figureLoc);
% 				if iSub<= 21
% 					topoplot_bcl(frePointSNR_sorted(iFre,:,iSub),'EEGChans',standard_eName,'maplimitsDouble',[0 4],'colormap',myColormap);
% 				end
% 				%-----------------------\	

% 				for iWave = 1:waveNum
% 					figureLoc  = figureLoc+1;
% 					subplot(subNumPerFig,2*(waveNum+1),figureLoc);

% 					[AX,H1,H2] = plotyy(freBins,MaxChanSNR_Power_sorted(:,iWave,iSub,iFre),freBins,MaxChanSNR_SNR_sorted(:,iWave,iSub,iFre),'plot');
					
% 					if iSub<=21
% 						if iWave<=maxChanNum
% 							title([num2str(sortSubSNR(iFre,iSub)),standard_eName(squeeze(chooseChanSNR_subSame(iFre,iWave)))]);
% 						else
% 							title(['eAvg']);					
% 						end
% 					else
% 						title(['allSubAvg']);
% 					end
					
% 					box off;

% 					set(AX(1),'XColor','k','YColor','b','Ylim',[0,0.2],'ytick',[0:0.05:0.2],'Xlim',freRange);		 				
% 					set(AX(2),'XColor','k','YColor','r','Ylim',[0,4],'ytick',[0:1:4],'Xlim',freRange);
% 					set(H1,'LineStyle','-','color','b');
% 					set(H2,'LineStyle','-','color','r');
% 					line(repmat(realFre(iFre,:),1,2),[0,.2],'linestyle',':','linewidth',2, 'color',[0 0 0 ]); %
% 					line([5 25],[0.05,0.05],'linestyle',':','linewidth',1, 'color',[0 0 0 ]); %
% 				end
% 			end
% 		end
% 	end
% 	legend([H1,H2],{'power';'SNR'},'Location','North');
% 	legend('boxoff');
% 	set(gcf,'color','w');
% 	% saveas(gcf,['S',num2str(cSub(1)),'~',num2str(cSub(end)),'_',fftType,'_SNRtopo.pdf']);
% end

% %%%%%%%%%%%%%%  have ploted sorted SNR, the end %%%%%%%%%%%%%%%%%%%%%%%%%%%










% load('D:\expData\SSVEP_8loc\Exp1_20170328\Data\scanData\20180407\average\fft_4con_allTimeRange.mat')




%%%%%%%%% 3Range  %%%%%%%%%%
SNR_8hz20hz_3Range   = [];
Power_8hz20hz_3Range = [];

for iFre = 1:size(allFrePointSNRMinusBase_4con_3Range,1)  % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-TimeRange


			SNR_8hz20hz_3Range(iFre,:,:,:,:) = allFrePointSNRMinusBase_4con_3Range(iFre,chooseChanSNR_subSame(iFre,:),:,:,:); 

			Power_8hz20hz_3Range(iFre,:,:,:,:) = allFrePointPowerMinusBase_4con_3Range(iFre,chooseChanPower_subSame(iFre,:),:,:,:); 

end

SNR_8hz20hz_3Range(3,:,:,:,:)                                                 = mean(SNR_8hz20hz_3Range,1);
SNR_8hz20hz_3Range(:,size(chooseChanSNR_subSame,2)+1,:,:,:)                   = mean(SNR_8hz20hz_3Range,2);
SNR_8hz20hz_3Range(:,:,end+1,:,:)                                             = mean(SNR_8hz20hz_3Range,3);
SNR_8hz20hz_3Range(:,:,:,:,size(allFrePointSNRMinusBase_4con_3Range,5)+1)     = mean(SNR_8hz20hz_3Range,5);
Power_8hz20hz_3Range(3,:,:,:,:)                                               = mean(Power_8hz20hz_3Range,1);
Power_8hz20hz_3Range(:,size(chooseChanPower_subSame,2)+1,:,:,:)               = mean(Power_8hz20hz_3Range,2);
Power_8hz20hz_3Range(:,:,end+1,:,:)                                           = mean(Power_8hz20hz_3Range,3);
Power_8hz20hz_3Range(:,:,:,:,size(allFrePointPowerMinusBase_4con_3Range,5)+1) = mean(Power_8hz20hz_3Range,5);


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
		end
	end
end

maxChanNum = 5;
%%-------------- start to plot -------------------/
ngroups    = size(SNR_8hz20hz_3Range,2);
nbars      = size(SNR_8hz20hz_3Range,4);
groupwidth = min(0.8, nbars/(nbars+1.5));	
barColor   = [235 95 52; 254 158 15;63 96 169; 0 0 200]/255;

%%%%%%   bar of SNR %%%%%%%%
pptTitleName = ['SNR_max', num2str(maxChanNum), 'chanSubSame_3Range'];
mkdir(pptTitleName);
[opt,ppt] = figPPT(pptTitleName);

plotSNRRange = [-3 3; -1  1];
freStr = {'8hz','20hz','mean of 8+20'};
timeBinStr = {'201To700','701To1200','mean of 2Range'};

for iSub = 1:size(allFrePointSNRMinusBase_4con_3Range,3)+1
				
	figure;
    set (gcf,'Position',[50,50,1800,900], 'color','w');

	for iRange = 1:size(SNR_8hz20hz_3Range,5)		

		for iFre = 1:size(SNR_8hz20hz_3Range,1)  % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-TimeRange
		
			subplot(3,3,3*(iRange-1)+iFre);

			cSNR = squeeze(SNR_8hz20hz_3Range(iFre,:,iSub,:,iRange)); 
			% cError = SNR_8hz20hz_3Range(iFre,:,iSub+1,:,iRange);			
			H=bar(cSNR,'grouped','EdgeColor','none');
		    ch = get(H,'children');
	        set(ch{1},'Facecolor',barColor(1,:));
	        set(ch{2},'Facecolor',barColor(2,:));
	        set(ch{3},'Facecolor',barColor(3,:));
	        set(ch{4},'Facecolor',barColor(4,:));
			hold on;box off;
			ylim(plotSNRRange(1,:));

			if iSub==size(allFrePointSNRMinusBase_4con_3Range,3)+1
				subStr = 'allSubMean';
				cError = squeeze(errorBar_SNR_8hz20hz_3Range(iFre,:,:,iRange));
				for i = 1:nbars
				    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
				    errorbar(x,cSNR(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
				    hold on;
					ylim(plotSNRRange(2,:));
				end

			elseif iSub<=size(allFrePointSNRMinusBase_4con_3Range,3)
					subStr = ['Sub:',num2str(iSub)];
			end

			title([subStr,freStr(iFre),timeBinStr(iRange)]);
			clear cSNR cError;
		end
	end
	l = legend([ch{1} ch{2} ch{3} ch{4} ],'sepCue--Uncue','ajCue--Uncue','ajCue--Cue','SepCue--Cue','Location',[0.908333333333333 0.248888888888889 0.0816666666666667 0.0881481481481481]);
	legend('boxoff');
	cfigFileName = fullfile(fullfile(pwd,pptTitleName),[pptTitleName,'_',subStr([1:3,5:end]),'.bmp']);
	saveas(gcf,cfigFileName);
	ppt = ppt.addImageSlide(['SNR_',subStr(1:3),'_',subStr(5:end)],cfigFileName,opt);
	close all;
end

ppt.saveAs(fullfile(pwd,pptTitleName));

ppt.close;




%%%%%%   bar of Power %%%%%%%%
pptTitleName = ['Power_max', num2str(maxChanNum), 'chanSubSame_3Range'];
mkdir(pptTitleName);
[opt,ppt] = figPPT(pptTitleName);

plotPowerRange = [-0.3 0.3; -0.1  0.1];
freStr = {'8hz','20hz','mean of 8+20'};
timeBinStr = {'201To700','701To1200','mean of 2Range'};

for iSub = 1: size(allFrePointPowerMinusBase_4con_3Range,3)+1
				
	figure;
    set (gcf,'Position',[50,50,1800,900], 'color','w');

	for iRange = 1:size(Power_8hz20hz_3Range,5)		

		for iFre = 1:size(Power_8hz20hz_3Range,1)  % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-TimeRange
		
			subplot(3,3,3*(iRange-1)+iFre);

			cPower = squeeze(Power_8hz20hz_3Range(iFre,:,iSub,:,iRange)); 
			% cError = Power_8hz20hz_3Range(iFre,:,iSub+1,:,iRange);			
			H=bar(cPower,'grouped','EdgeColor','none');
		    ch = get(H,'children');
	        set(ch{1},'Facecolor',barColor(1,:));
	        set(ch{2},'Facecolor',barColor(2,:));
	        set(ch{3},'Facecolor',barColor(3,:));
	        set(ch{4},'Facecolor',barColor(4,:));
			hold on;box off;
			ylim(plotPowerRange(1,:));

			if iSub==size(allFrePointPowerMinusBase_4con_3Range,3)+1
				subStr = 'allSubMean';
				cError = squeeze(errorBar_Power_8hz20hz_3Range(iFre,:,:,iRange));
				for i = 1:nbars
				    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
				    errorbar(x,cPower(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
				    hold on;
					ylim(plotPowerRange(2,:));
				end

			elseif iSub<=size(allFrePointPowerMinusBase_4con_3Range,3)
					subStr = ['Sub:',num2str(iSub)];
			end

			title([subStr,freStr(iFre),timeBinStr(iRange)]);
			clear cPower cError;
		end
	end
	l = legend([ch{1} ch{2} ch{3} ch{4} ],'sepCue--Uncue','ajCue--Uncue','ajCue--Cue','SepCue--Cue','Location',[0.908333333333333 0.248888888888889 0.0816666666666667 0.0881481481481481]);
	legend('boxoff');
	cfigFileName = fullfile(fullfile(pwd,pptTitleName),[pptTitleName,'_',subStr([1:3,5:end]),'.bmp']);
	saveas(gcf,cfigFileName);
	ppt = ppt.addImageSlide(['Power_',subStr(1:3),'_',subStr(5:end)],cfigFileName,opt);
	close all;
end

ppt.saveAs(fullfile(pwd,pptTitleName));

ppt.close;




%%%%%%%%% 2Range  %%%%%%%%%%
SNR_8hz20hz_2Range   = [];
Power_8hz20hz_2Range = [];

for iFre = 1:size(allFrePointSNRMinusBase_4con_2Range,1)  % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-TimeRange

	for iSub = 1:size(allFrePointSNRMinusBase_4con_2Range,3) 
			cSNR = allFrePointSNRMinusBase_4con_2Range(iFre,chooseChanSNR_subSame(iFre,:),iSub,:,:); 
			SNR_8hz20hz_2Range(iFre,:,iSub,:,:) = cSNR;

			cPower = allFrePointPowerMinusBase_4con_2Range(iFre,chooseChanPower_subSame(iFre,:),iSub,:,:); 
			Power_8hz20hz_2Range(iFre,:,iSub,:,:) = cPower;

	end
end

SNR_8hz20hz_2Range(3,:,:,:,:)                                   = mean(SNR_8hz20hz_2Range,1);
SNR_8hz20hz_2Range(:,size(chooseChanSNR_subSame,2)+1,:,:,:)     = mean(SNR_8hz20hz_2Range,2);
SNR_8hz20hz_2Range(:,:,22,:,:)                                  = mean(SNR_8hz20hz_2Range,3);

Power_8hz20hz_2Range(3,:,:,:,:)                                 = mean(Power_8hz20hz_2Range,1);
Power_8hz20hz_2Range(:,size(chooseChanPower_subSame,2)+1,:,:,:) = mean(Power_8hz20hz_2Range,2);
Power_8hz20hz_2Range(:,:,22,:,:)                                = mean(Power_8hz20hz_2Range,3);



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
		end
	end
end

maxChanNum = 5;
%%-------------- start to plot -------------------/
ngroups    = size(SNR_8hz20hz_2Range,2);
nbars      = size(SNR_8hz20hz_2Range,4);
groupwidth = min(0.8, nbars/(nbars+1.5));	
barColor   = [235 95 52; 254 158 15;63 96 169; 0 0 200]/255;



%%%%%%   bar of SNR %%%%%%%%
pptTitleName = ['SNR_max', num2str(maxChanNum), 'chanSubSame_2Range'];
mkdir(pptTitleName);
[opt,ppt] = figPPT(pptTitleName);

plotSNRRange = [-3 3; -1  1];
freStr = {'8hz','20hz','mean of 8+20'};
timeBinStr = {'201To1200'};

for iSub = 1:size(allFrePointSNRMinusBase_4con_2Range,3)+1
				
	figure;
    set (gcf,'Position',[50,50,1800,900], 'color','w');

	for iRange = 1:size(SNR_8hz20hz_2Range,5)		

		for iFre = 1:size(SNR_8hz20hz_2Range,1)  % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-TimeRange
		
			subplot(3,3,3*(iRange-1)+iFre);

			cSNR = squeeze(SNR_8hz20hz_2Range(iFre,:,iSub,:,iRange)); 
			% cError = SNR_8hz20hz_2Range(iFre,:,iSub+1,:,iRange);			
			H=bar(cSNR,'grouped','EdgeColor','none');
		    ch = get(H,'children');
	        set(ch{1},'Facecolor',barColor(1,:));
	        set(ch{2},'Facecolor',barColor(2,:));
	        set(ch{3},'Facecolor',barColor(3,:));
	        set(ch{4},'Facecolor',barColor(4,:));
			hold on;box off;
			ylim(plotSNRRange(1,:));

			if iSub==size(allFrePointSNRMinusBase_4con_2Range,3)+1
				subStr = 'allSubMean';
				cError = squeeze(errorBar_SNR_8hz20hz_2Range(iFre,:,:,iRange));
				for i = 1:nbars
				    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
				    errorbar(x,cSNR(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
				    hold on;
					ylim(plotSNRRange(2,:));
				end

			elseif iSub<=size(allFrePointSNRMinusBase_4con_2Range,3)
					subStr = ['Sub:',num2str(iSub)];
			end

			title([subStr,freStr(iFre),timeBinStr(iRange)]);
			clear cSNR cError;
		end
	end
	l = legend([ch{1} ch{2} ch{3} ch{4} ],'sepCue--Uncue','ajCue--Uncue','ajCue--Cue','SepCue--Cue','Location',[0.908333333333333 0.248888888888889 0.0816666666666667 0.0881481481481481]);
	legend('boxoff');
	cfigFileName = fullfile(fullfile(pwd,pptTitleName),[pptTitleName,'_',subStr([1:3,5:end]),'.bmp']);
	saveas(gcf,cfigFileName);
	ppt = ppt.addImageSlide(['SNR_',subStr(1:3),'_',subStr(5:end)],cfigFileName,opt);
	close all;
end

ppt.saveAs(fullfile(pwd,pptTitleName));

ppt.close;




%%%%%%   bar of Power %%%%%%%%
pptTitleName = ['Power_max', num2str(maxChanNum), 'chanSubSame_2Range'];
mkdir(pptTitleName);
[opt,ppt] = figPPT(pptTitleName);

plotPowerRange = [-0.3 0.3; -0.1  0.1];
freStr = {'8hz','20hz','mean of 8+20'};
timeBinStr = {'201To1200'};

for iSub = 1:size(allFrePointPowerMinusBase_4con_2Range,3)+1
				
	figure;
    set (gcf,'Position',[50,50,1800,900], 'color','w');

	for iRange = 1:size(Power_8hz20hz_2Range,5)		

		for iFre = 1:size(Power_8hz20hz_2Range,1)  % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-TimeRange
		
			subplot(3,3,3*(iRange-1)+iFre);

			cPower = squeeze(Power_8hz20hz_2Range(iFre,:,iSub,:,iRange)); 
			% cError = Power_8hz20hz_2Range(iFre,:,iSub+1,:,iRange);			
			H=bar(cPower,'grouped','EdgeColor','none');
		    ch = get(H,'children');
	        set(ch{1},'Facecolor',barColor(1,:));
	        set(ch{2},'Facecolor',barColor(2,:));
	        set(ch{3},'Facecolor',barColor(3,:));
	        set(ch{4},'Facecolor',barColor(4,:));
			hold on;box off;
			ylim(plotPowerRange(1,:));

			if iSub==size(allFrePointPowerMinusBase_4con_2Range,3)+1
				subStr = 'allSubMean';
				cError = squeeze(errorBar_Power_8hz20hz_2Range(iFre,:,:,iRange));
				for i = 1:nbars
				    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
				    errorbar(x,cPower(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
				    hold on;
					ylim(plotPowerRange(2,:));
				end

			elseif iSub<=size(allFrePointPowerMinusBase_4con_2Range,3)
					subStr = ['Sub:',num2str(iSub)];
			end

			title([subStr,freStr(iFre),timeBinStr(iRange)]);
			clear cPower cError;
		end
	end
	l = legend([ch{1} ch{2} ch{3} ch{4} ],'sepCue--Uncue','ajCue--Uncue','ajCue--Cue','SepCue--Cue','Location',[0.908333333333333 0.248888888888889 0.0816666666666667 0.0881481481481481]);
	legend('boxoff');
	cfigFileName = fullfile(fullfile(pwd,pptTitleName),[pptTitleName,'_',subStr([1:3,5:end]),'.bmp']);
	saveas(gcf,cfigFileName);
	ppt = ppt.addImageSlide(['Power_',subStr(1:3),'_',subStr(5:end)],cfigFileName,opt);
	close all;
end

ppt.saveAs(fullfile(pwd,pptTitleName));

ppt.close;









%%%%%%%%% 2Range800  %%%%%%%%%%
SNR_8hz20hz_2Range800   = [];
Power_8hz20hz_2Range800 = [];

for iFre = 1:size(allFrePointSNRMinusBase_4con_2Range800,1)  % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-TimeRange

	for iSub = 1:size(allFrePointSNRMinusBase_4con_2Range800,3) 
			cSNR = allFrePointSNRMinusBase_4con_2Range800(iFre,chooseChanSNR_subSame(iFre,:),iSub,:,:); 
			SNR_8hz20hz_2Range800(iFre,:,iSub,:,:) = cSNR;

			cPower = allFrePointPowerMinusBase_4con_2Range800(iFre,chooseChanPower_subSame(iFre,:),iSub,:,:); 
			Power_8hz20hz_2Range800(iFre,:,iSub,:,:) = cPower;

	end
end

SNR_8hz20hz_2Range800(3,:,:,:,:)                                   = mean(SNR_8hz20hz_2Range800,1);
SNR_8hz20hz_2Range800(:,size(chooseChanSNR_subSame,2)+1,:,:,:)     = mean(SNR_8hz20hz_2Range800,2);
SNR_8hz20hz_2Range800(:,:,22,:,:)                                  = mean(SNR_8hz20hz_2Range800,3);

Power_8hz20hz_2Range800(3,:,:,:,:)                                 = mean(Power_8hz20hz_2Range800,1);
Power_8hz20hz_2Range800(:,size(chooseChanPower_subSame,2)+1,:,:,:) = mean(Power_8hz20hz_2Range800,2);
Power_8hz20hz_2Range800(:,:,22,:,:)                                = mean(Power_8hz20hz_2Range800,3);



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
		end
	end
end

maxChanNum = 5;
%%-------------- start to plot -------------------/
ngroups    = size(SNR_8hz20hz_2Range800,2);
nbars      = size(SNR_8hz20hz_2Range800,4);
groupwidth = min(0.8, nbars/(nbars+1.5));	
barColor   = [235 95 52; 254 158 15;63 96 169; 0 0 200]/255;



%%%%%%   bar of SNR %%%%%%%%
pptTitleName = ['SNR_max', num2str(maxChanNum), 'chanSubSame_2Range800'];
mkdir(pptTitleName);
[opt,ppt] = figPPT(pptTitleName);

plotSNRRange = [-3 3; -1  1];
freStr = {'8hz','20hz','mean of 8+20'};
timeBinStr = {'401To1200'};

for iSub = 1:size(allFrePointSNRMinusBase_4con_2Range800,3)+1
				
	figure;
    set (gcf,'Position',[50,50,1800,900], 'color','w');

	for iRange = 1:size(SNR_8hz20hz_2Range800,5)		

		for iFre = 1:size(SNR_8hz20hz_2Range800,1)  % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-TimeRange
		
			subplot(3,3,3*(iRange-1)+iFre);

			cSNR = squeeze(SNR_8hz20hz_2Range800(iFre,:,iSub,:,iRange)); 
			% cError = SNR_8hz20hz_2Range800(iFre,:,iSub+1,:,iRange);			
			H=bar(cSNR,'grouped','EdgeColor','none');
		    ch = get(H,'children');
	        set(ch{1},'Facecolor',barColor(1,:));
	        set(ch{2},'Facecolor',barColor(2,:));
	        set(ch{3},'Facecolor',barColor(3,:));
	        set(ch{4},'Facecolor',barColor(4,:));
			hold on;box off;
			ylim(plotSNRRange(1,:));

			if iSub==size(allFrePointSNRMinusBase_4con_2Range800,3)+1
				subStr = 'allSubMean';
				cError = squeeze(errorBar_SNR_8hz20hz_2Range800(iFre,:,:,iRange));
				for i = 1:nbars
				    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
				    errorbar(x,cSNR(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
				    hold on;
					ylim(plotSNRRange(2,:));
				end

			elseif iSub<=size(allFrePointSNRMinusBase_4con_2Range800,3)
					subStr = ['Sub:',num2str(iSub)];
			end

			title([subStr,freStr(iFre),timeBinStr(iRange)]);
			clear cSNR cError;
		end
	end
	l = legend([ch{1} ch{2} ch{3} ch{4} ],'sepCue--Uncue','ajCue--Uncue','ajCue--Cue','SepCue--Cue','Location',[0.908333333333333 0.248888888888889 0.0816666666666667 0.0881481481481481]);
	legend('boxoff');
	cfigFileName = fullfile(fullfile(pwd,pptTitleName),[pptTitleName,'_',subStr([1:3,5:end]),'.bmp']);
	saveas(gcf,cfigFileName);
	ppt = ppt.addImageSlide(['SNR_',subStr(1:3),'_',subStr(5:end)],cfigFileName,opt);
	close all;
end

ppt.saveAs(fullfile(pwd,pptTitleName));

ppt.close;




%%%%%%   bar of Power %%%%%%%%
pptTitleName = ['Power_max', num2str(maxChanNum), 'chanSubSame_2Range800'];
mkdir(pptTitleName);
[opt,ppt] = figPPT(pptTitleName);

plotPowerRange = [-0.3 0.3; -0.1  0.1];
freStr = {'8hz','20hz','mean of 8+20'};
timeBinStr = {'401To1200'};

for iSub = 1:size(allFrePointPowerMinusBase_4con_2Range800,3)+1
				
	figure;
    set (gcf,'Position',[50,50,1800,900], 'color','w');

	for iRange = 1:size(Power_8hz20hz_2Range800,5)		

		for iFre = 1:size(Power_8hz20hz_2Range800,1)  % size: 1-freBins; 2-chan; 3-sub; 4-condition;5-TimeRange
		
			subplot(3,3,3*(iRange-1)+iFre);

			cPower = squeeze(Power_8hz20hz_2Range800(iFre,:,iSub,:,iRange)); 
			% cError = Power_8hz20hz_2Range800(iFre,:,iSub+1,:,iRange);			
			H=bar(cPower,'grouped','EdgeColor','none');
		    ch = get(H,'children');
	        set(ch{1},'Facecolor',barColor(1,:));
	        set(ch{2},'Facecolor',barColor(2,:));
	        set(ch{3},'Facecolor',barColor(3,:));
	        set(ch{4},'Facecolor',barColor(4,:));
			hold on;box off;
			ylim(plotPowerRange(1,:));

			if iSub==size(allFrePointPowerMinusBase_4con_2Range800,3)+1
				subStr = 'allSubMean';
				cError = squeeze(errorBar_Power_8hz20hz_2Range800(iFre,:,:,iRange));
				for i = 1:nbars
				    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
				    errorbar(x,cPower(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
				    hold on;
					ylim(plotPowerRange(2,:));
				end

			elseif iSub<=size(allFrePointPowerMinusBase_4con_2Range800,3)
					subStr = ['Sub:',num2str(iSub)];
			end

			title([subStr,freStr(iFre),timeBinStr(iRange)]);
			clear cPower cError;
		end
	end
	l = legend([ch{1} ch{2} ch{3} ch{4} ],'sepCue--Uncue','ajCue--Uncue','ajCue--Cue','SepCue--Cue','Location',[0.908333333333333 0.248888888888889 0.0816666666666667 0.0881481481481481]);
	legend('boxoff');
	cfigFileName = fullfile(fullfile(pwd,pptTitleName),[pptTitleName,'_',subStr([1:3,5:end]),'.bmp']);
	saveas(gcf,cfigFileName);
	ppt = ppt.addImageSlide(['Power_',subStr(1:3),'_',subStr(5:end)],cfigFileName,opt);
	close all;
end

ppt.saveAs(fullfile(pwd,pptTitleName));

ppt.close;












clear SNR_merge820_3Range pvalue_SNR_merge820_3Range
clear SNR_merge820_2Range pvalue_SNR_merge820_2Range
clear SNR_merge820_2Range800 pvalue_SNR_merge820_2Range800


SNR_merge820_3Range = SNR_8hz20hz_3Range(3,:,[1:21],:,:);

for iTimeRange = 1:size(SNR_merge820_3Range,5)
 	
 	for  iChan = 1:size(SNR_merge820_3Range,2)

 		for iCmp = 1:3

 			[H,pvalue_SNR_merge820_3Range(iCmp,iChan,iTimeRange),CI] = ttest(squeeze(SNR_merge820_3Range(1,iChan,:,1,iTimeRange)), squeeze(SNR_merge820_3Range(1,iChan,:,iCmp+1,iTimeRange)));

 		end
 	end
end




SNR_merge820_2Range = SNR_8hz20hz_2Range(3,:,[1:21],:,:);

for iTimeRange = 1:size(SNR_merge820_2Range,5)
 	
 	for  iChan = 1:size(SNR_merge820_2Range,2)

 		for iCmp = 1:3

 			[H,pvalue_SNR_merge820_2Range(iCmp,iChan,iTimeRange),CI] = ttest(squeeze(SNR_merge820_2Range(1,iChan,:,1,iTimeRange)), squeeze(SNR_merge820_2Range(1,iChan,:,iCmp+1,iTimeRange)));

 		end
 	end
end


SNR_merge820_2Range800 = SNR_8hz20hz_2Range800(3,:,[1:21],:,:);

for iTimeRange = 1:size(SNR_merge820_2Range800,5)
 	
 	for  iChan = 1:size(SNR_merge820_2Range800,2)

 		for iCmp = 1:3

 			[H,pvalue_SNR_merge820_2Range800(iCmp,iChan,iTimeRange),CI] = ttest(squeeze(SNR_merge820_2Range800(1,iChan,:,1,iTimeRange)), squeeze(SNR_merge820_2Range800(1,iChan,:,iCmp+1,iTimeRange)));

 		end
 	end
end





clear Power_merge820_3Range pvalue_Power_merge820_3Range
clear Power_merge820_2Range pvalue_Power_merge820_2Range
clear Power_merge820_2Range800 pvalue_Power_merge820_2Range800


Power_merge820_3Range = Power_8hz20hz_3Range(3,:,[1:21],:,:);

for iTimeRange = 1:size(Power_merge820_3Range,5)
 	
 	for  iChan = 1:size(Power_merge820_3Range,2)

 		for iCmp = 1:3

 			[H,pvalue_Power_merge820_3Range(iCmp,iChan,iTimeRange),CI] = ttest(squeeze(Power_merge820_3Range(1,iChan,:,1,iTimeRange)), squeeze(Power_merge820_3Range(1,iChan,:,iCmp+1,iTimeRange)));

 		end
 	end
end




Power_merge820_2Range = Power_8hz20hz_2Range(3,:,[1:21],:,:);

for iTimeRange = 1:size(Power_merge820_2Range,5)
 	
 	for  iChan = 1:size(Power_merge820_2Range,2)

 		for iCmp = 1:3

 			[H,pvalue_Power_merge820_2Range(iCmp,iChan,iTimeRange),CI] = ttest(squeeze(Power_merge820_2Range(1,iChan,:,1,iTimeRange)), squeeze(Power_merge820_2Range(1,iChan,:,iCmp+1,iTimeRange)));

 		end
 	end
end


Power_merge820_2Range800 = Power_8hz20hz_2Range800(3,:,[1:21],:,:);

for iTimeRange = 1:size(Power_merge820_2Range800,5)
 	
 	for  iChan = 1:size(Power_merge820_2Range800,2)

 		for iCmp = 1:3

 			[H,pvalue_Power_merge820_2Range800(iCmp,iChan,iTimeRange),CI] = ttest(squeeze(Power_merge820_2Range800(1,iChan,:,1,iTimeRange)), squeeze(Power_merge820_2Range800(1,iChan,:,iCmp+1,iTimeRange)));

 		end
 	end
end

