


% directory = 'D:\expData\SSVEP_8loc\Exp1_20170328\Data\scanData\20180201\cnt';
% allFile   = dir(fullfile(directory,'*BA3.eeg'));
% fftType   = 'avgFFT'; % 'trialFFT'
% cd(directory);
% freRange  = [5,25];
% targetFre = [8;20]; % for loc2, loc3
% [P1eNo,standard_eName] = getElectrodeNo_bcl(1,[],[]);


% %---- do fft of [-500 1200] for all conditions, find the chan ----------/
% timeBin      = [-500 1200];
% stimTrigger  = [1:40];% 

% for iSub = 1:size(allFile,1)

% 	eegName = allFile(iSub).name;
% 	[iSubAvg,freBins,Power,SNR] = eegfftSNR_bcl(eegName,freRange,timeBin,stimTrigger,fftType); 
% 	allPower(:,:,iSub)          = Power; % size: 1-freBins; 2-chan; 3-sub; 4-condition; 5-timeRange(base)
% 	allSNR(:,:,iSub)            = SNR;
% 	clear  Power SNR iSubAvg;

% end 

% %-------find closest frePoint-------/
% for iFre = 1:size(targetFre,1)
% 	[X, freIndex]         = sort(abs(freBins-targetFre(iFre,1)),'ascend');
% 	freTimePoint(iFre,1)  = freIndex(1);
% 	realFre(iFre,1) 	  = freBins(freIndex(1));
% end
% %-----------------------------------\
% frePointSNR   = allSNR(freTimePoint(:),:,:);
% frePointPower = allPower(freTimePoint(:),:,:);

% for iSub = 1:size(frePointSNR,3)
% 	for iFre = 1:size(targetFre,1)
				
% 		[X, chanIndex]               = sort(frePointSNR(iFre,[1:size(standard_eName,2)-2],iSub),'descend');						
% 		chooseChanSNR(iSub,:,iFre)   = chanIndex(1:3);	
% 		clear X chanIndex
				
% 		[X, chanIndex]               = sort(frePointPower(iFre,[1:size(standard_eName,2)-2],iSub),'descend');						
% 		chooseChanPower(iSub,:,iFre) = chanIndex(1:3);	
% 		clear X  chanIndex
% 	end
% end

% save('mergeAllcon_-500To1200');

% %----------------  have found the chan, it's the end  -------------\




% %------plot topo of SNR to check chosen chan ---------/
% load('mergeAllcon_-500To1200');
% waveNum = 3;
% [myColormap] = makeColormap(20,0,1,'r');

% subNumPerFig = 5;
% figNum       = ceil(size(allFile,1)./subNumPerFig);
% [P1eNo,standard_eName]=getElectrodeNo_bcl(0,'haveM1',[]);

% for iPage = 1:figNum
% 	figure;
% 	set(gcf,'Position',get(0,'ScreenSize'),'color','w');

% 	suptitle([fftType,'SNR     -500To1200    chooseChan   left:    ', num2str(realFre(1)) 'Hz          right:    ',num2str(realFre(2))  'Hz']);
	
% 	cSub = [1:subNumPerFig]+(iPage-1)*subNumPerFig;

% 	for iSub = cSub
				

% 		for istim = 1:size(targetFre,1)
				
% 				if iSub<=size(chooseChanSNR,1)+1

% 					if iSub==size(chooseChanSNR,1)+1

% 						MaxChanSNR_Power(:,:,size(chooseChanSNR,1)+1,istim) = mean(MaxChanSNR_Power(:,:,:,istim),3);
% 						MaxChanSNR_SNR(:,:,size(chooseChanSNR,1)+1,istim)   = mean(MaxChanSNR_SNR(:,:,:,istim),3);
% 					end
					
% 					%-------plot topo-------/
% 					if istim==1
% 						figureLoc = 2*(waveNum+1)*(iSub-1-(iPage-1)*subNumPerFig)+1;
% 					elseif istim==2
% 						figureLoc = 2*(waveNum+1)*(iSub-1-(iPage-1)*subNumPerFig)+waveNum+2;
% 					end	

% 					subplot(subNumPerFig,2*(waveNum+1),figureLoc);
% 					if iSub<= size(chooseChanSNR,1)
% 						topoplot_bcl(frePointSNR(istim,:,iSub),'EEGChans',standard_eName,'maplimitsDouble',[0 3],'colormap',myColormap);
% 						chanIndex = chooseChanSNR(iSub,:,istim);
% 					end
% 					%-----------------------\	
					


% 					for iWave = 1:waveNum
% 						figureLoc  = figureLoc+1;
% 						subplot(subNumPerFig,2*(waveNum+1),figureLoc);

% 						if iSub<= size(chooseChanSNR,1)
% 							MaxChanSNR_Power(:,iWave,iSub,istim) =  allPower(:,chanIndex(iWave),iSub);
% 							MaxChanSNR_SNR(:,iWave,iSub,istim)   =  allSNR(:,chanIndex(iWave),iSub);
% 						end

% 						[AX,H1,H2] = plotyy(freBins,MaxChanSNR_Power(:,iWave,iSub,istim),freBins,MaxChanSNR_SNR(:,iWave,iSub,istim),'plot');
						
% 						if iSub<=size(chooseChanSNR,1)
% 							title([num2str(iSub),standard_eName(chanIndex(iWave))]);
% 						else
% 							title(['allSubAvg']);
% 						end
						
% 						box off;

% 						set(AX(1),'XColor','k','YColor','b','Ylim',[0,0.15],'ytick',[0:0.05:0.15],'Xlim',freRange);		 				
% 						set(AX(2),'XColor','k','YColor','r','Ylim',[0,5],'ytick',[0:1:5],'Xlim',freRange);
% 						set(H1,'LineStyle','-','color','b');
% 						set(H2,'LineStyle','-','color','r');
% 						line(repmat(realFre(istim,:),1,2),[0,.2],'linestyle',':','linewidth',2, 'color',[0 0 0 ]); %
% 					end
% 				end
% 		end
% 	end
% 	legend([H1,H2],{'power';'SNR'},'Location','North');
% 	legend('boxoff');
% 	set(gcf,'color','w');
% 	% saveas(gcf,['S',num2str(cSub(1)),'~',num2str(cSub(end)),'_',fftType,'_SNRtopo.pdf']);
% end
% %------plot topo of SNR to check chosen chan,it' the end----------\





% %------plot topo of Power to check chosen chan ---------/
% load('mergeAllcon_-500To1200');
% waveNum = 3;
% [myColormap] = makeColormap(20,0,1,'b');
% subNumPerFig = 5;
% figNum       = ceil(size(allFile,1)./subNumPerFig);
% [P1eNo,standard_eName]=getElectrodeNo_bcl(0,'haveM1',[]);

% for iPage = 1:figNum
% 	figure;
% 	set(gcf,'Position',get(0,'ScreenSize'),'color','w');

% 	suptitle([fftType,'Power    -500To1200    chooseChan   left:    ', num2str(realFre(1)) 'Hz          right:    ',num2str(realFre(2))  'Hz']);
	
% 	cSub = [1:subNumPerFig]+(iPage-1)*subNumPerFig;

% 	for iSub = cSub
				

% 		for istim = 1:size(targetFre,1)
				
% 				if iSub<=size(chooseChanPower,1)+1

% 					if iSub==size(chooseChanPower,1)+1

% 						MaxChanPower_Power(:,:,size(chooseChanPower,1)+1,istim) = mean(MaxChanPower_Power(:,:,:,istim),3);
% 						MaxChanPower_SNR(:,:,size(chooseChanPower,1)+1,istim)   = mean(MaxChanPower_SNR(:,:,:,istim),3);
% 					end
					
% 					%-------plot topo-------/
% 					if istim==1
% 						figureLoc = 2*(waveNum+1)*(iSub-1-(iPage-1)*subNumPerFig)+1;
% 					elseif istim==2
% 						figureLoc = 2*(waveNum+1)*(iSub-1-(iPage-1)*subNumPerFig)+waveNum+2;
% 					end	

% 					subplot(subNumPerFig,2*(waveNum+1),figureLoc);
% 					if iSub<= size(chooseChanPower,1)
% 						topoplot_bcl(frePointPower(istim,:,iSub),'EEGChans',standard_eName,'maplimitsDouble',[0 0.15],'colormap',myColormap);
% 						chanIndex = chooseChanPower(iSub,:,istim);
% 					end
% 					%-----------------------\	
					


% 					for iWave = 1:waveNum
% 						figureLoc  = figureLoc+1;
% 						subplot(subNumPerFig,2*(waveNum+1),figureLoc);

% 						if iSub<= size(chooseChanPower,1)
% 							MaxChanPower_Power(:,iWave,iSub,istim) =  allPower(:,chanIndex(iWave),iSub);
% 							MaxChanPower_SNR(:,iWave,iSub,istim)   =  allSNR(:,chanIndex(iWave),iSub);
% 						end

% 						[AX,H1,H2] = plotyy(freBins,MaxChanPower_Power(:,iWave,iSub,istim),freBins,MaxChanPower_SNR(:,iWave,iSub,istim),'plot');
						
% 						if iSub<=size(chooseChanPower,1)
% 							title([num2str(iSub),standard_eName(chanIndex(iWave))]);
% 						else
% 							title(['allSubAvg']);
% 						end
						
% 						box off;

% 						set(AX(1),'XColor','k','YColor','b','Ylim',[0,0.15],'ytick',[0:0.05:0.15],'Xlim',freRange);		 				
% 						set(AX(2),'XColor','k','YColor','r','Ylim',[0,5],'ytick',[0:1:5],'Xlim',freRange);
% 						set(H1,'LineStyle','-','color','b');
% 						set(H2,'LineStyle','-','color','r');
% 						line(repmat(realFre(istim,:),1,2),[0,.2],'linestyle',':','linewidth',2, 'color',[0 0 0 ]); %
% 					end
% 				end
% 		end
% 	end
% 	legend([H1,H2],{'Power';'SNR'},'Location','North');
% 	legend('boxoff');
% 	set(gcf,'color','w');
% 	% saveas(gcf,['S',num2str(cSub(1)),'~',num2str(cSub(end)),'_',fftType,'_SNRtopo.pdf']);
% end

% %------plot topo of Power to check chosen chan, the end------\







% %-------------------- do fft of baseline and post-baseline for 4 conditions,separately ------------------/
% timeBin      = [-500,-1;201,1200];
% stimTrigger  = [1:10;11:20;21:30;31:40];

% for iSub = 1:size(allFile,1)

% 	eegName = allFile(iSub).name;

% 	for iCondition = 1:size(stimTrigger,1)
	
% 		for iTimeRange = 1:size(timeBin,1)  % 
			
% 			[iSubAvg,freBins,Power,SNR] = eegfftSNR_bcl(eegName,freRange,timeBin(iTimeRange,:),stimTrigger(iCondition,:),fftType); 
% 			allPower_4con(:,:,iSub,iCondition,iTimeRange)  = Power; % size: 1-freBins; 2-chan; 3-sub; 4-condition; 5-timeRange(base)
% 			allSNR_4con(:,:,iSub,iCondition,iTimeRange)    = SNR;
% 			clear  Power SNR iSubAvg;
% 		end
% 	end
% end 

% allPowerMinusBase_4con         = squeeze(allPower_4con(:,:,:,:,2)-allPower_4con(:,:,:,:,1)); % size: 1-freBins; 2-chan; 3-sub; 4-condition;
% allFrePointPowerMinusBase_4con = allPowerMinusBase_4con(freTimePoint(:),:,:,:); % size: 1-freBins; 2-chan; 3-sub; 4-condition;
% allSNRMinusBase_4con           = squeeze(allSNR_4con(:,:,:,:,2)-allSNR_4con(:,:,:,:,1)); % size: 1-freBins; 2-chan; 3-sub; 4-condition;
% allFrePointSNRMinusBase_4con   = allSNRMinusBase_4con(freTimePoint(:),:,:,:); % size: 1-freBins; 2-chan; 3-sub; 4-condition;
% save('fft_4con_base&postBase');
% %---------- It's the end..  Do fft  for 4 conditions,separately. It's the end---------------\




% %----- do fft of baseline and post-baseline for 2 merged conditions(sepCue & ajCue),separately ------------------/
% timeBin      = [-500,-1;201,1200];
% stimTrigger  = [1:10,31:40;11:20,21:30;];% 1:sepCue; 2:ajCue.

% for iSub = 1:size(allFile,1)

% 	eegName = allFile(iSub).name;

% 	for iCondition = 1:size(stimTrigger,1)
	
% 		for iTimeRange = 1:size(timeBin,1)  % 
			
% 			[iSubAvg,freBins,Power,SNR] = eegfftSNR_bcl(eegName,freRange,timeBin(iTimeRange,:),stimTrigger(iCondition,:),fftType); 
% 			allPower_2con(:,:,iSub,iCondition,iTimeRange)  = Power; % size: 1-freBins; 2-chan; 3-sub; 4-condition; 5-timeRange(base)
% 			allSNR_2con(:,:,iSub,iCondition,iTimeRange)    = SNR;
% 			clear  Power SNR iSubAvg;
% 		end
% 	end
% end 

% allPowerMinusBase_2con         = squeeze(allPower_2con(:,:,:,:,2)-allPower_2con(:,:,:,:,1)); % size: 1-freBins; 2-chan; 3-sub; 4-condition;
% allFrePointPowerMinusBase_2con = allPowerMinusBase_2con(freTimePoint(:),:,:,:); % size: 1-freBins; 2-chan; 3-sub; 4-condition;
% allSNRMinusBase_2con           = squeeze(allSNR_2con(:,:,:,:,2)-allSNR_2con(:,:,:,:,1)); % size: 1-freBins; 2-chan; 3-sub; 4-condition;
% allFrePointSNRMinusBase_2con   = allSNRMinusBase_2con(freTimePoint(:),:,:,:); % size: 1-freBins; 2-chan; 3-sub; 4-condition;
% save('fft_2con_base&postBase');
% %-------------It's the end. Do fft  for 2 merged conditions(sepCue & ajCue),separately. It's the end -----------\


clear all;
load('mergeAllcon_-500To1200');
load('fft_2con_base&postBase.mat');
load('fft_4con_base&postBase');
[P1eNo,standard_eName] = getElectrodeNo_bcl(1,[],[]);
subNumPerFig    = 5;
subfigNumPerRow = 9;
figNum          = ceil(size(allPowerMinusBase_4con,3)./subNumPerFig);

%%%%%%   bar of SNR %%%%%%%%

plotSNRRange = [-2 3];
for iPage = 1:figNum
	figure;
	set(gcf,'Position',get(0,'ScreenSize'),'color','w');
	set(gcf,'color','white')
	
	cSub = [1:subNumPerFig]+(iPage-1)*subNumPerFig;
	suptitle('SNR');
	for iSub = cSub
		if iSub<=size(allSNRMinusBase_4con,3)+1

			if iSub==size(allSNRMinusBase_4con,3)+1 % compute mean 

				for iChan = 1:size(SNR_con13_8hz,1)
  						[gm,gsem] = grpstats(squeeze(SNR_con13_8hz(iChan,:,[1:size(allSNRMinusBase_4con,3)]))',{},{'mean','sem'});
  						SNR_con13_8hz(iChan,:,size(allSNRMinusBase_4con,3)+1) = gm;
  						SNR_con13_8hz(iChan,:,size(allSNRMinusBase_4con,3)+2) = gsem;

  						
  						[gm,gsem] = grpstats(squeeze(SNR_con42_20hz(iChan,:,[1:size(allSNRMinusBase_4con,3)]))',{},{'mean','sem'});
  						SNR_con42_20hz(iChan,:,size(allSNRMinusBase_4con,3)+1) = gm;
  						SNR_con42_20hz(iChan,:,size(allSNRMinusBase_4con,3)+2) = gsem;
  						
  						[gm,gsem] = grpstats(squeeze(SNR_con13_8hz_con42_20hz(iChan,:,[1:size(allSNRMinusBase_4con,3)]))',{},{'mean','sem'});
  						SNR_con13_8hz_con42_20hz(iChan,:,size(allSNRMinusBase_4con,3)+1) = gm;
  						SNR_con13_8hz_con42_20hz(iChan,:,size(allSNRMinusBase_4con,3)+2) = gsem;

  						[gm,gsem] = grpstats(squeeze(SNR_all_8hz(iChan,:,[1:size(allSNRMinusBase_4con,3)]))',{},{'mean','sem'});
  						SNR_all_8hz(iChan,:,size(allSNRMinusBase_4con,3)+1) = gm;
  						SNR_all_8hz(iChan,:,size(allSNRMinusBase_4con,3)+2) = gsem;

  						[gm,gsem] = grpstats(squeeze(SNR_all_20hz(iChan,:,[1:size(allSNRMinusBase_4con,3)]))',{},{'mean','sem'});
  						SNR_all_20hz(iChan,:,size(allSNRMinusBase_4con,3)+1) = gm;
  						SNR_all_20hz(iChan,:,size(allSNRMinusBase_4con,3)+2) = gsem;

  						[gm,gsem] = grpstats(squeeze(SNR_all(iChan,:,[1:size(allSNRMinusBase_4con,3)]))',{},{'mean','sem'});
  						SNR_all(iChan,:,size(allSNRMinusBase_4con,3)+1) = gm;
  						SNR_all(iChan,:,size(allSNRMinusBase_4con,3)+2) = gsem;

  						[gm,gsem] = grpstats(squeeze(SNR_8hz_2con(iChan,:,[1:size(allSNRMinusBase_4con,3)]))',{},{'mean','sem'});
  						SNR_8hz_2con(iChan,:,size(allSNRMinusBase_4con,3)+1) = gm;
  						SNR_8hz_2con(iChan,:,size(allSNRMinusBase_4con,3)+2) = gsem;  						  						

  						[gm,gsem] = grpstats(squeeze(SNR_20hz_2con(iChan,:,[1:size(allSNRMinusBase_4con,3)]))',{},{'mean','sem'});
  						SNR_20hz_2con(iChan,:,size(allSNRMinusBase_4con,3)+1) = gm;
  						SNR_20hz_2con(iChan,:,size(allSNRMinusBase_4con,3)+2) = gsem;

  						[gm,gsem] = grpstats(squeeze(SNR_all_2con(iChan,:,[1:size(allSNRMinusBase_4con,3)]))',{},{'mean','sem'});
  						SNR_all_2con(iChan,:,size(allSNRMinusBase_4con,3)+1) = gm;
  						SNR_all_2con(iChan,:,size(allSNRMinusBase_4con,3)+2) = gsem;  	

						ngroups    = 3;
						nbars      = 2;
						groupwidth = min(0.8, nbars/(nbars+1.5));	

  				end

			end
			

			%----- 1: Cue13(8hz) Vs. Cue34(8hz) ------------
			% size: 1-freBins; 2-chan; 3-sub; 4-condition;
			subplot(subNumPerFig,subfigNumPerRow,subfigNumPerRow*(iSub-1-(iPage-1)*subNumPerFig)+1);
			if iSub<=size(allSNRMinusBase_4con,3)
					cSNR = squeeze(allFrePointSNRMinusBase_4con(1,chooseChanSNR(iSub,:,1),iSub,[1 3])); 
					SNR_con13_8hz(:,:,iSub) = cSNR;
			else
					cSNR = SNR_con13_8hz(:,:,iSub);
					cError = SNR_con13_8hz(:,:,iSub+1);
			end
			H=bar(cSNR,'BarWidth',0.3);
			hold on;
			% axis([0.7 3.3 -0.2 0.2]);
			box off;
			ylim(plotSNRRange);

			if iSub==size(allSNRMinusBase_4con,3)+1
				subStr = 'allSubMean';
				for i = 1:nbars
				    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
				    errorbar(x,cSNR(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
				    hold on;
				end
			elseif iSub<=size(allSNRMinusBase_4con,3)
				subStr = ['Sub',num2str(iSub)];
				set(gca,'XTickLabel',standard_eName(chooseChanSNR(iSub,:,1)));
			end
			title([subStr,' cue13&cue34: 8hz']);
			clear cSNR cError;



			%----- 2: Cue42(20hz) Vs. Cue21(20hz) ------------
			% size: 1-freBins; 2-chan; 3-sub; 4-condition;
			subplot(subNumPerFig,subfigNumPerRow,subfigNumPerRow*(iSub-1-(iPage-1)*subNumPerFig)+2);
			if iSub<=size(allSNRMinusBase_4con,3)
					cSNR = squeeze(allFrePointSNRMinusBase_4con(2,chooseChanSNR(iSub,:,2),iSub,[4 2])); 
					SNR_con42_20hz(:,:,iSub) = cSNR;
			else
				    cSNR = SNR_con42_20hz(:,:,iSub); 
				    cError = SNR_con42_20hz(:,:,iSub+1); 				    
			end	
			H=bar(cSNR,'BarWidth',0.3);
			hold on;box off;
			ylim(plotSNRRange);
			if iSub<=size(allSNRMinusBase_4con,3)
					set(gca,'XTickLabel',standard_eName(chooseChanSNR(iSub,:,2)));
			else
					for i = 1:nbars
					    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
					    errorbar(x,cSNR(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
					    hold on;
					end	
			end
			title(['  cue42&cue21: 20hz']);
			clear cSNR  cError;



			%----- 3: Cue13(8hz)+Cue42(20hz) Vs. Cue34(8hz)+Cue21(20hz) ------------
			% size: 1-freBins; 2-chan; 3-sub; 4-condition;
			subplot(subNumPerFig,subfigNumPerRow,subfigNumPerRow*(iSub-1-(iPage-1)*subNumPerFig)+3);
			cSNR = mean(cat(3,SNR_con13_8hz(:,:,iSub),SNR_con42_20hz(:,:,iSub)),3);
			SNR_con13_8hz_con42_20hz(:,:,iSub) = cSNR;
			if iSub>size(allSNRMinusBase_4con,3)

					cError = SNR_con13_8hz_con42_20hz(:,:,iSub+1);
					for i = 1:nbars
					    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
					    errorbar(x,cSNR(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
					    hold on;
					end	
			end
			H=bar(cSNR,'BarWidth',0.3);
			ylim(plotSNRRange);
			hold on;box off;
			title([' mean of fig1&2']);
			clear cSNR  cError;



			%----- 4: Cue13(8hz)+Cue42(8hz) Vs. Cue34(8hz)+Cue21(8hz) ------------
			% size: 1-freBins; 2-chan; 3-sub; 4-condition;
			subplot(subNumPerFig,subfigNumPerRow,subfigNumPerRow*(iSub-1-(iPage-1)*subNumPerFig)+4);
			if iSub<=size(allSNRMinusBase_4con,3)
					cSNR(:,1) = squeeze(mean(allFrePointSNRMinusBase_4con(1,chooseChanSNR(iSub,:,1),iSub,[1 4]),4)); 
					cSNR(:,2) = squeeze(mean(allFrePointSNRMinusBase_4con(1,chooseChanSNR(iSub,:,1),iSub,[2 3]),4)); 
					SNR_all_8hz(:,:,iSub) = cSNR;
			else
					cSNR = SNR_all_8hz(:,:,iSub);
					cError = SNR_all_8hz(:,:,iSub+1);

			end	
			H=bar(cSNR,'BarWidth',0.3);
			ylim(plotSNRRange);
			hold on;box off;

			if iSub<=size(allSNRMinusBase_4con,3)
					set(gca,'XTickLabel',standard_eName(chooseChanSNR(iSub,:,1)));
			else
					for i = 1:nbars
					    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
					    errorbar(x,cSNR(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
					    hold on;
					end	
			end
			title(['  sepCue&ajCue: 8hz']);
			clear cSNR  cError;



			%----- 5: Cue13(20hz)+Cue42(20hz) Vs. Cue34(20hz)+Cue21(20hz) ------------
			% size: 1-freBins; 2-chan; 3-sub; 4-condition;
			subplot(subNumPerFig,subfigNumPerRow,subfigNumPerRow*(iSub-1-(iPage-1)*subNumPerFig)+5);
			if iSub<=size(allSNRMinusBase_4con,3)
					cSNR(:,1) = squeeze(mean(allFrePointSNRMinusBase_4con(2,chooseChanSNR(iSub,:,2),iSub,[1 4]),4)); 
					cSNR(:,2) = squeeze(mean(allFrePointSNRMinusBase_4con(2,chooseChanSNR(iSub,:,2),iSub,[2 3]),4)); 
					SNR_all_20hz(:,:,iSub) = cSNR;
			else
					cSNR = SNR_all_20hz(:,:,iSub);
					cError = SNR_all_20hz(:,:,iSub+1);
			end
			H=bar(cSNR,'BarWidth',0.3);
			ylim(plotSNRRange);
			hold on;box off;
			if iSub<=size(allSNRMinusBase_4con,3)
					set(gca,'XTickLabel',standard_eName(chooseChanSNR(iSub,:,2)));
			else
					for i = 1:nbars
					    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
					    errorbar(x,cSNR(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
					    hold on;
					end				
			end
			title(['  sepCue&ajCue: 20hz']);
			clear cSNR  cError;



			%----- 6: Cue13(20hz)+Cue42(8hz) Vs. Cue34(8hz)+Cue21(8hz) ------------
			% size: 1-freBins; 2-chan; 3-sub; 4-condition;
			subplot(subNumPerFig,subfigNumPerRow,subfigNumPerRow*(iSub-1-(iPage-1)*subNumPerFig)+6);
			cSNR = mean(cat(3,SNR_all_8hz(:,:,iSub),SNR_all_20hz(:,:,iSub)),3);
			SNR_all(:,:,iSub) = cSNR;
			
			if iSub>size(allSNRMinusBase_4con,3)

					cError = SNR_all(:,:,iSub+1);
					for i = 1:nbars
					    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
					    errorbar(x,cSNR(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
					    hold on;
					end	
			end

			H=bar(cSNR,'BarWidth',0.3);
			ylim(plotSNRRange);
			hold on;box off;
			title(['   mean of fig4&5']);
			clear cSNR  cError;



			%----- 7: mergeCondition(to imporove SNR),then do fft...SepCue(8hz) Vs. ajCue(8hz) ------------
			% size: 1-freBins; 2-chan; 3-sub; 4-condition;
			subplot(subNumPerFig,subfigNumPerRow,subfigNumPerRow*(iSub-1-(iPage-1)*subNumPerFig)+7);
			if iSub<=size(allSNRMinusBase_2con,3)
					cSNR = squeeze(allFrePointSNRMinusBase_2con(1,chooseChanSNR(iSub,:,1),iSub,:)); 
					SNR_8hz_2con(:,:,iSub) = cSNR;
			else
					cSNR = SNR_8hz_2con(:,:,iSub);
					cError = SNR_8hz_2con(:,:,iSub+1);
			end	
			H=bar(cSNR,'BarWidth',0.3);
			ylim(plotSNRRange);
			hold on;box off;

			if iSub<=size(allSNRMinusBase_2con,3)
					set(gca,'XTickLabel',standard_eName(chooseChanSNR(iSub,:,1)));
			else
					for i = 1:nbars
					    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
					    errorbar(x,cSNR(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
					    hold on;
					end				
			end
			title([' sepCue&ajCue fft: 8hz']);
			clear cSNR  cError;



			%----- 8: mergeCondition(to imporove SNR),then do fft...SepCue(20hz) Vs. ajCue(20hz) ------------
			% size: 1-freBins; 2-chan; 3-sub; 4-condition;
			subplot(subNumPerFig,subfigNumPerRow,subfigNumPerRow*(iSub-1-(iPage-1)*subNumPerFig)+8);
			if iSub<=size(allSNRMinusBase_2con,3)
					cSNR = squeeze(allFrePointSNRMinusBase_2con(2,chooseChanSNR(iSub,:,2),iSub,:)); 
					SNR_20hz_2con(:,:,iSub) = cSNR;
			else
					cSNR = SNR_20hz_2con(:,:,iSub);
					cError = SNR_20hz_2con(:,:,iSub+1);
			end	

			H=bar(cSNR,'BarWidth',0.3);
			ylim(plotSNRRange);
			hold on;box off;
			if iSub<=size(allSNRMinusBase_2con,3)
					set(gca,'XTickLabel',standard_eName(chooseChanSNR(iSub,:,2)));
			else
					for i = 1:nbars
					    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
					    errorbar(x,cSNR(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
					    hold on;
					end				
			end
			title(['  sepCue&ajCue fft: 20hz']);
			clear cSNR  cError;



			%----- 9: mergeCondition(to imporove SNR),then do fft...SepCue(8+20hz) Vs. ajCue(8+20hz) ------------
			% size: 1-freBins; 2-chan; 3-sub; 4-condition;
			subplot(subNumPerFig,subfigNumPerRow,subfigNumPerRow*(iSub-1-(iPage-1)*subNumPerFig)+9);
			cSNR = mean(cat(3,SNR_8hz_2con(:,:,iSub),SNR_20hz_2con(:,:,iSub)),3);
			SNR_all_2con(:,:,iSub) = cSNR;
			
			if iSub>size(allSNRMinusBase_4con,3)

					cError = SNR_all_2con(:,:,iSub+1);
					for i = 1:nbars
					    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
					    errorbar(x,cSNR(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
					    hold on;
					end	
			end			
			H=bar(cSNR,'BarWidth',0.3);
			ylim(plotSNRRange);
			hold on;box off;
			title(['   mean of fig7&8']);
			clear cSNR  cError;
		end
	end
	legend('sepCue','ajCue','Location','North');
	legend('boxoff');
end







%%%%%%   plot bar of Power %%%%%%%%

plotPowerRange = [-0.2 0.1];
for iPage = 1:figNum
	figure;
	set(gcf,'Position',get(0,'ScreenSize'),'color','w');
	set(gcf,'color','white')
	suptitle('Power');
	cSub = [1:subNumPerFig]+(iPage-1)*subNumPerFig;

	for iSub = cSub
		if iSub<=size(allPowerMinusBase_4con,3)+1

			if iSub==size(allPowerMinusBase_4con,3)+1 % compute mean 

				for iChan = 1:size(Power_con13_8hz,1)
  						[gm,gsem] = grpstats(squeeze(Power_con13_8hz(iChan,:,[1:size(allPowerMinusBase_4con,3)]))',{},{'mean','sem'});
  						Power_con13_8hz(iChan,:,size(allPowerMinusBase_4con,3)+1) = gm;
  						Power_con13_8hz(iChan,:,size(allPowerMinusBase_4con,3)+2) = gsem;

  						
  						[gm,gsem] = grpstats(squeeze(Power_con42_20hz(iChan,:,[1:size(allPowerMinusBase_4con,3)]))',{},{'mean','sem'});
  						Power_con42_20hz(iChan,:,size(allPowerMinusBase_4con,3)+1) = gm;
  						Power_con42_20hz(iChan,:,size(allPowerMinusBase_4con,3)+2) = gsem;
  						
  						[gm,gsem] = grpstats(squeeze(Power_con13_8hz_con42_20hz(iChan,:,[1:size(allPowerMinusBase_4con,3)]))',{},{'mean','sem'});
  						Power_con13_8hz_con42_20hz(iChan,:,size(allPowerMinusBase_4con,3)+1) = gm;
  						Power_con13_8hz_con42_20hz(iChan,:,size(allPowerMinusBase_4con,3)+2) = gsem;

  						[gm,gsem] = grpstats(squeeze(Power_all_8hz(iChan,:,[1:size(allPowerMinusBase_4con,3)]))',{},{'mean','sem'});
  						Power_all_8hz(iChan,:,size(allPowerMinusBase_4con,3)+1) = gm;
  						Power_all_8hz(iChan,:,size(allPowerMinusBase_4con,3)+2) = gsem;

  						[gm,gsem] = grpstats(squeeze(Power_all_20hz(iChan,:,[1:size(allPowerMinusBase_4con,3)]))',{},{'mean','sem'});
  						Power_all_20hz(iChan,:,size(allPowerMinusBase_4con,3)+1) = gm;
  						Power_all_20hz(iChan,:,size(allPowerMinusBase_4con,3)+2) = gsem;

  						[gm,gsem] = grpstats(squeeze(Power_all(iChan,:,[1:size(allPowerMinusBase_4con,3)]))',{},{'mean','sem'});
  						Power_all(iChan,:,size(allPowerMinusBase_4con,3)+1) = gm;
  						Power_all(iChan,:,size(allPowerMinusBase_4con,3)+2) = gsem;

  						[gm,gsem] = grpstats(squeeze(Power_8hz_2con(iChan,:,[1:size(allPowerMinusBase_4con,3)]))',{},{'mean','sem'});
  						Power_8hz_2con(iChan,:,size(allPowerMinusBase_4con,3)+1) = gm;
  						Power_8hz_2con(iChan,:,size(allPowerMinusBase_4con,3)+2) = gsem;  						  						

  						[gm,gsem] = grpstats(squeeze(Power_20hz_2con(iChan,:,[1:size(allPowerMinusBase_4con,3)]))',{},{'mean','sem'});
  						Power_20hz_2con(iChan,:,size(allPowerMinusBase_4con,3)+1) = gm;
  						Power_20hz_2con(iChan,:,size(allPowerMinusBase_4con,3)+2) = gsem;

  						[gm,gsem] = grpstats(squeeze(Power_all_2con(iChan,:,[1:size(allPowerMinusBase_4con,3)]))',{},{'mean','sem'});
  						Power_all_2con(iChan,:,size(allPowerMinusBase_4con,3)+1) = gm;
  						Power_all_2con(iChan,:,size(allPowerMinusBase_4con,3)+2) = gsem;  	

						ngroups    = 3;
						nbars      = 2;
						groupwidth = min(0.8, nbars/(nbars+1.5));	

  				end

			end
			

			%----- 1: Cue13(8hz) Vs. Cue34(8hz) ------------
			% size: 1-freBins; 2-chan; 3-sub; 4-condition;
			subplot(subNumPerFig,subfigNumPerRow,subfigNumPerRow*(iSub-1-(iPage-1)*subNumPerFig)+1);
			if iSub<=size(allPowerMinusBase_4con,3)
					cPower = squeeze(allFrePointPowerMinusBase_4con(1,chooseChanPower(iSub,:,1),iSub,[1 3])); 
					Power_con13_8hz(:,:,iSub) = cPower;
			else
					cPower = Power_con13_8hz(:,:,iSub);
					cError = Power_con13_8hz(:,:,iSub+1);
			end
			H=bar(cPower,'BarWidth',0.3);
			hold on;
			% axis([0.7 3.3 -0.2 0.2]);
			box off;
			ylim(plotPowerRange);

			if iSub==size(allPowerMinusBase_4con,3)+1
				subStr = 'allSubMean';
				for i = 1:nbars
				    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
				    errorbar(x,cPower(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
				    hold on;
				end
			elseif iSub<=size(allPowerMinusBase_4con,3)
				subStr = ['Sub',num2str(iSub)];
				set(gca,'XTickLabel',standard_eName(chooseChanPower(iSub,:,1)));
			end
			title([subStr,'  cue13&cue34: 8hz']);
			clear cPower cError;



			%----- 2: Cue42(20hz) Vs. Cue21(20hz) ------------
			% size: 1-freBins; 2-chan; 3-sub; 4-condition;
			subplot(subNumPerFig,subfigNumPerRow,subfigNumPerRow*(iSub-1-(iPage-1)*subNumPerFig)+2);
			if iSub<=size(allPowerMinusBase_4con,3)
					cPower = squeeze(allFrePointPowerMinusBase_4con(2,chooseChanPower(iSub,:,2),iSub,[4 2])); 
					Power_con42_20hz(:,:,iSub) = cPower;
			else
				    cPower = Power_con42_20hz(:,:,iSub); 
				    cError = Power_con42_20hz(:,:,iSub+1); 				    
			end	
			H=bar(cPower,'BarWidth',0.3);
			hold on;box off;
			ylim(plotPowerRange);
			if iSub<=size(allPowerMinusBase_4con,3)
					set(gca,'XTickLabel',standard_eName(chooseChanPower(iSub,:,2)));
			else
					for i = 1:nbars
					    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
					    errorbar(x,cPower(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
					    hold on;
					end	
			end
			title(['  cue42&cue21: 20hz']);
			clear cPower  cError;



			%----- 3: Cue13(8hz)+Cue42(20hz) Vs. Cue34(8hz)+Cue21(20hz) ------------
			% size: 1-freBins; 2-chan; 3-sub; 4-condition;
			subplot(subNumPerFig,subfigNumPerRow,subfigNumPerRow*(iSub-1-(iPage-1)*subNumPerFig)+3);
			cPower = mean(cat(3,Power_con13_8hz(:,:,iSub),Power_con42_20hz(:,:,iSub)),3);
			Power_con13_8hz_con42_20hz(:,:,iSub) = cPower;
			if iSub>size(allPowerMinusBase_4con,3)

					cError = Power_con13_8hz_con42_20hz(:,:,iSub+1);
					for i = 1:nbars
					    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
					    errorbar(x,cPower(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
					    hold on;
					end	
			end
			H=bar(cPower,'BarWidth',0.3);
			ylim(plotPowerRange);
			hold on;box off;
			title([' mean of fig1&2']);
			clear cPower  cError;



			%----- 4: Cue13(8hz)+Cue42(8hz) Vs. Cue34(8hz)+Cue21(8hz) ------------
			% size: 1-freBins; 2-chan; 3-sub; 4-condition;
			subplot(subNumPerFig,subfigNumPerRow,subfigNumPerRow*(iSub-1-(iPage-1)*subNumPerFig)+4);
			if iSub<=size(allPowerMinusBase_4con,3)
					cPower(:,1) = squeeze(mean(allFrePointPowerMinusBase_4con(1,chooseChanPower(iSub,:,1),iSub,[1 4]),4)); 
					cPower(:,2) = squeeze(mean(allFrePointPowerMinusBase_4con(1,chooseChanPower(iSub,:,1),iSub,[2 3]),4)); 
					Power_all_8hz(:,:,iSub) = cPower;
			else
					cPower = Power_all_8hz(:,:,iSub);
					cError = Power_all_8hz(:,:,iSub+1);

			end	
			H=bar(cPower,'BarWidth',0.3);
			ylim(plotPowerRange);
			hold on;box off;

			if iSub<=size(allPowerMinusBase_4con,3)
					set(gca,'XTickLabel',standard_eName(chooseChanPower(iSub,:,1)));
			else
					for i = 1:nbars
					    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
					    errorbar(x,cPower(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
					    hold on;
					end	
			end
			title(['  sepCue&ajCue: 8hz']);
			clear cPower  cError;



			%----- 5: Cue13(20hz)+Cue42(20hz) Vs. Cue34(20hz)+Cue21(20hz) ------------
			% size: 1-freBins; 2-chan; 3-sub; 4-condition;
			subplot(subNumPerFig,subfigNumPerRow,subfigNumPerRow*(iSub-1-(iPage-1)*subNumPerFig)+5);
			if iSub<=size(allPowerMinusBase_4con,3)
					cPower(:,1) = squeeze(mean(allFrePointPowerMinusBase_4con(2,chooseChanPower(iSub,:,2),iSub,[1 4]),4)); 
					cPower(:,2) = squeeze(mean(allFrePointPowerMinusBase_4con(2,chooseChanPower(iSub,:,2),iSub,[2 3]),4)); 
					Power_all_20hz(:,:,iSub) = cPower;
			else
					cPower = Power_all_20hz(:,:,iSub);
					cError = Power_all_20hz(:,:,iSub+1);
			end
			H=bar(cPower,'BarWidth',0.3);
			ylim(plotPowerRange);
			hold on;box off;
			if iSub<=size(allPowerMinusBase_4con,3)
					set(gca,'XTickLabel',standard_eName(chooseChanPower(iSub,:,2)));
			else
					for i = 1:nbars
					    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
					    errorbar(x,cPower(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
					    hold on;
					end				
			end
			title(['  sepCue&ajCue: 20hz']);
			clear cPower  cError;



			%----- 6: Cue13(20hz)+Cue42(8hz) Vs. Cue34(8hz)+Cue21(8hz) ------------
			% size: 1-freBins; 2-chan; 3-sub; 4-condition;
			subplot(subNumPerFig,subfigNumPerRow,subfigNumPerRow*(iSub-1-(iPage-1)*subNumPerFig)+6);
			cPower = mean(cat(3,Power_all_8hz(:,:,iSub),Power_all_20hz(:,:,iSub)),3);
			Power_all(:,:,iSub) = cPower;
			
			if iSub>size(allPowerMinusBase_4con,3)

					cError = Power_all(:,:,iSub+1);
					for i = 1:nbars
					    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
					    errorbar(x,cPower(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
					    hold on;
					end	
			end

			H=bar(cPower,'BarWidth',0.3);
			ylim(plotPowerRange);
			hold on;box off;
			title(['   mean of fig4&5']);
			clear cPower  cError;



			%----- 7: mergeCondition(to imporove SNR),then do fft...SepCue(8hz) Vs. ajCue(8hz) ------------
			% size: 1-freBins; 2-chan; 3-sub; 4-condition;
			subplot(subNumPerFig,subfigNumPerRow,subfigNumPerRow*(iSub-1-(iPage-1)*subNumPerFig)+7);
			if iSub<=size(allPowerMinusBase_2con,3)
					cPower = squeeze(allFrePointPowerMinusBase_2con(1,chooseChanPower(iSub,:,1),iSub,:)); 
					Power_8hz_2con(:,:,iSub) = cPower;
			else
					cPower = Power_8hz_2con(:,:,iSub);
					cError = Power_8hz_2con(:,:,iSub+1);
			end	
			H=bar(cPower,'BarWidth',0.3);
			ylim(plotPowerRange);
			hold on;box off;

			if iSub<=size(allPowerMinusBase_2con,3)
					set(gca,'XTickLabel',standard_eName(chooseChanPower(iSub,:,1)));
			else
					for i = 1:nbars
					    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
					    errorbar(x,cPower(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
					    hold on;
					end				
			end
			title([' sepCue&ajCue fft: 8hz']);
			clear cPower  cError;



			%----- 8: mergeCondition(to imporove SNR),then do fft...SepCue(20hz) Vs. ajCue(20hz) ------------
			% size: 1-freBins; 2-chan; 3-sub; 4-condition;
			subplot(subNumPerFig,subfigNumPerRow,subfigNumPerRow*(iSub-1-(iPage-1)*subNumPerFig)+8);
			if iSub<=size(allPowerMinusBase_2con,3)
					cPower = squeeze(allFrePointPowerMinusBase_2con(2,chooseChanPower(iSub,:,2),iSub,:)); 
					Power_20hz_2con(:,:,iSub) = cPower;
			else
					cPower = Power_20hz_2con(:,:,iSub);
					cError = Power_20hz_2con(:,:,iSub+1);
			end	

			H=bar(cPower,'BarWidth',0.3);
			ylim(plotPowerRange);
			hold on;box off;
			if iSub<=size(allPowerMinusBase_2con,3)
					set(gca,'XTickLabel',standard_eName(chooseChanPower(iSub,:,2)));
			else
					for i = 1:nbars
					    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
					    errorbar(x,cPower(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
					    hold on;
					end				
			end
			title(['  sepCue&ajCue fft: 20hz']);
			clear cPower  cError;



			%----- 9: mergeCondition(to imporove SNR),then do fft...SepCue(8+20hz) Vs. ajCue(8+20hz) ------------
			% size: 1-freBins; 2-chan; 3-sub; 4-condition;
			subplot(subNumPerFig,subfigNumPerRow,subfigNumPerRow*(iSub-1-(iPage-1)*subNumPerFig)+9);
			cPower = mean(cat(3,Power_8hz_2con(:,:,iSub),Power_20hz_2con(:,:,iSub)),3);
			Power_all_2con(:,:,iSub) = cPower;
			
			if iSub>size(allPowerMinusBase_4con,3)

					cError = Power_all_2con(:,:,iSub+1);
					for i = 1:nbars
					    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
					    errorbar(x,cPower(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
					    hold on;
					end	
			end			
			H=bar(cPower,'BarWidth',0.3);
			ylim(plotPowerRange);
			hold on;box off;
			title(['   mean of fig7&8']);
			clear cPower  cError;
		end
	end
	legend('sepCue','ajCue','Location','South');
	legend('boxoff');
end


 

% for iSub = 1:size(SNR_con13_8hz,3)-2
% 	cSubData            = squeeze(SNR_con13_8hz(:,:,iSub))';
% 	spssfor8hz_SNR(iSub,:)  =  cSubData(:)';
% 	clear cSubData
% 	cSubData            = squeeze(SNR_con42_20hz(:,:,iSub))';
% 	spssfor20hz_SNR(iSub,:) =  cSubData(:)';
% 	clear cSubData
% 	cSubData            = squeeze(SNR_con13_8hz_con42_20hz(:,:,iSub))';
% 	spssforBoth_SNR(iSub,:) =  cSubData(:)';
% 	clear cSubData
% end





