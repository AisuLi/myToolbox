

%%%%%%%%%%%%%%%%%%%%%%%%%		20180326  find the max chan		%%%%%%%%%%%%%%%%%%%%%%%%%%%		



load('D:\expData\SSVEP_8loc\Exp1_20170328\Data\scanData\20180201\cnt\20180325\mergeAllcon_-500To1200.mat')
clearvars -except allPower allSNR freBins

directory = 'D:\expData\SSVEP_8loc\Exp1_20170328\Data\scanData\20180201\cnt';
allFile   = dir(fullfile(directory,'*BA3.eeg'));
fftType   = 'avgFFT'; % 'trialFFT'
cd(directory);
freRange  = [5,25];
targetFre = [12;15]; % for loc2, loc3
[P1eNo,standard_eName] = getElectrodeNo_bcl(1,[],[]);


%-------find closest frePoint-------/
for iFre = 1:size(targetFre,1)
	[X, freIndex]         = sort(abs(freBins-targetFre(iFre,1)),'ascend');
	freTimePoint(iFre,1)  = freIndex(1);
	realFre(iFre,1) 	  = freBins(freIndex(1));
end
%-----------------------------------\
frePointSNR   = allSNR(freTimePoint(:),:,:);
frePointPower = allPower(freTimePoint(:),:,:);


save('mergeAllcon_-500To1200');

%----------------  have found the chan, it's the end  -------------\




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




% clear all;
% load('mergeAllcon_-500To1200');
% % load('fft_2con_base&postBase.mat');
% load('fft_4con_base&postBase');
% [P1eNo,standard_eName] = getElectrodeNo_bcl(1,[],[]);
% subNumPerFig    = 5;
% subfigNumPerRow = 6;
% figNum          = ceil(size(allPowerMinusBase_4con,3)./subNumPerFig);


% barColor = repmat([0.9,1,0.1,0.3]',1,3);
% %%%%%%   bar of SNR %%%%%%%%

% plotSNRRange = [-2 3; -0.8 0.8];
% for iPage = 1:figNum
% 	figure;
% 	set(gcf,'Position',get(0,'ScreenSize'),'color','w');
% 	set(gcf,'color','white')
	
% 	cSub = [1:subNumPerFig]+(iPage-1)*subNumPerFig;
% 	% suptitle('SNR');
% 	for iSub = cSub
% 		if iSub<=size(allSNRMinusBase_4con,3)+1

% 			if iSub==size(allSNRMinusBase_4con,3)+1 % compute mean 

% 				for iChan = 1:size(SNR_con13_8hz,1)
%   						[gm,gsem] = grpstats(squeeze(SNR_con13_8hz(iChan,:,[1:size(allSNRMinusBase_4con,3)]))',{},{'mean','sem'});
%   						SNR_con13_8hz(iChan,:,size(allSNRMinusBase_4con,3)+1) = gm;
%   						SNR_con13_8hz(iChan,:,size(allSNRMinusBase_4con,3)+2) = gsem;

  						
%   						[gm,gsem] = grpstats(squeeze(SNR_con42_20hz(iChan,:,[1:size(allSNRMinusBase_4con,3)]))',{},{'mean','sem'});
%   						SNR_con42_20hz(iChan,:,size(allSNRMinusBase_4con,3)+1) = gm;
%   						SNR_con42_20hz(iChan,:,size(allSNRMinusBase_4con,3)+2) = gsem;
  						
%   						[gm,gsem] = grpstats(squeeze(SNR_con13_8hz_con42_20hz(iChan,:,[1:size(allSNRMinusBase_4con,3)]))',{},{'mean','sem'});
%   						SNR_con13_8hz_con42_20hz(iChan,:,size(allSNRMinusBase_4con,3)+1) = gm;
%   						SNR_con13_8hz_con42_20hz(iChan,:,size(allSNRMinusBase_4con,3)+2) = gsem;

% 						ngroups    = 3;
% 						nbars      = 4;
% 						groupwidth = min(0.8, nbars/(nbars+1.5));	

%   				end

% 			end
			

% 			%----- 1: Cue13(8hz) Vs. Cue34(8hz) ------------
% 			% size: 1-freBins; 2-chan; 3-sub; 4-condition;
% 			subplot(subNumPerFig,subfigNumPerRow,subfigNumPerRow*(iSub-1-(iPage-1)*subNumPerFig)+1);
% 			if iSub<=size(allSNRMinusBase_4con,3)
% 					cSNR = squeeze(allFrePointSNRMinusBase_4con(1,chooseChanSNR(iSub,:,1),iSub,[1 3 2 4])); 
% 					SNR_con13_8hz(:,:,iSub) = cSNR;
% 			else
% 					cSNR = SNR_con13_8hz(:,:,iSub);
% 					cError = SNR_con13_8hz(:,:,iSub+1);
% 			end
% 			H=bar(cSNR,'grouped');
% 		    ch = get(H,'children');
%             set(ch{1},'Facecolor',barColor(1,:));
%             set(ch{2},'Facecolor',barColor(2,:));
%             set(ch{3},'Facecolor',barColor(3,:));
%             set(ch{4},'Facecolor',barColor(4,:));

% 			hold on;
% 			% axis([0.7 3.3 -0.2 0.2]);
% 			box off;
% 			ylim(plotSNRRange(1,:));

% 			if iSub==size(allSNRMinusBase_4con,3)+1
% 				subStr = 'allSubMean';
% 				for i = 1:nbars
% 				    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
% 				    errorbar(x,cSNR(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
% 				    hold on;
% 					ylim(plotSNRRange(2,:));
% 				end
% 			elseif iSub<=size(allSNRMinusBase_4con,3)
% 				subStr = ['Sub',num2str(iSub)];
% 				set(gca,'XTickLabel',standard_eName(chooseChanSNR(iSub,:,1)));
% 			end
% 			title([subStr,'  8hz']);
% 			clear cSNR cError;



% 			%----- 2: Cue42(20hz) Vs. Cue21(20hz) ------------
% 			% size: 1-freBins; 2-chan; 3-sub; 4-condition;
% 			subplot(subNumPerFig,subfigNumPerRow,subfigNumPerRow*(iSub-1-(iPage-1)*subNumPerFig)+2);
% 			if iSub<=size(allSNRMinusBase_4con,3)
% 					cSNR = squeeze(allFrePointSNRMinusBase_4con(2,chooseChanSNR(iSub,:,2),iSub,[4 2 3 1])); 
% 					SNR_con42_20hz(:,:,iSub) = cSNR;
% 			else
% 				    cSNR = SNR_con42_20hz(:,:,iSub); 
% 				    cError = SNR_con42_20hz(:,:,iSub+1); 				    
% 			end	
% 			H=bar(cSNR,'grouped');
% 		    ch = get(H,'children');
%             set(ch{1},'Facecolor',barColor(1,:));
%             set(ch{2},'Facecolor',barColor(2,:));
%             set(ch{3},'Facecolor',barColor(3,:));
%             set(ch{4},'Facecolor',barColor(4,:));

% 			hold on;box off;
% 			ylim(plotSNRRange(1,:));
% 			if iSub<=size(allSNRMinusBase_4con,3)
% 					set(gca,'XTickLabel',standard_eName(chooseChanSNR(iSub,:,2)));
% 			else
% 					for i = 1:nbars
% 					    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
% 					    errorbar(x,cSNR(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
% 					    hold on;
% 					  	ylim(plotSNRRange(2,:));
% 					end	
% 			end
% 			title(['   20hz']);
% 			clear cSNR  cError;



% 			%----- 3: Cue13(8hz)+Cue42(20hz) Vs. Cue34(8hz)+Cue21(20hz) ------------
% 			% size: 1-freBins; 2-chan; 3-sub; 4-condition;
% 			subplot(subNumPerFig,subfigNumPerRow,subfigNumPerRow*(iSub-1-(iPage-1)*subNumPerFig)+3);
% 			cSNR = mean(cat(3,SNR_con13_8hz(:,:,iSub),SNR_con42_20hz(:,:,iSub)),3);
% 			SNR_con13_8hz_con42_20hz(:,:,iSub) = cSNR;
% 			H=bar(cSNR,'grouped');
% 		    ch = get(H,'children');
%             set(ch{1},'Facecolor',barColor(1,:));
%             set(ch{2},'Facecolor',barColor(2,:));
%             set(ch{3},'Facecolor',barColor(3,:));
%             set(ch{4},'Facecolor',barColor(4,:));

% 			ylim(plotSNRRange(1,:));
% 			hold on;box off;
% 			if iSub>size(allSNRMinusBase_4con,3)

% 					cError = SNR_con13_8hz_con42_20hz(:,:,iSub+1);
% 					for i = 1:nbars
% 					    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
% 					    errorbar(x,cSNR(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
% 					    hold on;
% 					end	
% 					ylim(plotSNRRange(2,:));
% 			end
% 			title([' mean of 8hz and 20hz']);
% 			clear cSNR  cError;
% 		end
% 	end
% 	l = legend([ch{1} ch{2} ch{3} ch{4} ],'sepCue--Uncue','ajCue--Uncue','ajCue--Cue','SepCue--Cue');
%     legend('boxoff');
% 	set(l,'Position',[0.511631944444441 0.664700098328417 0.0765625 0.0780072107505736]);
% end





% %%%%%%   bar of Power %%%%%%%%

% plotPowerRange = [-0.3 0.1; -0.15 0];
% for iPage = 1:figNum
% 	figure;
% 	set(gcf,'Position',get(0,'ScreenSize'),'color','w');
% 	set(gcf,'color','white')
	
% 	cSub = [1:subNumPerFig]+(iPage-1)*subNumPerFig;
% 	% suptitle('Power');
% 	for iSub = cSub
% 		if iSub<=size(allPowerMinusBase_4con,3)+1

% 			if iSub==size(allPowerMinusBase_4con,3)+1 % compute mean 

% 				for iChan = 1:size(Power_con13_8hz,1)
%   						[gm,gsem] = grpstats(squeeze(Power_con13_8hz(iChan,:,[1:size(allPowerMinusBase_4con,3)]))',{},{'mean','sem'});
%   						Power_con13_8hz(iChan,:,size(allPowerMinusBase_4con,3)+1) = gm;
%   						Power_con13_8hz(iChan,:,size(allPowerMinusBase_4con,3)+2) = gsem;

  						
%   						[gm,gsem] = grpstats(squeeze(Power_con42_20hz(iChan,:,[1:size(allPowerMinusBase_4con,3)]))',{},{'mean','sem'});
%   						Power_con42_20hz(iChan,:,size(allPowerMinusBase_4con,3)+1) = gm;
%   						Power_con42_20hz(iChan,:,size(allPowerMinusBase_4con,3)+2) = gsem;
  						
%   						[gm,gsem] = grpstats(squeeze(Power_con13_8hz_con42_20hz(iChan,:,[1:size(allPowerMinusBase_4con,3)]))',{},{'mean','sem'});
%   						Power_con13_8hz_con42_20hz(iChan,:,size(allPowerMinusBase_4con,3)+1) = gm;
%   						Power_con13_8hz_con42_20hz(iChan,:,size(allPowerMinusBase_4con,3)+2) = gsem;

% 						ngroups    = 3;
% 						nbars      = 4;
% 						groupwidth = min(0.8, nbars/(nbars+1.5));	

%   				end

% 			end
			

% 			%----- 1: Cue13(8hz) Vs. Cue34(8hz) ------------
% 			% size: 1-freBins; 2-chan; 3-sub; 4-condition;
% 			subplot(subNumPerFig,subfigNumPerRow,subfigNumPerRow*(iSub-1-(iPage-1)*subNumPerFig)+1);
% 			if iSub<=size(allPowerMinusBase_4con,3)
% 					cPower = squeeze(allFrePointPowerMinusBase_4con(1,chooseChanPower(iSub,:,1),iSub,[1 3 2 4])); 
% 					Power_con13_8hz(:,:,iSub) = cPower;
% 			else
% 					cPower = Power_con13_8hz(:,:,iSub);
% 					cError = Power_con13_8hz(:,:,iSub+1);
% 			end
% 			H=bar(cPower,'grouped');
% 		    ch = get(H,'children');
%             set(ch{1},'Facecolor',barColor(1,:));
%             set(ch{2},'Facecolor',barColor(2,:));
%             set(ch{3},'Facecolor',barColor(3,:));
%             set(ch{4},'Facecolor',barColor(4,:));

% 			hold on;
% 			% axis([0.7 3.3 -0.2 0.2]);
% 			box off;
% 			ylim(plotPowerRange(1,:));

% 			if iSub==size(allPowerMinusBase_4con,3)+1
% 				subStr = 'allSubMean';
% 				for i = 1:nbars
% 				    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
% 				    errorbar(x,cPower(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
% 				    hold on;
% 					ylim(plotPowerRange(2,:));
% 				end
% 			elseif iSub<=size(allPowerMinusBase_4con,3)
% 				subStr = ['Sub',num2str(iSub)];
% 				set(gca,'XTickLabel',standard_eName(chooseChanPower(iSub,:,1)));
% 			end
% 			title([subStr,'  8hz']);
% 			clear cPower cError;



% 			%----- 2: Cue42(20hz) Vs. Cue21(20hz) ------------
% 			% size: 1-freBins; 2-chan; 3-sub; 4-condition;
% 			subplot(subNumPerFig,subfigNumPerRow,subfigNumPerRow*(iSub-1-(iPage-1)*subNumPerFig)+2);
% 			if iSub<=size(allPowerMinusBase_4con,3)
% 					cPower = squeeze(allFrePointPowerMinusBase_4con(2,chooseChanPower(iSub,:,2),iSub,[4 2 3 1])); 
% 					Power_con42_20hz(:,:,iSub) = cPower;
% 			else
% 				    cPower = Power_con42_20hz(:,:,iSub); 
% 				    cError = Power_con42_20hz(:,:,iSub+1); 				    
% 			end	
% 			H=bar(cPower,'grouped');
% 		    ch = get(H,'children');
%             set(ch{1},'Facecolor',barColor(1,:));
%             set(ch{2},'Facecolor',barColor(2,:));
%             set(ch{3},'Facecolor',barColor(3,:));
%             set(ch{4},'Facecolor',barColor(4,:));

% 			hold on;box off;
% 			ylim(plotPowerRange(1,:));
% 			if iSub<=size(allPowerMinusBase_4con,3)
% 					set(gca,'XTickLabel',standard_eName(chooseChanPower(iSub,:,2)));
% 			else
% 					for i = 1:nbars
% 					    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
% 					    errorbar(x,cPower(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
% 					    hold on;
% 					  	ylim(plotPowerRange(2,:));
% 					end	
% 			end
% 			title(['   20hz']);
% 			clear cPower  cError;



% 			%----- 3: Cue13(8hz)+Cue42(20hz) Vs. Cue34(8hz)+Cue21(20hz) ------------
% 			% size: 1-freBins; 2-chan; 3-sub; 4-condition;
% 			subplot(subNumPerFig,subfigNumPerRow,subfigNumPerRow*(iSub-1-(iPage-1)*subNumPerFig)+3);
% 			cPower = mean(cat(3,Power_con13_8hz(:,:,iSub),Power_con42_20hz(:,:,iSub)),3);
% 			Power_con13_8hz_con42_20hz(:,:,iSub) = cPower;
% 			H=bar(cPower,'grouped');
% 		    ch = get(H,'children');
%             set(ch{1},'Facecolor',barColor(1,:));
%             set(ch{2},'Facecolor',barColor(2,:));
%             set(ch{3},'Facecolor',barColor(3,:));
%             set(ch{4},'Facecolor',barColor(4,:));

% 			ylim(plotPowerRange(1,:));
% 			hold on;box off;
% 			if iSub>size(allPowerMinusBase_4con,3)

% 					cError = Power_con13_8hz_con42_20hz(:,:,iSub+1);
% 					for i = 1:nbars
% 					    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
% 					    errorbar(x,cPower(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
% 					    hold on;
% 					end	
% 					ylim(plotPowerRange(2,:));
% 			end
% 			title([' mean of 8hz and 20hz']);
% 			clear cPower  cError;
% 		end
% 	end
% 	l = legend([ch{1} ch{2} ch{3} ch{4} ],'sepCue--Uncue','ajCue--Uncue','ajCue--Cue','SepCue--Cue');
%     legend('boxoff');
% 	set(l,'Position',[0.511631944444441 0.664700098328417 0.0765625 0.0780072107505736]);
% end








%%%%%%%%%%%%%%%%%%%%%%%%%	20180328   avg of allSub, then find the max 3 chan		%%%%%%%%%%%%%%%%%%%%%%%%%%%		

%%%--------- do fft of [-500 1200] for 4 conditions,separately -------/
directory              = 'D:\expData\SSVEP_8loc\Exp1_20170328\Data\scanData\20180201\cnt';
allFile                = dir(fullfile(directory,'*BA3.eeg'));
fftType                = 'avgFFT'; % 'trialFFT'
cd(directory);
timeBin                = [-500 1200];

stimTrigger            = [1:10;11:20;21:30;31:40]; % 8Hz: cue13;21;34;42

targetFre              = [8;20]; % for loc2, loc3
freRange               = [5,25];
[P1eNo,standard_eName] = getElectrodeNo_bcl(1,[],[]);

for iSub = 1:size(allFile,1)

	eegName = allFile(iSub).name;

	for iCondition = 1:size(stimTrigger,1)
	
		for iTimeRange = 1:size(timeBin,1)  % 
			
			[iSubAvg,freBins,Power,SNR] = eegfftSNR_bcl(eegName,freRange,timeBin(iTimeRange,:),stimTrigger(iCondition,:),fftType); 
			allPower_wholeEpoch_4con(:,:,iSub,iCondition,iTimeRange)  = Power; % size: 1-freBins; 2-chan; 3-sub; 4-condition; 5-timeRange(base)
			allSNR_wholeEpoch_4con(:,:,iSub,iCondition,iTimeRange)    = SNR;
			clear  Power SNR iSubAvg;
		end
	end
end 

%-------  Do fft for 4 conditions,separately. It's the end---------\



%-------find closest frePoint-------/
for iFre = 1:size(targetFre,1)
	[X, freIndex]         = sort(abs(freBins-targetFre(iFre,1)),'ascend');
	freTimePoint(iFre,1)  = freIndex(1);
	realFre(iFre,1) 	  = freBins(freIndex(1));
end
%-------have found, the end---------\


%%%%%%%%%%%%%%%%%%%%%%%  adjust the conditions	%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allFrePointSNR_wholeEpoch_4con   = cat(1,allSNR_wholeEpoch_4con(freTimePoint(1),:,:,[1 3 2 4]),allSNR_wholeEpoch_4con(freTimePoint(2),:,:,[4 2 3 1])); % size: 1-freBins; 2-chan; 3-sub; 4-condition;
allFrePointPower_wholeEpoch_4con = cat(1,allPower_wholeEpoch_4con(freTimePoint(1),:,:,[1 3 2 4]),allPower_wholeEpoch_4con(freTimePoint(2),:,:,[4 2 3 1])); % size: 1-freBins; 2-chan; 3-sub; 4-condition;

allFrePointSNR_wholeEpoch_4con(:,:,size(allSNR_wholeEpoch_4con,3)+1,:)     = mean(allFrePointSNR_wholeEpoch_4con,3);
allFrePointPower_wholeEpoch_4con(:,:,size(allPower_wholeEpoch_4con,3)+1,:) = mean(allFrePointPower_wholeEpoch_4con,3);
%%%%%%%%%%%%%%%%%%%%%%%  have adjusted the condition, the end %%%%%%%%%%%%%

save('fft_4con_-500To1200');


%%%%---------plot topo of SNR & Power of 4conditions,left:8hz right:20hz-----------/
waveNum = 3;
load('fft_4con_-500To1200');
[myColormap] = makeColormap(20,0,1,'r');
subNumPerFig = 5;
figNum       = ceil(size(allFile,1)./subNumPerFig);
[P1eNo,standard_eName]=getElectrodeNo_bcl(0,'haveM1',[]);


% %%------------plot topo of SNR------------/
% for iPage = 1:figNum
% 	figure;
% 	set(gcf,'Position',get(0,'ScreenSize'),'color','w');
% 	suptitle(['			left:  '  num2str(realFre(1))  'Hz       cue13 cue34 cue21 Cue42                        ',   '       SNR     -500To1200                '     	   '   				       right:    ',num2str(realFre(2))  'Hz       cue42 cue34 cue21 Cue13   ']);
% 	cSub = [1:subNumPerFig]+(iPage-1)*subNumPerFig; 

% 	for iSub = cSub
		
% 		if iSub<=size(allFrePointSNR_wholeEpoch_4con,3)
% 			for istim = 1:size(targetFre,1)

% 				for iCondition = 1:size(allFrePointSNR_wholeEpoch_4con,4)

% 					if istim==1
% 						figureLoc = 9*(iSub-1-(iPage-1)*subNumPerFig)+iCondition;
% 						subplot(subNumPerFig,9,9*(iSub-1-(iPage-1)*subNumPerFig)+5);

% 						if iSub<size(allFrePointSNR_wholeEpoch_4con,3)
% 							title(['sub:' num2str(iSub)]);
% 						else
% 							title('allSubMean');
% 						end
		
% 					elseif istim==2
% 						figureLoc = 9*(iSub-1-(iPage-1)*subNumPerFig)+5+iCondition;
% 					end
% 					subplot(subNumPerFig,9,figureLoc);
% 					topoplot_bcl(allFrePointSNR_wholeEpoch_4con(istim,:,iSub,iCondition),'EEGChans',standard_eName,'maplimitsDouble',[0 3],'colormap',myColormap);
					
% 				end
% 			end
% 		end
% 	end
% 	set(gcf,'color','white');
% end
% %%----plot topo of SNR,it's the end------\


% %%--------------plot topo of Power----------------/
% [myColormap] = makeColormap(20,0,1,'b');
% for iPage = 1:figNum
% 	figure;
% 	set(gcf,'Position',get(0,'ScreenSize'),'color','w');
% 	set(gcf,'color','white');
% 	suptitle(['			left:  '  num2str(realFre(1))  'Hz       cue13 cue21 cue34 Cue42                        ',   '       Power     -500To1200                '     	   '   				       right:    ',num2str(realFre(2))  'Hz       cue13 cue21 cue34 Cue42   ']);
% 	cSub = [1:subNumPerFig]+(iPage-1)*subNumPerFig; 

% 	for iSub = cSub
		
% 		if iSub<=size(allFrePointPower_wholeEpoch_4con,3)
% 			for istim = 1:size(targetFre,1)

% 				for iCondition = 1:size(allFrePointPower_wholeEpoch_4con,4)

% 					if istim==1
% 						figureLoc = 9*(iSub-1-(iPage-1)*subNumPerFig)+iCondition;
% 						subplot(subNumPerFig,9,9*(iSub-1-(iPage-1)*subNumPerFig)+5);

% 						if iSub<size(allFrePointPower_wholeEpoch_4con,3)
% 							title(['sub:' num2str(iSub)]);
% 						else
% 							title('allSubMean');
% 						end
		
% 					elseif istim==2
% 						figureLoc = 9*(iSub-1-(iPage-1)*subNumPerFig)+5+iCondition;
% 					end
% 					subplot(subNumPerFig,9,figureLoc);
% 					topoplot_bcl(allFrePointPower_wholeEpoch_4con(istim,:,iSub,iCondition),'EEGChans',standard_eName,'maplimitsDouble',[0 0.2],'colormap',myColormap);
					
% 				end
% 			end
% 		end
% 	end
% 	set(gcf,'color','white');
% end
% %%-----plot topo of Power,it's the end----------\





%%%%%----- find the max 3 chans of [mean of allSub's fft]--------/ 
load('fft_4con_-500To1200');
load('D:\expData\SSVEP_8loc\Exp1_20170328\Data\scanData\20180201\cnt\20180325\fft_4con_base&postBase.mat');

for iFre = 1:size(targetFre,1)
	for iCondition = 1:size(allFrePointSNR_wholeEpoch_4con,4)
		[X, chanIndex]                          = sort(allFrePointSNR_wholeEpoch_4con(iFre,[1:size(standard_eName,2)-2],end,iCondition),'descend');						
		chooseChanSNR_4con(:,iFre,iCondition)   = chanIndex(1:3);	%%% size: 1-e, 2-fre, 3-conditions	
		clear X chanIndex
				
		[X, chanIndex]                          = sort(allFrePointPower_wholeEpoch_4con(iFre,[1:size(standard_eName,2)-2],end,iCondition),'descend');						
		chooseChanPower_4con(:,iFre,iCondition) = chanIndex(1:3); %%%% size: 1-e, 2-fre, 3-conditions	
		clear X  chanIndex
	end
end
%---------- have found max 3 chans, it's the end-------------\
		

%%%%%%%%  adjust the condition %%%%%%%%%%
%--------------fft of 4 conditions for [-500 -1] & [201 1200]----------------/
allSNRDivideBase_4con           = squeeze(allSNR_4con(:,:,:,:,2)./allSNR_4con(:,:,:,:,1)); % size: 1-freBins; 2-chan; 3-sub; 4-condition;
allFrePointSNRDivideBase_4con   = cat(1,allSNRDivideBase_4con(freTimePoint(1),:,:,[1 3 2 4]),allSNRDivideBase_4con(freTimePoint(2),:,:,[4 2 3 1]));
allPowerDivideBase_4con         = squeeze(allPower_4con(:,:,:,:,2)./allPower_4con(:,:,:,:,1)); % size: 1-freBins; 2-chan; 3-sub; 4-condition;
allFrePointPowerDivideBase_4con = cat(1,allPowerDivideBase_4con(freTimePoint(1),:,:,[1 3 2 4]),allPowerDivideBase_4con(freTimePoint(2),:,:,[4 2 3 1]));
%%%%%%%%  have adjusted the condition %%%%%%%%%%
% size: 1-freBins; 2-chan; 3-sub; 4-condition



%%%% 
for iFre = 1:size(chooseChanSNR_4con,2)
	for iCondition = 1:size(chooseChanSNR_4con,3)
		SNR_4con_8hz20hz(iFre,:,:,iCondition)   = allFrePointSNRDivideBase_4con(iFre,chooseChanSNR_4con(:,iFre,iCondition),:,iCondition);
		Power_4con_8hz20hz(iFre,:,:,iCondition) = allFrePointPowerDivideBase_4con(iFre,chooseChanPower_4con(:,iFre,iCondition),:,iCondition);
	end
end


SNR_4con_8hz20hz(3,:,:,:)                                   = mean(SNR_4con_8hz20hz,1); 
SNR_4con_8hz20hz(:,size(SNR_4con_8hz20hz,2)+1,:,:)          = mean(SNR_4con_8hz20hz,2);
SNR_4con_8hz20hz(:,:,size(allSNRDivideBase_4con,3)+1,:)     = mean(SNR_4con_8hz20hz,3); 

Power_4con_8hz20hz(3,:,:,:)                                 = mean(Power_4con_8hz20hz,1); 
Power_4con_8hz20hz(:,size(Power_4con_8hz20hz,2)+1,:,:)      = mean(Power_4con_8hz20hz,2);
Power_4con_8hz20hz(:,:,size(allPowerDivideBase_4con,3)+1,:) = mean(Power_4con_8hz20hz,3); 

for iFre = 1:size(SNR_4con_8hz20hz,1)
	for iCondition = 1:size(chooseChanSNR_4con,3)

		[gm,gsem] = grpstats(squeeze(SNR_4con_8hz20hz(iFre,:,[1:size(allSNRDivideBase_4con,3)],iCondition))',{},{'mean','sem'});
		SNR_4con_8hz20hz_sem(iFre,:,iCondition) = gsem;
		clear gm gsem
		[gm,gsem] = grpstats(squeeze(Power_4con_8hz20hz(iFre,:,[1:size(allPowerDivideBase_4con,3)],iCondition))',{},{'mean','sem'});
		Power_4con_8hz20hz_sem(iFre,:,iCondition) = gsem;
		clear gm gsem
	end
end


%%%%%%   bar of SNR %%%%%%%%
barColor               = repmat([0.85,1,0.1,0.3]',1,3);

[P1eNo,standard_eName] = getElectrodeNo_bcl(1,[],[]);
subNumPerFig           = 5;
subfigNumPerRow        = 6;
figNum                 = ceil(size(allSNRDivideBase_4con,3)./subNumPerFig);
plotSNRRange           = [0 6; 0 4];

ngroups                = size(SNR_4con_8hz20hz,2);
nbars                  = size(SNR_4con_8hz20hz,4);
groupwidth             = min(0.8, nbars/(nbars+1.5));	


for iPage = 1:figNum
	figure;
	set(gcf,'Position',get(0,'ScreenSize'),'color','w');
	set(gcf,'color','white')
	
	cSub = [1:subNumPerFig]+(iPage-1)*subNumPerFig;
	suptitle('SNR');
	for iSub = cSub
		if iSub<=size(allSNRDivideBase_4con,3)+1

			for iFre = 1:size(SNR_4con_8hz20hz,1)
					%----- 1: Cue13(8hz) Vs. Cue34(8hz) ------------
					% size: 1-freBins; 2-chan; 3-sub; 4-condition;
					subplot(subNumPerFig,subfigNumPerRow,subfigNumPerRow*(iSub-1-(iPage-1)*subNumPerFig)+iFre);
					cSNR = squeeze(SNR_4con_8hz20hz(iFre,:,iSub,:));

					H  = bar(cSNR,'grouped');
					ch = get(H,'children');
		            set(ch{1},'Facecolor',barColor(1,:));
		            set(ch{2},'Facecolor',barColor(2,:));
		            set(ch{3},'Facecolor',barColor(3,:));
		            set(ch{4},'Facecolor',barColor(4,:));
					hold on;box off;

					if iSub==size(allSNRDivideBase_4con,3)+1  
						subStr = 'allSubMean';
						ylim(plotSNRRange(2,:));
						cError = squeeze(SNR_4con_8hz20hz_sem(iFre,:,:));   %%% size: 1-freBins; 2-chan; 3-condition;
						for i = 1:nbars
						    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
						    errorbar(x,cSNR(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
						    hold on;
							ylim(plotSNRRange(2,:));
						end

					elseif iSub<=size(allSNRDivideBase_4con,3)
						ylim(plotSNRRange(1,:));
						subStr = ['Sub',num2str(iSub)];
					end

					switch iFre
						case 1 
							freStr = '  8hz';
						case 2 
							freStr = '  20hz';						
						case 3 
							freStr = '  mean of 8hz & 20hz';		
					end	

					title([subStr,freStr]);
					clear cSNR cError;
			end
		end
	end 

	l = legend([ch{1} ch{2} ch{3} ch{4} ],'sepCue--Uncue','ajCue--Uncue','ajCue--Cue','SepCue--Cue');
    legend('boxoff');
	set(l,'Position',[0.511631944444441 0.664700098328417 0.0765625 0.0780072107505736]);
end



%%%  plot topo of Power  %%%%%%%%
plotPowerRange           = [0 4.5; 0 2.5];
for iPage = 1:figNum
	figure;
	set(gcf,'Position',get(0,'ScreenSize'),'color','w');
	set(gcf,'color','white')
	
	cSub = [1:subNumPerFig]+(iPage-1)*subNumPerFig;
	suptitle('Power');
	for iSub = cSub
		if iSub<=size(allPowerDivideBase_4con,3)+1

			for iFre = 1:size(Power_4con_8hz20hz,1)
					%----- 1: Cue13(8hz) Vs. Cue34(8hz) ------------
					% size: 1-freBins; 2-chan; 3-sub; 4-condition;
					subplot(subNumPerFig,subfigNumPerRow,subfigNumPerRow*(iSub-1-(iPage-1)*subNumPerFig)+iFre);
					cPower = squeeze(Power_4con_8hz20hz(iFre,:,iSub,:));

					H  = bar(cPower,'grouped');
					ch = get(H,'children');
		            set(ch{1},'Facecolor',barColor(1,:));
		            set(ch{2},'Facecolor',barColor(2,:));
		            set(ch{3},'Facecolor',barColor(3,:));
		            set(ch{4},'Facecolor',barColor(4,:));
					hold on;box off;

					if iSub==size(allPowerDivideBase_4con,3)+1  
						subStr = 'allSubMean';
						ylim(plotPowerRange(2,:));
						cError = squeeze(Power_4con_8hz20hz_sem(iFre,:,:));   %%% size: 1-freBins; 2-chan; 3-condition;
						for i = 1:nbars
						    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
						    errorbar(x,cPower(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
						    hold on;
							ylim(plotPowerRange(2,:));
						end

					elseif iSub<=size(allPowerDivideBase_4con,3)
						ylim(plotPowerRange(1,:));
						subStr = ['Sub',num2str(iSub)];
					end

					switch iFre
						case 1 
							freStr = '  8hz';
						case 2 
							freStr = '  20hz';						
						case 3 
							freStr = '  mean of 8hz & 20hz';		
					end	

					title([subStr,freStr]);
					clear cPower cError;
			end
		end
	end 

	l = legend([ch{1} ch{2} ch{3} ch{4} ],'sepCue--Uncue','ajCue--Uncue','ajCue--Cue','SepCue--Cue');
    legend('boxoff');
	set(l,'Position',[0.511631944444441 0.664700098328417 0.0765625 0.0780072107505736]);
end









%%%%    prepare data for spss %%%%%%
% size: 1-freBins; 2-chan; 3-sub; 4-condition;

PowerForSpss_4conChooseChan = [];

for iFre = 1:size(Power_4con_8hz20hz,1)
	for iChan = 1:size(Power_4con_8hz20hz,2)
		cData = squeeze(Power_4con_8hz20hz(iFre,iChan,[1:size(Power_4con_8hz20hz,3)-1],:));
		PowerForSpss_4conChooseChan = [PowerForSpss_4conChooseChan cData];
	end
end

SNRForSpss_4conChooseChan = [];

for iFre = 1:size(SNR_4con_8hz20hz,1)
	for iChan = 1:size(SNR_4con_8hz20hz,2)
		cData = squeeze(SNR_4con_8hz20hz(iFre,iChan,[1:size(SNR_4con_8hz20hz,3)-1],:));
		SNRForSpss_4conChooseChan = [SNRForSpss_4conChooseChan cData];
	end
end

