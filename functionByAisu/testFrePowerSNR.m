

load('a.mat');
clearvars -except chooseChanSNR chooseChanPower;

directory    = 'D:\expData\SSVEP_8loc\Exp1_20170328\Data\scanData\20180201\cnt';
allFile      = dir(fullfile(directory,'*BA2.eeg'));
freRange     = [4 12;4 12]; % 16 24  4 12
timeBin      = [-500 -1;201 1200]; % 
stimTrigger  = [1:10,31:40;11:30];  



isSortchan   = false; % if true ,sort maxPower automatically; if false,choose given chans 

if freRange(1,2)<15
	stimFre(1,1)   = [8];
else
	stimFre(1,1)   = [20];
end

if freRange(2,1)>15
	stimFre(1,2)   = [20];
else
	stimFre(1,2)   = [8];
end



% if triggerLoc(:)<11
% 	cueCon = ['cueLoc1&3  compareLoc2:  ',num2str(stimFre(1)),'Hz'];
% elseif triggerLoc(:)>30
% 	cueCon = ['cueLoc2&4  compareLoc3:  ',num2str(stimFre(1)),'Hz'];   
% end	   


fftType      = 'avgFFT'; % 'trialFFT'


cd(directory);

for iSub = 1:size(allFile,1)

	eegName                                  = allFile(iSub).name;

	[iSubAvg,freBins_8Hz,Powermean,SNRmean]  = eegfftSNR_bcl(eegName,freRange(1,:),timeBin(1,:),stimTrigger(1,:),fftType);
	allPower_8Hz(:,:,iSub)                   = Powermean;
	allSNR_8Hz(:,:,iSub)                     = SNRmean;
	avgData_8Hz(:,:,iSub)					 = iSubAvg;
	clear  Powermean SNRmean iSubAvg;
	
	[iSubAvg,freBins_20Hz,Powermean,SNRmean] = eegfftSNR_bcl(eegName,freRange(2,:),timeBin(2,:),stimTrigger(1,:),fftType);
	allPower_20Hz(:,:,iSub)                  = Powermean;
	allSNR_20Hz(:,:,iSub)                    = SNRmean;
	avgData_20Hz(:,:,iSub)					 = iSubAvg;
	clear  Powermean SNRmean iSubAvg;
end


if strcmp(fftType,'avgFFT')
	[avgData,freBins,allPower_8Hz(:,:,size(allPower_8Hz,3)+1),allSNR_8Hz(:,:,size(allSNR_8Hz,3)+1)]     = eegfftSNR_bcl(avgData_8Hz,freRange(1,:),timeBin(1,:),stimTrigger(:),fftType);
	[avgData,freBins,allPower_20Hz(:,:,size(allPower_20Hz,3)+1),allSNR_20Hz(:,:,size(allSNR_20Hz,3)+1)] = eegfftSNR_bcl(avgData_20Hz,freRange(2,:),timeBin(2,:),stimTrigger(:),fftType);
	clear avgData;
end


a = min([size(allSNR_8Hz,1) size(allSNR_20Hz,1)]); % keep the number of frebins is equal
allSNR(:,:,:,1)   =  allSNR_8Hz([1:a],:,:);
allSNR(:,:,:,2)   =  allSNR_20Hz([1:a],:,:);
allPower(:,:,:,1) =  allPower_8Hz([1:a],:,:);
allPower(:,:,:,2) =  allPower_20Hz([1:a],:,:);
freBins           =  [freBins_8Hz([1:a]);freBins_20Hz([1:a])];


[none,standard_eName] = getElectrodeNo_bcl('noM1',[]);

for istim = 1:size(stimFre,2)

	%-------find closest 1 or 2 or 3 frePoint-------/
	% eg. targetFre=8hz, frebins=[...,7.95 8.03,]
	[X, freIndex]         = sort(abs(freBins(istim,:)-stimFre(istim)),'ascend');
	freTimePoint(:,istim) = freIndex(1)';
	realFre(istim,:) 	  = freBins(istim,freIndex(1));
	%----------------------------------------------\
end

frePointSNR(1,:,:) = allSNR_8Hz(freTimePoint(:,1),:,:); % 8Hz
frePointSNR(2,:,:) = allSNR_20Hz(freTimePoint(:,2),:,:);% 20Hz


frePointPower(1,:,:) = allPower_8Hz(freTimePoint(:,1),:,:); % 8Hz
frePointPower(2,:,:) = allPower_20Hz(freTimePoint(:,2),:,:);% 20Hz



waveNum = 3;

[myColormap] = makeColormap(20,0,1,'r');

subNumPerFig = 5;
figNum       = ceil(size(allFile,1)./subNumPerFig);


%-------------------------plot topo of SNR----------------------
for iPage = 1:figNum
	figure;
	set(gcf,'Position',get(0,'ScreenSize'),'color','w');
	% suptitle([fftType,'    ',num2str(realFre(1,:)),'Hz      ',num2str(timeBin(1,1)),'To',num2str(timeBin(1,2)),'         SNR-topo:   left:   sepCue           right: ajCue ' ]);
	% suptitle([fftType,'      ',num2str(timeBin(1,1)),'To',num2str(timeBin(1,2)),'      AllCondition   SNR-topo:   left: ',num2str(realFre(1,:)),'Hz           right: ',num2str(realFre(2,:)),'Hz']);
	% suptitle([fftType, num2str(timeBin(1,1)),'To',num2str(timeBin(1,2)), '     SNR-topo left: cue1&3: ', num2str(realFre(1,:)), '          right:   cue2&4  ',num2str(realFre(2,:))]); 
	suptitle([fftType,'   ',num2str(realFre(1,:)),' Cue24    SNR-topo left: ',   num2str(timeBin(1,1)),'To',num2str(timeBin(1,2)),  '          right:    ',num2str(timeBin(2,1)),'To',num2str(timeBin(2,2))]); 

	cSub = [1:subNumPerFig]+(iPage-1)*subNumPerFig;

	for iSub = cSub
		for istim = 1:size(stimFre,2)

			if iSub<=size(allSNR_8Hz,3)
				%-------plot topo-------/
				if istim==1
					figureLoc = 2*(waveNum+1)*(iSub-1-(iPage-1)*subNumPerFig)+1;
				elseif istim==2
					figureLoc = 2*(waveNum+1)*(iSub-1-(iPage-1)*subNumPerFig)+waveNum+2;
				end	
				subplot(subNumPerFig,2*(waveNum+1),figureLoc);
				topoplot_bcl(frePointSNR(istim,:,iSub),'EEGChans',standard_eName,'maplimitsDouble',[0 5],'colormap',myColormap);
				%-----------------------\
				if isSortchan
					[X, chanIndex] = sort(frePointSNR(istim,[1:size(standard_eName,2)-2],iSub),'descend');						
					chooseChanSNR(iSub,:,istim)= chanIndex(1:3);	
				else
					chanIndex = chooseChanSNR(iSub,:,1);
				end


				for iWave = 1:waveNum
					figureLoc  = figureLoc+1;
					subplot(subNumPerFig,2*(waveNum+1),figureLoc);
					[AX,H1,H2] = plotyy(freBins(istim,:),allPower(:,chanIndex(iWave),iSub,istim),freBins(istim,[3:end-2]),allSNR([3:end-2],chanIndex(iWave),iSub,istim),'plot');
					if iSub<size(allSNR_8Hz,3)
						title([num2str(iSub),standard_eName(chanIndex(iWave))]);
					else
						title(['allSubAvg',standard_eName(chanIndex(iWave))]);
					end
					box off;
					if  istim ==1
						set(AX(1),'XColor','k','YColor','b','Ylim',[0,0.3],'ytick',[0:0.1:0.3],'Xlim',freRange(1,:));
					elseif istim ==2
						set(AX(1),'XColor','k','YColor','b','Ylim',[0,0.3],'ytick',[0:0.1:0.3],'Xlim',freRange(2,:));						
					end
					
					set(AX(2),'XColor','k','YColor','r','Ylim',[0,5],'ytick',[0:1:5]);
					set(H1,'LineStyle','-','color','b');
					set(H2,'LineStyle','-','color','r');
					line(repmat(realFre(istim,:),1,2),[0,.2],'linestyle',':','linewidth',2, 'color',[0 0 0 ]); %
				end
			end
		end
	end
	legend([H1,H2],{'power';'SNR'});
	legend('boxoff');
	set(gcf,'color','w');
	% saveas(gcf,['S',num2str(cSub(1)),'~',num2str(cSub(end)),'_',fftType,'_SNRtopo.pdf']);
end
%------------------------------\\




%-------------------------plot topo of power----------------------


[myColormap] = makeColormap(20,0,1,'b');
for iPage = 1:figNum
	figure;
	set(gcf,'Position',get(0,'ScreenSize'),'color','w');
	suptitle([fftType,num2str(realFre(1,:)),' cue1&3    SNR-topo left: ',   num2str(timeBin(1,1)),'To',num2str(timeBin(1,2)),  '          right:    ',num2str(timeBin(1,1)),'To',num2str(timeBin(1,2))]); 

% 	suptitle([fftType,'     ',num2str(realFre(1,:)),'Hz     ',num2str(timeBin(1,1)),'To',num2str(timeBin(1,2)),'         Power-topo:   left:   sepCue           right: ajCue ' ]);
	% suptitle([fftType,'      ',num2str(timeBin(1,1)),'To',num2str(timeBin(1,2)),'      AllCondition   Power-topo:   left: ',num2str(realFre(1,:)),'Hz           right: ',num2str(realFre(2,:)),'Hz']);
% 	suptitle([fftType,num2str(realFre(1,:)),'Hz        Power-topo',cueCon,'  left:    ', num2str(timeBin(1,1)),'To',num2str(timeBin(1,2)),'          right:    ',num2str(timeBin(2,1)),'To',num2str(timeBin(2,2)),]);
	% suptitle([fftType, '   ',num2str(timeBin(1,1)),'To',num2str(timeBin(1,2)), '     Power-topo left: cue1&3: ', num2str(realFre(1,:)), '          right:   cue2&4  ',num2str(realFre(2,:))]); 
	suptitle([fftType,'   ',num2str(realFre(1,:)),' Cue24    Power-topo left: ',   num2str(timeBin(1,1)),'To',num2str(timeBin(1,2)),  '          right:    ',num2str(timeBin(2,1)),'To',num2str(timeBin(2,2))]); 

	cSub = [1:subNumPerFig]+(iPage-1)*subNumPerFig;
	for iSub = cSub
		for istim = 1:size(stimFre,2)
			if iSub<=size(allSNR_8Hz,3)
				%-------plot topo-------/
				if istim==1
					figureLoc = 2*(waveNum+1)*(iSub-1-(iPage-1)*subNumPerFig)+1;
				elseif istim==2
					figureLoc = 2*(waveNum+1)*(iSub-1-(iPage-1)*subNumPerFig)+waveNum+2;
				end	
				subplot(subNumPerFig,2*(waveNum+1),figureLoc);
				topoplot_bcl(frePointPower(istim,:,iSub),'EEGChans',standard_eName,'maplimitsDouble',[0 0.2],'colormap',myColormap);
				%-----------------------\
				if isSortchan
					[X, chanIndex] = sort(frePointPower(istim,[1:size(standard_eName,2)-2],iSub),'descend');						
					chooseChanPower(iSub,:,istim)= chanIndex(1:3);		
				else			
					chanIndex = chooseChanPower(iSub,:,1);
				end
					
				for iWave = 1:waveNum
					figureLoc  = figureLoc+1;
					subplot(subNumPerFig,2*(waveNum+1),figureLoc);
					[AX,H1,H2] = plotyy(freBins(istim,:),allPower(:,chanIndex(iWave),iSub,istim),freBins(istim,[3:end-2]),allSNR([3:end-2],chanIndex(iWave),iSub,istim),'plot');
					
					if iSub<size(allSNR_8Hz,3)
						title([num2str(iSub),standard_eName(chanIndex(iWave))]);
					else
						title(['allSubAvg',standard_eName(chanIndex(iWave))]);
					end

					box off;

					if istim ==1
						set(AX(1),'XColor','k','YColor','b','Ylim',[0,0.3],'ytick',[0:0.1:0.3],'Xlim',freRange(1,:));
					elseif istim==2
						set(AX(1),'XColor','k','YColor','b','Ylim',[0,0.3],'ytick',[0:0.1:0.3],'Xlim',freRange(2,:));
					end

					set(AX(2),'XColor','k','YColor','r','Ylim',[0,5],'ytick',[0:1:5]);
					set(H1,'LineStyle','-','color','b');
					set(H2,'LineStyle','-','color','r');
					line(repmat(realFre(istim,:),1,2),[0,.2],'linestyle',':','linewidth',2, 'color',[0 0 0 ]); %
				end
			end
		end
	end
	legend([H1,H2],{'power';'SNR'});
	legend('boxoff');
	set(gcf,'color','w');
	% saveas(gcf,['S',num2str(cSub(1)),'~',num2str(cSub(end)),'_',fftType,'_Powertopo.pdf']);
end
%------------------------------\\















% %--------------------- bbbbkkk--------------------------

% % for iPage = 1:3
% % 	figure;
% % 	ifigure = 0;
% % 	for iSub = [1:5]+(iPage-1)*5
% % 		for istim = 1:2
% % 				ifigure=ifigure+1;
% % 				subplot(5,6,ifigure);
% % 				topoplot_bcl(allSNR(3*(istim-1)+iSNR,:,iSub),'EEGChans',standard_eName,'maplimitsDouble',[0 3]);
% % 			end
% % 		end
% % 	end
% % end
% % %------------------------------\\

% % % %------ plot each subject's power & SNR ----------------/

% % for ie = 1:size(allY_20Hz,2)

% % 	set (gcf,'Position',get(0,'ScreenSize'), 'color','w');
% % 	suptitle([standard_eName(ie),'-20Hz']);
% % 	for isub = 1:size(allY_20Hz,3)
% % 		subplot(3,5,isub);
% % 		[AX,H1,H2] = plotyy(freBins_20Hz,allY_20Hz(:,ie,isub),freBins_20Hz(3:end-2),allSNR_20Hz([3:end-2],ie,isub),'plot');
% % 		set(AX(1),'XColor','k','YColor','b','Ylim',[0,0.2],'ytick',[0:0.05:0.2]);
% % 		set(AX(2),'XColor','k','YColor','r','Ylim',[0.5,4.5],'ytick',[0.5:1:4.5]);
% % 		set(H1,'LineStyle','-','color','b');
% % 		set(H2,'LineStyle','-','color','r');
% % 		line([20 20],[0,.1],'linestyle',':','linewidth',1, 'color',[0 0 0 ]); %
% %     end
% %     hold on;
% % 	legend([H1,H2],{'power';'SNR'});
% % 	KbWait;
% % 	close 
% % end
% % %-----------------------------------------------------\
