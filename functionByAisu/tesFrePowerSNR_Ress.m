
%%%%%%%    20180330  ress has only one chan, do fft of 4 conditions %%%%%%%%%%%%%%



% %%%-------read RESS data of each conditon of each fre of each subject ----------/


% directory     = 'D:\expData\SSVEP_8loc\Exp1_20170328\Data\scanData\20180325\cnt';
% allFile       = dir(fullfile(directory,'*BA3.eeg'));
% stimFre       = [8 20];
% cond2use      = [1:10;11:20;21:30;31:40];% 1:cue13; 2:cue21; 3:cue34; 4:cue42.
% caculateRange = [-500,1200];
% fftType       = 'avgFFT';


% for iSub = 1:size(allFile,1)

% 	eegName = allFile(iSub).name;
	
%  	for iCondition = 1:size(cond2use,1)
% 		[EEG]   = RESS_Yang_su(eegName,fftType,stimFre,cond2use(iCondition,:),caculateRange,0);

% 		allRess_ts(:,iSub,iCondition,1) = EEG.ssvep(1,1).ress_ts; % 8hz
% 		allRess_ts(:,iSub,iCondition,2) = EEG.ssvep(1,2).ress_ts; % 20hz
% 	end
% end

% save('RawRessEEG_4con'); %  do not do any shift 


load('RawRessEEG_4con'); 
% %----------- adjust the conditions-----------/
allRess_ts_merge =  cat(4,allRess_ts(:,:,[1 3 2 4],1),allRess_ts(:,:,[4 2 3 1],2)); 

% %---have adjusted the conditions,the end-----\

% %----it's the end, have read data of each conditon of each fre of each subject----\




%%%%%%%   do fft of wholeEpoch to see SNR %%%%%
timeBin      = [-500 1200];
freRange     = [5 25];
targetFre    = [8;20];

for iFre = 1:size(allRess_ts,4)

	for iCondition = 1:size(allRess_ts,3)
	
		for iSub = 1:size(allRess_ts,2)

			for iTimeBin = 1:size(timeBin,1)  % 
				
				timePoint(iTimeBin,:) = timeWinToPoint(EEG.xmin,timeBin(iTimeBin,:));	

				eegName = allRess_ts([timePoint(iTimeBin,1):timePoint(iTimeBin,2)],iSub,iCondition,iFre);

				[iSubAvg,freBins,Power,SNR] = eegfftSNR_bcl(eegName',freRange,[],[],fftType); 
				allPower_4con_wholeEpoch_ress(:,iSub,iCondition,iFre,iTimeBin)  = Power; % size: 1-freBins;  2-sub; 3-condition; 4-freStim; 5-timeRange(base,postBase)
				allSNR_4con_wholeEpoch_ress(:,iSub,iCondition,iFre,iTimeBin)    = SNR;
				clear  eegName Power SNR iSubAvg;
			end
		end
	end
end 

allSNR_4con_wholeEpoch_ress(:,size(allRess_ts,2)+1,:,:)   = mean(allSNR_4con_wholeEpoch_ress,2);
allPower_4con_wholeEpoch_ress(:,size(allRess_ts,2)+1,:,:) = mean(allPower_4con_wholeEpoch_ress,2);



targetFre = [8; 20]; 

%%-------find closest frePoint-------/
for iFre = 1:size(targetFre,1)
	[X, freIndex]         = sort(abs(freBins-targetFre(iFre,1)),'ascend');
	freTimePoint(iFre,1)  = freIndex(1);
	realFre(iFre,1) 	  = freBins(freIndex(1));
end
%%-------have found, the end--------\

subNumPerFig = 5;
figNum       = ceil(size(allSNR_4con_wholeEpoch_ress,2)./subNumPerFig);

%-------------------------plot  SNR & Power----------------------
for iPage = 1:figNum
	figure;
	set(gcf,'Position',get(0,'ScreenSize'),'color','w');
	suptitle(['RESS     left:    ',num2str(realFre(1,:)),'Hz     con1		con2   con3   con4     ',   num2str(timeBin(1,1)),'To',num2str(timeBin(1,2)),  '          right:  '   num2str(realFre(2,:))  'Hz     con1		con2   con3   con4     ']); 

	cSub = [1:subNumPerFig]+(iPage-1)*subNumPerFig;

	for iSub = cSub

		for iFre = 1:size(stimFre,2)

			for iCondition = 1:size(allSNR_4con_wholeEpoch_ress,3)
				
				if iSub<=size(allSNR_4con_wholeEpoch_ress,2)
					%-------plot topo-------/
					if iFre==1
						figureLoc = 9*(iSub-1-(iPage-1)*subNumPerFig)+iCondition;
					elseif iFre==2
						figureLoc = 9*(iSub-1-(iPage-1)*subNumPerFig)+iCondition+5;
					end	
					subplot(subNumPerFig,9,figureLoc);

					%%%%%%size: 1-freBins;  2-sub; 3-condition; 4-freStim; 
					[AX,H1,H2] = plotyy(freBins,allPower_4con_wholeEpoch_ress(:,iSub,iCondition,iFre),freBins,allSNR_4con_wholeEpoch_ress(:,iSub,iCondition,iFre),'plot');

					box off;

					set(AX(1),'XColor','k','YColor','b','Ylim',[0,0.1],'ytick',[0:0.05:0.1],'Xlim',freRange);
					set(AX(2),'XColor','k','YColor','r','Ylim',[0,4],'ytick',[0:1:3],'Xlim',freRange);
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
		if iSub<size(allPower_4con_wholeEpoch_ress,2)
				title(['sub:' num2str(iSub)]);
		elseif iSub==size(allPower_4con_wholeEpoch_ress,2)
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









%%%%%%%%%%%%%%%%%%%%%%%%%%% do fft of [-500,-1;201,700; 701 1200] %%%%%%%%%%%%%%%%%
%%%------ do fft of baseline and post-baseline for 4 conditions of 2 fre,separately -------/
timeBin      = [-500,-1;201,700; 701 1200];
freRange     = [5 25];
targetFre    = [8;20];

for iFre = 1:size(allRess_ts,4)

	for iCondition = 1:size(allRess_ts,3)
	
		for iSub = 1:size(allRess_ts,2)

			for iTimeBin = 1:size(timeBin,1)  % 
				
				timePoint(iTimeBin,:) = timeWinToPoint(EEG.xmin,timeBin(iTimeBin,:));	

				eegName = allRess_ts([timePoint(iTimeBin,1):timePoint(iTimeBin,2)],iSub,iCondition,iFre);

				[iSubAvg,freBins,Power,SNR] = eegfftSNR_bcl(eegName',freRange,[],[],fftType); 
				allPower_4con_ress(:,iSub,iCondition,iFre,iTimeBin)  = Power; % size: 1-freBins;  2-sub; 3-condition; 4-freStim; 5-timeRange(base,postBase)
				allSNR_4con_ress(:,iSub,iCondition,iFre,iTimeBin)    = SNR;
				clear  eegName Power SNR iSubAvg;
			end
		end
	end
end 
%%---do fft of baseline & post-baseline for 4 conditions of 2 fre,separately, the end-----\

allSNR_4con_ress(:,end+1,:,:,:)   = mean(allSNR_4con_ress,2);
allPower_4con_ress(:,end+1,:,:,:) = mean(allPower_4con_ress,2);





pptTitleName = 'RESS_SNR_2range';
mkdir(pptTitleName);
[opt,ppt] = figPPT(pptTitleName);

for iSub = 1:size(allSNR_4con_ress,2)

	figure;
    set (gcf,'Position',[50,50,1800,600], 'color','w');
	suptitle(['RESS     left:    ',num2str(realFre(1,:)),'Hz     con1		con2   con3   con4     ',     '          right:  '   num2str(realFre(2,:))  'Hz     con1		con2   con3   con4     ']); 

	for iTimeBin = 1:size(allSNR_4con_ress,5)
			
		for iFre = 1:size(stimFre,2)

			for iCondition = 1:size(allSNR_4con_ress,3)
				if iFre==1
					figureLoc = 9*(iTimeBin-1)+iCondition;
				elseif iFre==2
					figureLoc = 9*(iTimeBin-1)+iCondition+5;
				end	
				subplot(size(allSNR_4con_ress,5),9,figureLoc);

				%%%%%%size: 1-freBins;  2-sub; 3-condition; 4-freStim; 
				[AX,H1,H2] = plotyy(freBins,allPower_4con_ress(:,iSub,iCondition,iFre,iTimeBin),freBins,allSNR_4con_ress(:,iSub,iCondition,iFre,iTimeBin),'plot');

				box off;

				set(AX(1),'XColor','k','YColor','b','Ylim',[0,0.1],'ytick',[0:0.05:0.1],'Xlim',freRange);
				set(AX(2),'XColor','k','YColor','r','Ylim',[0,4],'ytick',[0:1:3],'Xlim',freRange);
				set(H1,'LineStyle','-','color','b');
				set(H2,'LineStyle','-','color','r');
				line(repmat(realFre(iFre,:),1,2),[0,.2],'linestyle',':','linewidth',2, 'color',[0 0 0 ]); %
			end
		end
	end

	if  iSub<size(allSNR_4con_ress,2)
			subStr = ['Sub:',num2str(iSub)];
	else
			subStr = 'allSubMean';
	end			
	cfigFileName = fullfile(fullfile(pwd,pptTitleName),[pptTitleName,'_',subStr([1:3,5:end]),'.bmp']);
	saveas(gcf,cfigFileName);
	ppt = ppt.addImageSlide(['RESS_',subStr(1:3),'_',subStr(5:end)],cfigFileName,opt);
	close all;
end


ppt.saveAs(fullfile(pwd,pptTitleName));

ppt.close;

















%%-------find closest frePoint-------/
for iFre = 1:size(targetFre,1)
	[X, freIndex]         = sort(abs(freBins-targetFre(iFre,1)),'ascend');
	freTimePoint(iFre,1)  = freIndex(1);
	realFre(iFre,1) 	  = freBins(freIndex(1));
end
%%-------have found, the end--------\


%%-----divide base------/
allSNRDivideBase_4con_ress           = squeeze(allSNR_4con_ress(:,:,:,:,2)./allSNR_4con_ress(:,:,:,:,1)); % size: size: 1-freBins;  2-sub; 3-condition; 4-freStim; 
allFrePointSNRDivideBase_4con_ress   = cat(1,allSNRDivideBase_4con_ress(freTimePoint(1),:,:,1),allSNRDivideBase_4con_ress(freTimePoint(2),:,:,2)); % 1-freBins;  2-sub; 3-condition
allPowerDivideBase_4con_ress         = squeeze(allPower_4con_ress(:,:,:,:,2)./allPower_4con_ress(:,:,:,:,1)); % size: size: 1-freBins;  2-sub; 3-condition; 4-freStim; 
allFrePointPowerDivideBase_4con_ress = cat(1,allPowerDivideBase_4con_ress(freTimePoint(1),:,:,1),allPowerDivideBase_4con_ress(freTimePoint(2),:,:,2)); % 1-freBins;  2-sub; 3-condition
%%---divide base,end---\



allFrePointSNRDivideBase_4con_ress(3,:,:)                                      = mean(allFrePointSNRDivideBase_4con_ress,1); 
allFrePointSNRDivideBase_4con_ress(:,size(allSNRDivideBase_4con_ress,2)+1,:)   = mean(allFrePointSNRDivideBase_4con_ress,2);
allFrePointPowerDivideBase_4con_ress(3,:,:)                                    = mean(allFrePointPowerDivideBase_4con_ress,1); 
allFrePointPowerDivideBase_4con_ress(:,size(allSNRDivideBase_4con_ress,2)+1,:) = mean(allFrePointPowerDivideBase_4con_ress,2);

%%--------compute sem ---------/
for iFre = 1:size(allFrePointSNRDivideBase_4con_ress,1)
	for iCondition = 1:size(allFrePointSNRDivideBase_4con_ress,3)

		[gm,gsem] = grpstats(squeeze(allFrePointSNRDivideBase_4con_ress(iFre,[1:size(allSNRDivideBase_4con_ress,2)],iCondition))',{},{'mean','sem'});
		allFrePointSNRDivideBase_4con_ress_sem(iFre,iCondition) = gsem;
		clear gm gsem
		[gm,gsem] = grpstats(squeeze(allFrePointPowerDivideBase_4con_ress(iFre,[1:size(allPowerDivideBase_4con_ress,2)],iCondition))',{},{'mean','sem'});
		allFrePointPowerDivideBase_4con_ress_sem(iFre,:,iCondition) = gsem;
		clear gm gsem
	end
end
%%-----compute sem,the end-----\

%%%%%%  start to plot bar of SNR & Power  %%%%%%%%
barColor               = repmat([0.85,1,0.1,0.3]',1,3);
ngroups                = size(allFrePointSNRDivideBase_4con_ress,1);
nbars                  = size(allFrePointSNRDivideBase_4con_ress,3);
groupwidth             = min(0.8, nbars/(nbars+1.5));	

%%% ------- SNR start ------------/
plotSNRRange           = [0 4; 0 4];
figure;
set(gcf,'Position',get(0,'ScreenSize'),'color','w');
set(gcf,'color','white')
suptitle('SNR');
for iSub = 1:size(allFrePointSNRDivideBase_4con_ress,2)
		subplot(5,5,iSub);

		cSNR = squeeze(allFrePointSNRDivideBase_4con_ress(:,iSub,:));

		H  = bar(cSNR,'grouped');
		ch = get(H,'children');
        set(ch{1},'Facecolor',barColor(1,:));
        set(ch{2},'Facecolor',barColor(2,:));
        set(ch{3},'Facecolor',barColor(3,:));
        set(ch{4},'Facecolor',barColor(4,:));
		hold on;box off;
		set(gca,'XTickLabel',{'8hz','20hz','mean'});
		subStr = ['Sub',num2str(iSub)];
		ylim(plotSNRRange(1,:));

		if iSub==size(allSNRDivideBase_4con_ress,2)+1  
			subStr = 'allSubMean';
			cError = allFrePointSNRDivideBase_4con_ress_sem;   %%% size: 1-freBins; 2-chan; 3-condition;
			for i = 1:nbars
			    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
			    errorbar(x,cSNR(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
			    hold on;
				ylim(plotSNRRange(2,:));
			end
		end
		title([subStr]);
end 
l = legend([ch{1} ch{2} ch{3} ch{4} ],'sepCue--Uncue','ajCue--Uncue','ajCue--Cue','SepCue--Cue');
legend('boxoff');
set(l,'Position',[0.434548611111108 0.142414860681115 0.0765625 0.0818713450292398]);
%%% ------- SNR, the end -----------\



%%% ------- Power start ------------/
plotPowerRange           = [0 4; 0 4];
figure;
set(gcf,'Position',get(0,'ScreenSize'),'color','w');
set(gcf,'color','white')
suptitle('Power');
for iSub = 1:size(allFrePointPowerDivideBase_4con_ress,2)
		subplot(5,5,iSub);

		cPower = squeeze(allFrePointPowerDivideBase_4con_ress(:,iSub,:));

		H  = bar(cPower,'grouped');
		ch = get(H,'children');
        set(ch{1},'Facecolor',barColor(1,:));
        set(ch{2},'Facecolor',barColor(2,:));
        set(ch{3},'Facecolor',barColor(3,:));
        set(ch{4},'Facecolor',barColor(4,:));
		hold on;box off;
		set(gca,'XTickLabel',{'8hz','20hz','mean '});
		subStr = ['Sub',num2str(iSub)];
		ylim(plotPowerRange(1,:));

		if iSub==size(allPowerDivideBase_4con_ress,2)+1  
			subStr = 'allSubMean';
			cError = allFrePointPowerDivideBase_4con_ress_sem;   %%% size: 1-freBins; 2-chan; 3-condition;
			for i = 1:nbars
			    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
			    errorbar(x,cPower(:,i),cError(:,i),'LineStyle','none','color',[0 0 0],'linewidth',1);
			    hold on;
				ylim(plotPowerRange(2,:));
			end
		end
		title([subStr]);
end 
l = legend([ch{1} ch{2} ch{3} ch{4} ],'sepCue--Uncue','ajCue--Uncue','ajCue--Cue','SepCue--Cue');
legend('boxoff');
set(l,'Position',[0.434548611111108 0.142414860681115 0.0765625 0.0818713450292398]);
%%% ------- Power, the end -----------\




% SNRForSpss_ress = [];

% for iFre = 1:size(allFrePointPowerDivideBase_4con_ress,1)

% 		cData = squeeze(allFrePointPowerDivideBase_4con_ress(iFre,[1:size(allFrePointPowerDivideBase_4con_ress,2)-1],:));
% 		SNRForSpss_ress = [SNRForSpss_ress cData];

% end



% [H,P_allavg(1,1),CI] = ttest(allavgCutBase(:,2,1),allavgCutBase(:,2,3));


% PowerForSpss_ress = [];

% for iFre = 1:size(allFrePointPowerDivideBase_4con_ress,1)

% 		cData = squeeze(allFrePointPowerDivideBase_4con_ress(iFre,[1:size(allFrePointPowerDivideBase_4con_ress,2)-1],:));
% 		PowerForSpss_ress = [PowerForSpss_ress cData];

% end




























% directory     = 'D:\expData\SSVEP_8loc\Exp1_20170328\Data\scanData\20180201\cnt';
% allFile       = dir(fullfile(directory,'*BA3.eeg'));
% stimTrigger = [1:10;11:20;21:30;31:40];

% for iSub = 1:size(allFile,1)

% 	eegName = allFile(iSub).name;

% 	[signal, accept, typeeeg, rt, response, chan_names, pnts, ...
% 	ntrials, srate, xmin, xmax] = loadeeg(eegName);

% 	trialNum(iSub,1) =  ntrials;

% 		for iCondition = 1:size(stimTrigger,1)
% 			trialNum(iSub,iCondition+1) =  sum(ismember(typeeeg,stimTrigger(iCondition,:)));		

% 		end
% end 