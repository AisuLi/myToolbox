
%----------20180413, do AdaptiveRLS_foravg-------------

directory     = 'D:\expData\SSVEP_8loc\Exp1_20170328\Data\scanData\20180407\average';
f0            = [8,20];
chanNameList  = {'all'}; 
isDetrend     = 0;
orders        = 1; 
isAdaptiveRLS = 1;
condition     = {'allcue13_2~45Hz','allcue21_2~45Hz','allcue34_2~45Hz','allcue42_2~45Hz'};
subNum        = [1:21];


for iCondition = 1:size(condition,2)

		cd(fullfile(directory,condition{iCondition}));

		for iSub = subNum

				iSubFileName = [num2str(iSub),'-',condition{iCondition},'.avg'];

				fprintf('%s\n',condition{iCondition} );
				fprintf('%d\n',iSub );	
							
				for iFre = 1:size(f0,2)  
					
					[amp hsns hcns]  = AdaptiveRLS_foravg(iSubFileName,f0(iFre),chanNameList,isDetrend,orders,isAdaptiveRLS);					
					
					allAmp(:,:,iSub,iCondition,iFre)  = amp;  % size : 1-chan, 2-timeBin, 3-sub, 4-condition, 5-Fre
					% allHsns(:,:,iSub,iCondition,iFre) = hsns;
					% allHcns(:,:,iSub,iCondition,iFre) = hcns;
				end
				save('D:\expData\SSVEP_8loc\Exp1_20170328\Data\scanData\20180407\average\adapative.mat');
		end

		cd('..');
end 


%----------------shift conditions ----------------/
allAmp(:,:,:,:,1) = allAmp(:,:,:,[2 3 4 1],1);
allAmp(:,:,:,:,2) = allAmp(:,:,:,[3 2 1 4],2);
%-------------shift conditions, end---------------\



[P1eNo,standard_eName]=getElectrodeNo_bcl('haveM1',[]);
%--------minus base--------/
baseTimeRange   = [-500,-1];
xmin            = -1; 
baseTimePoint   = timeWinToPoint(xmin,baseTimeRange);
allAmpMinusBase = allAmp - repmat(mean(allAmp(:,[baseTimePoint(1):baseTimePoint(2)],:,:,:),2),[1 size(allAmp,2) 1 1 1]);
%------minus base,end------\

allAmpMinusBaseAvg = squeeze(mean(allAmpMinusBase,3));

allAmpMinusBaseAvg = allAmpMinusBaseAvg(:,[plotTimePoint(1):plotTimePoint(2)],:,:,:);

close all;
for iFig = 1:8 
	figure;
	set(gcf,'Position',get(0,'ScreenSize'),'color','w');
	i=0;
	for iChan = iFig*8-7:iFig*8
		i = i+1;
		subplot(4,2,i);
		H=plot([plotTimeRange(1):plotTimeRange(2)],squeeze(allAmpMinusBaseAvg(iChan,:,1,2)),'k','lineWidth',1);
		hold on;
		H=plot([plotTimeRange(1):plotTimeRange(2)],squeeze(allAmpMinusBaseAvg(iChan,:,2,2)),':k','lineWidth',1);
		H=plot([plotTimeRange(1):plotTimeRange(2)],squeeze(allAmpMinusBaseAvg(iChan,:,3,2)),'color',[0.3 0.3 0.3],'lineWidth',1);
		H=plot([plotTimeRange(1):plotTimeRange(2)],squeeze(allAmpMinusBaseAvg(iChan,:,4,2)),'r');
		title(standard_eName{iChan});
		box off;

		xlabel('time bins');
		ylabel('power (uV)');
		ylim([-0.1,0.3]);
		set(gca,'ytick',[-0.1:0.1:0.3]);
		xlim([-500,1200]);
		set(gca,'xtick',[-500,0,200,700,1200]);
		legend;
	end
end









close all;
allAmpMinusBaseAvg = squeeze(mean(allAmpMinusBase,3));
for iFig = 1:8 
	figure;
	set(gcf,'Position',get(0,'ScreenSize'),'color','w');
	i=0;
	for iChan = iFig*8-7:iFig*8
		i = i+1;
		subplot(4,2,i);
		H=plot([-1000:2000],squeeze(allAmpMinusBaseAvg(iChan,:,1,1)),'k','lineWidth',1);
		hold on;
		H=plot([-1000:2000],squeeze(allAmpMinusBaseAvg(iChan,:,2,1)),':k','lineWidth',1);
		H=plot([-1000:2000],squeeze(allAmpMinusBaseAvg(iChan,:,3,1)),'k','lineWidth',1);
		H=plot([-1000:2000],squeeze(allAmpMinusBaseAvg(iChan,:,4,1)),'r');
		title(standard_eName{iChan});
		box off;

		xlabel('time bins');
		ylabel('power (uV)');
		ylim([-0.2,1.4]);
		set(gca,'ytick',[-0.2:0.2:1.4]);
		xlim([-800,2000]);
		legend;
		% set(gca,'xtick',[-500,0,200,700,1200]);
	end
end













plotTimeRange    = [-500,1200];
plotTimePoint    = timeWinToPoint(xmin,plotTimeRange);

for iSub = 1:size(allAmpMinusBase,3)

	 for iFre = 1:size(allAmpMinusBase,5)

		maxChan_AmpMinusBase(:,:,iSub,:,iFre) = allAmpMinusBase(chooseChanSNR(iSub,:,iFre),:,iSub,:,iFre);
	 end
end



maxChan_AmpMinusBase = maxChan_AmpMinusBase(:,[plotTimePoint(1):plotTimePoint(2)],:,:,:);

maxChan_AmpMinusBase(6,:,:,:,:)  = mean(maxChan_AmpMinusBase,1);                                            
maxChan_AmpMinusBase(:,:,:,:,3)  = mean(maxChan_AmpMinusBase,5);
maxChan_AmpMinusBase(:,:,22,:,:) = mean(maxChan_AmpMinusBase,3);

H=plot([plotTimeRange(1):plotTimeRange(2)],squeeze(maxChan_AmpMinusBase(end,:,end,1,3)),'-b');
hold on;
H=plot([plotTimeRange(1):plotTimeRange(2)],squeeze(maxChan_AmpMinusBase(end,:,end,2,3)),':b');

H=plot([plotTimeRange(1):plotTimeRange(2)],squeeze(maxChan_AmpMinusBase(end,:,end,3,3)),'-r');

H=plot([plotTimeRange(1):plotTimeRange(2)],squeeze(maxChan_AmpMinusBase(end,:,end,4,3)),':r');
legend;