
clear all;

[myColormap] = makeColormap(30,0,1,0);	

monitorNum   = 1; % 1,2,3,4 for ulmb ,od, 278,CRT

onoffsetBin  = [34 50; 34 116; 34,116; 35,40]; % ulmb ,od, 278,CRT

 

nframes    = 4;
ifi        = 1000/120;	
cntNum     = [1:9];

cDirectory = pwd;
cntFiles   = dir(fullfile(cDirectory,'*.cnt'));


figure;
set (gcf,'Position',[100,100,800,600], 'color','w');
iLoc      = 0;


% randColor    = rand(50,3);

rectColor    = [255,255,102]/255;
transParency = 0.4;

for i = cntNum
	icntFile          = cntFiles(i).name;
	underline         = strfind(icntFile,'_');
	monitorStrNum     = underline(1)+1:underline(2)-1;
	modeStrNum        = underline(2)+1:underline(2)+2;
	stimulsTypeStrNum = underline(3)+1:underline(3)+4;
	locStrNum         = underline(4)+1:underline(4)+3;
	iLoc              = iLoc+1;
	[eventBins,foundMarkers,epochedTriggerBins,beplotedX,epochedData,randColor] = monitorTest2(icntFile,nframes,ifi,iLoc);

	alleventBins(iLoc,:)          = eventBins;
	allFoundMarkers(iLoc,:)       = foundMarkers;
	allepochedTriggerBins(iLoc,:) = epochedTriggerBins;
	allbeplotedX(iLoc,:)          = beplotedX;
	allepochedData(:,:,iLoc)      = epochedData;
	hold on;
	cMeanEpochedData              = mean(epochedData,2);

YellowRect = [20 50 50 20;8 31 31 8;12 42 42 12;50 85 85 50];

	for iRect = 1:3
		frameBin(iRect,:) = [onoffsetBin(monitorNum,:)]+round(ifi*(iRect-1)*10);
		rectX(iRect,:)    = [frameBin(iRect,1) frameBin(iRect,1) frameBin(iRect,2) frameBin(iRect,2)]-(abs(beplotedX(1))+1);
		iRectFrameMaxMin  = [max(cMeanEpochedData(frameBin(iRect,1):frameBin(iRect,2))),min(cMeanEpochedData(frameBin(iRect,1):frameBin(iRect,2)))];
		% rectY(iRect,:)    = [iRectFrameMaxMin(2)-5 iRectFrameMaxMin(1)+5 iRectFrameMaxMin(1)+5 iRectFrameMaxMin(2)-5];
		rectY(iRect,:)    = YellowRect(monitorNum,:); %  ulmb:[20 50 50 20]; od: [8 31 31 8] 278:[12 42 42 12]; crt:[50 85 85 50];

	end
		
	rectX  = rectX';
	rectY  = rectY'/85.175;
	H_rect = patch(rectX,rectY,rectColor);

	% patch([Xleft, Xleft, Xright, Xright], [Ydown, Yup, Yup, Ydown],color);

	set(H_rect,'EdgeColor',rectColor,'EdgeAlpha',transParency,'FaceAlpha',transParency);

	titleName = [icntFile(modeStrNum),'-',icntFile(stimulsTypeStrNum),'-',icntFile(locStrNum)];
    title(titleName);
	set(gcf,'color','w');
	if i<9
		clear cMeanEpochedData frameBin rectX iRectFrameMaxMin  rectY;
	end
end

figureName = [icntFile(monitorStrNum),'-',icntFile(modeStrNum),'-',icntFile(stimulsTypeStrNum)];

print('-depsc','-painters',fullfile(cDirectory,[figureName,'_fig1.eps']));
saveas(gcf,[figureName,'_fig1_1207.bmp']); 




meanEpochedData        = squeeze(mean(allepochedData,2));

meanEpochedTriggerBins = mean(allepochedTriggerBins,2);

% save([figureName,'.mat'],'alleventBins','allFoundMarkers','allbeplotedX','allepochedData','allepochedTriggerBins','meanEpochedData','meanEpochedTriggerBins','ifi','figureName');



% firstType: Onset to offset;



firstFrameData = mean(meanEpochedData([frameBin(1,1):frameBin(1,2)],:),1);	
frame2Data     = mean(meanEpochedData([frameBin(2,1):frameBin(2,2)],:),1);
frame3Data     = mean(meanEpochedData([frameBin(3,1):frameBin(3,2)],:),1);
midFrameData   = mean([frame2Data;frame3Data],1);
framePercent   = firstFrameData./midFrameData;

%%% max of CRT %%%%%


frame1max       = squeeze(max(allepochedData([34:34+82],:,:)));
frame2max       = squeeze(max(allepochedData([34+83:34+82+83],:,:)));
frame3max       = squeeze(max(allepochedData([34+83*2:34+82+83*2],:,:)));
mid2frameMax    = (frame2max+frame3max)/2;
framePercentMax = mean(frame1max./mid2frameMax);



for i=1:3;
	nineLocFramePercentMax(i,:)     = framePercentMax([3*i-2:3*i]);
	nineLocMeanEpochedData(i,[1:3]) = mean(meanEpochedData([34:end],[3*i-2:3*i]),1);
	delay(i,:)                      = meanEpochedTriggerBins([3*i-2:3*i])'/10;
	nineLocFramePercent(i,:)        = framePercent([3*i-2:3*i]);
end



if monitorNum==4
	figure;
	imagesc(nineLocFramePercentMax);
	box off;
	colorbar
	caxis([0.7,1]);
	colormap(myColormap);	
	title([figureName,'-framePercent']);
	print('-dpdf','-painters',[figureName,'_fig4framePercent_max.pdf']);
	saveas(gcf,[figureName,'_fig4max.bmp']);
end


nineLocMeanEpochedData = nineLocMeanEpochedData/nineLocMeanEpochedData(2,2);
nineLocMeanDelay       = delay/delay(2,2);

figure;
imagesc(abs(delay));
box off;
colorbar;
caxis([0,20]);
colormap(myColormap);
title([figureName,'-delay']);
print('-dpdf','-painters',[figureName,'_fig2delay.pdf']);
saveas(gcf,[figureName,'_fig2.bmp']);

figure;
imagesc(nineLocMeanEpochedData);
box off;
colorbar;
caxis([0.7,1.1]);
colormap(myColormap);	
title([figureName,'-allmeanVoltage']);
print('-dpdf','-painters',[figureName,'_fig3allmeanVoltage.pdf']);
saveas(gcf,[figureName,'_fig3.bmp']);



figure;
imagesc(nineLocFramePercent);
box off;
colorbar
caxis([0.7,1.1]);
colormap(myColormap);	
title([figureName,'-framePercent']);
print('-dpdf','-painters',[figureName,'_fig4framePercent.pdf']);
saveas(gcf,[figureName,'_fig4.bmp']);

save([figureName,'.mat']);





save([figureName,'.mat']);


% % %%%%%%%%%%%%%%%%%%%%%%%%%%  luminance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all;
% load('luminance.mat');
% allMode      = {'278Q-ulmb','278Q-OD','278','CRT'};

% [myColormap] = makeColormap(30,1);	

% % %%%%% spatial consistency %%%%
% for imode = 1:4
% 	figure;
% 	imodeValue = luminance(:,:,imode);
% 	verPer(:,:,imode)  = [reshape(imodeValue(:,1)./imodeValue(:,2),3,3)]';
% 	oblqPer(:,:,imode) = [reshape(imodeValue(:,3)./imodeValue(:,2),3,3)]';
% 	imagesc(verPer(:,:,imode));
% 	box off;
% 	colorbar;
% 	caxis([0.7,1.1]);
% 	colormap(myColormap);
% 	title([allMode{imode},'-ver']);
% 	% print('-depsc','-painters',fullfile(pwd,[allMode{imode},'_verPer.eps']));
% 	print('-dpdf','-painters',fullfile(pwd,[allMode{imode},'_verPer.pdf']));

% 	figure;
% 	imagesc(oblqPer(:,:,imode));
% 	box off;
% 	colorbar;	
% 	caxis([0.7,1.1]);
% 	colormap(myColormap);	
% 	title([allMode{imode},'-oblq']);
% 	% print('-depsc','-painters',fullfile(pwd,[allMode{imode},'_oblqPer.eps']));
% 	print('-dpdf','-painters',fullfile(pwd,[allMode{imode},'_oblqPer.pdf']));
% end

% %%%%%%%%%%%%% bar %%%%%%%%%%%
% clear all;
% load('luminance.mat');
% ylimRange = [0.7,1.1];
% allMode   = {'278Q-ulmb','278Q-OD','278','CRT'};

% for imode = 1:4
% 	imodeValue = luminance(:,:,imode);
% 	imodePer(:,1,imode) = imodeValue(:,1)./imodeValue(:,2);
% 	imodePer(:,2,imode) = imodeValue(:,3)./imodeValue(:,2);	

% 	figure;
% 	for iLoc=1:9
% 			subplot(3,3,iLoc);
% 			bar(imodePer(iLoc,:,imode));
% 			ylim(ylimRange);
% 			set(gca,'ytick',[ylimRange(1):0.1:ylimRange(2)]);
% 			set(gca,'XTickLabel',{'ver','oblq'});
% 			hold on;
% 			box off;
% 	end
% 	monitorMode = allMode{imode};
% 	title(monitorMode);
% 	print('-depsc','-painters',fullfile(pwd,[monitorMode,'_lumiPer.eps']));
% end

% figure;
% for iLoc=1:9
% 	subplot(3,3,iLoc);
% 	bar(luminance(iLoc,[1:3]));
% 	ylim([10,70]);
% 	set(gca,'ytick',[10:20:70]);
% 	set(gca,'XTickLabel',{'ver','hor','obliq'});
% 	hold on;
% 	box off;
% end
% monitorMode='278Q-ulmb';
% title(monitorMode);
% print('-depsc','-painters',fullfile(pwd,[monitorMode,'_luminance.eps']));



% figure;
% for iLoc=1:9
% 	subplot(3,3,iLoc);
% 	bar(luminance(iLoc,[5:7]));
% 	ylim([120,180]);
% 	set(gca,'ytick',[120:20:180]);
% 	set(gca,'XTickLabel',{'ver','hor','obliq'});
% 	hold on;
% 	box off;
% end
% monitorMode='278Q-OD';
% title(monitorMode);
% print('-depsc','-painters',fullfile(pwd,[monitorMode,'_luminance.eps']));



% figure;
% for iLoc=1:9
% 	subplot(3,3,iLoc);
% 	bar(luminance(iLoc,[9:11]));
% 	ylim([110,170]);
% 	set(gca,'ytick',[110:20:170]);
% 	set(gca,'XTickLabel',{'ver','hor','obliq'});
% 	hold on;
% 	box off;
% end
% monitorMode='278-stan';
% title(monitorMode);
% print('-depsc','-painters',fullfile(pwd,[monitorMode,'_luminance.eps']));


% figure;
% for iLoc=1:9
% 	subplot(3,3,iLoc);
% 	bar(luminance(iLoc,[13:15]));
% 	ylim([0,60]);
% 	set(gca,'ytick',[0:20:60]);
% 	set(gca,'XTickLabel',{'ver','hor','obliq'});
% 	hold on;
% 	box off;
% end
% monitorMode='CRT';
% title(monitorMode);
% print('-depsc','-painters',fullfile(pwd,[monitorMode,'_luminance.eps']));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%