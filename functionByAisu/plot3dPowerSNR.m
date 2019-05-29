
function  plot3dPowerSNR(nSub,fre,cntfilename,isEpoched,moreBinNum,stimTrigger,freRange,fftType,ldrType)


%% written by Li Aisu, 08-31-2017
%% it is used for comparing the power and SNR of different stimulus

if ~exist('nSub','var')
	error(' how many subjects ?');
end 

% nSub = 7;
% fre = [8 12 15 120/7];
% cntfilename = {'white-2frame.cnt','white-halframe.cnt','checker-2frame.cnt','checker-halframe.cnt','counterphase-2frame.cnt','counterphase-halframe.cnt'}; % 
% isEpoched = 0;
% moreBinNum = 0;
% stimTrigger = [1:40];
% freRange = [5, 20];
% fftType = 'cntFFT';
% ldrType = 'noLDR';


if 	~isEpoched
	fftType = 'cntFFT';
end

if ~exist('fre','var')
	error(' which  Frequencies ?');
end 



if ~exist('freRange','var')
	freRange = [5,20];
end

if ~exist('ldrType','var') || strcmpi('noLDR',ldrType)
    
    ldrType='noLDR';
	eType = 'noM1';
	eName = {'C5','C3','C1','CZ','C2','C4','C6',...
             'CP5','CP3','CP1','CPZ','CP2','CP4','CP6','M2',...
             'P7','P5','P3','P1','PZ','P2','P4','P6','P8',...
             'PO7','PO5','PO3','POZ','PO4','PO6','PO8',...
             'O1','OZ','O2'};
 	% eName   = {'C2','C4','C6','T8',...
		%        'CP2','CP4','CP6','TP8','M2',...
	 %           'P2','P4','P6','P8',...
	 %           'PO4','PO6','PO8',...
	 %           'O2','CB2'};

elseif strcmpi('m1m2',ldrType)
	
	eType = 'haveM1';
	eName = {'C5','C3','C1','CZ','C2','C4','C6','M1',...
             'CP5','CP3','CP1','CPZ','CP2','CP4','CP6','M2',...          
             'P7','P5','P3','P1','PZ','P2','P4','P6','P8',...          
             'PO7','PO5','PO3','POZ','PO4','PO6','PO8',...          
             'O1','OZ','O2'};          
                        

elseif strcmpi('m1m2Lps',ldrType)
	
	eType = 'haveM1';
    eName = {'CP5','CP3','CP1','CPZ','CP2','CP4','CP6','TP8','M2',...          
             'P7','P5','P3','P1','PZ','P2','P4','P6','P8',...          
             'PO7','PO5','PO3','POZ','PO4','PO6','PO8',...          
             'O1','OZ','O2'};  
end

[eNo]       = getElectrodeNo_bcl(eType,eName);
snrskipNum  = 1;


% pptTitleName = ['_power_',ldrType];
% [opt,ppt] = figPPT(pptTitleName);

allFolder = dir;

for iSub = 1:nSub

    isubName=allFolder(iSub+2).name;

    cd(fullfile(fullfile(pwd,isubName),'DATA'));

		for icnt = 1:size(cntfilename,2)

				icntfilename         =  cntfilename{icnt};

			    [freBins,y,trialSNR] = cntfft_bcl(icntfilename,isEpoched,moreBinNum,stimTrigger,freRange,[],fftType,ldrType);

				% %%%%%%%%%%%%%%%%%%%%%   compute SNR  %%%%%%%%%%%%%%%%%
				if ~strcmpi('trialFFT',fftType)
					for ie  = 1:size(y,2)       
						for iFrebin = snrskipNum+2:size(y,1)-(snrskipNum+1)
							SNR(iFrebin-(snrskipNum+1),ie) = y(iFrebin,ie)/mean(y([iFrebin-(snrskipNum+1),iFrebin+(snrskipNum+1)],ie));
						end
					end	
				else
					SNR = squeeze(mean(trialSNR,3)); 
				end
				% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				allSNR(:,:,icnt,iSub)     = SNR(:,eNo);
				allFreBins(:,:,icnt,iSub) = freBins;
				ally(:,:,icnt,iSub)       = y(:,eNo);
		end

    cd('..');
    cd('..');    
end



%%%%%%%%%%    find closest 2 frePoint  %%%%%%%

for iFre =1:size(fre,2)

		[minValue minPosition]  = min(abs(freBins-fre(iFre)));
		freBins(minPosition)    = 0;
		[minValue minPosition2] = min(abs(freBins-fre(iFre)));
		freTimePoint(1,[2*iFre-1,2*iFre]) = [minPosition, minPosition2];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% freBinStepNum = 15;
% freBinStep    = diff(freRange)/freBinStepNum;


colorbarRange = [0,0.2];


[myColormap]  = makeColormap(30,0,1,1);	

% ytick = round([0:size(freBins,2)/freBinStepNum:size(freBins,2)]);
% ytick(13) =  ytick(13)+(120/7-17)*size(freBins,2)/freBinStepNum;
% ytickLabel = [freRange(1):freBinStep:freRange(2)];

% ytickLabel(13) = 17.1;

doubleFrePower = ally(freTimePoint,:,:,:);
doubleFreSNR   = allSNR(freTimePoint-snrskipNum-1,:,:,:);



%%%%%%%%%%%%%%%%%% bug  %%%%%%%%%%%%%%%%%
if nSub>1
	for iFre = 1:size(fre,2)
		  frePower(iFre,:,:,:) = mean(doubleFrePower([2*iFre-1,2*iFre],:,:,:),1) ;
		  freSNR(iFre,:,:,:) = mean(doubleFreSNR([2*iFre-1,2*iFre],:,:,:),1) ;
	end
	frePower(:,:,:,size(ally,4)+1) = mean(frePower,4);
	freSNR(:,:,:,size(ally,4)+1) = mean(freSNR,4);	
else
	for iFre = 1:size(fre,2)
		  frePower(iFre,:,:) = mean(doubleFrePower([2*iFre-1,2*iFre],:,:),1) ;
		  freSNR(iFre,:,:) = mean(doubleFreSNR([2*iFre-1,2*iFre],:,:),1) ;
	end
	frePower(:,:,:,2) = frePower;
	freSNR(:,:,:,2) = freSNR;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subfigWidth  = 200/1920;
% subFigHeight = 150/1080;
% horStart     = 10/1920;
% verStart     = 10/1080;

H=figure;
set(H, 'position', get(0,'ScreenSize'));

nSub=nSub+1;

for iSub = 1 :nSub

		for iFre = 1 :size(fre,2)

			cFre = fre(iFre);


			if iSub<nSub
				cSubcFrePower = squeeze(frePower(iFre,:,:,iSub));
				cSubcFrePower = sort(cSubcFrePower,'descend');
				newFrePower(iFre,:,:,iSub) = cSubcFrePower;
			else 
				newFrePower(:,:,:,iSub)  = mean(newFrePower,4);
				cSubcFrePower            = squeeze(newFrePower(iFre,:,:,iSub));

			end
			hh=subplot(size(fre,2),nSub,iSub+nSub*(iFre-1));
			imagesc(cSubcFrePower);
			box off;
			caxis([colorbarRange]);
			% subFigRect = [horStart+1.1*(iSub-1)*subfigWidth verStart+1.1*(iFre-1)*subFigHeight subfigWidth subFigHeight];
% 			set(hh,'pos',subFigRect);
			% set(gca,'ytick',[1:1:size(eName,2)]);
			% set(gca,'ytickLabel',eName);
			% set(gca,'xTickLabel', cntfilename);
			view(-90,-90);
         	colormap(myColormap);

			if iSub<nSub
					title([num2str(iSub),'-power-',ldrType,'-',num2str(cFre),'Hz']);

			elseif iSub==nSub

					title(['maxMeanPower-',ldrType,'-',num2str(cFre),'Hz']);
            end
		end
end


hBar = colorbar;
rect = get(hBar,'Position');
set(hBar,'Position',[rect(1)+0.05 rect([2 3 4])]);
cfigFileName = [ldrType,'_',fftType,'_power.bmp'];
saveas(gcf,cfigFileName);

% ppt=ppt.addImageSlide([subNum,'-power-',ldrType],fullfile(pwd,cfigFileName),opt);
% close all;

% ppt.saveAs(fullfile(pwd,pptTitleName));

% ppt.close;

SNRcolorbarRange = [0,4];


H=figure;
set(H, 'position', get(0,'ScreenSize'));

for iSub = 1:nSub

		for iFre = 1:size(fre,2)

			cFre = fre(iFre);

				if iSub<nSub
					cSubcFreSNR = squeeze(freSNR(iFre,:,:,iSub));
					cSubcFreSNR = sort(cSubcFreSNR,'descend');
					newFreSNR(iFre,:,:,iSub) = cSubcFreSNR;
				else 
					newFreSNR(:,:,:,iSub)  = mean(newFreSNR,4);
					cSubcFreSNR            = squeeze(newFreSNR(iFre,:,:,iSub));
				end


			subplot(size(fre,2),nSub,iSub+nSub*(iFre-1));
			imagesc(cSubcFreSNR);
			box off;
			caxis([SNRcolorbarRange]);
			% set(gca,'ytick',[1:1:size(eName,2)]);
			% set(gca,'ytickLabel',eName);
			% set(gca,'xTickLabel', cntfilename);
			view(-90,-90);
            colormap(myColormap);
			if iSub<nSub
					title([num2str(iSub),'-SNR-',ldrType,'-',num2str(cFre),'Hz']);

            elseif iSub==nSub
					title(['maxMeanSNR-',ldrType,'-',num2str(cFre),'Hz']);
			end

		end
end

% tightfig;
hBar = colorbar;
rect = get(hBar,'Position');
set(hBar,'Position',[rect(1)+0.05 rect([2 3 4])]);
cfigFileName = [ldrType,'_',fftType,'_SNR.bmp'];
saveas(gcf,cfigFileName);



% close all;
% ppt.saveAs(fullfile(pwd,pptTitleName));

% ppt.close;



















% for ifig=1:figNum

% 	H=figure;
% 	set(H, 'position', get(0,'ScreenSize'));

% 	for figichan=1:eNumPerFig;

% 		ichan=eNumPerFig*(ifig-1)+figichan;

% 		subplot(1,eNumPerFig,figichan);

% 			for icnt = 1:size(cntfilename,2)
% 			    line([8 8],snrlineRange(:,:,icnt),'linestyle','--','linewidth',linewidth,'color',[lineColor transparency]);
% 				hold on;
% 			    line([120/7 120/7],snrlineRange(:,:,icnt),'linestyle','--','linewidth',linewidth,'color',[lineColor transparency]);
% 				line([12 12],snrlineRange(:,:,icnt),'linestyle','--','linewidth',linewidth,'color',[lineColor transparency]);
% 				line([15 15],snrlineRange(:,:,icnt),'linestyle','--','linewidth',linewidth,'color',[lineColor transparency]);
% 				plot(allFreBins(:,[snrskipNum+2:end-(snrskipNum+1)],icnt),allSNR(:,ichan,icnt),'color',figColor);
%                 text(mean(freRange),snrlineRange(:,end,icnt)+ytickStep,cntfilename{icnt});
% 			end
			
% 			box off;
% 			xlabel('Frequency Hz');
% 			ylabel('|y(f)|');
% 			title(eName{ichan});
% 			xlim(freRange);
% 		    ylim(snryRange);
% 		    set(gca,'xtick',xtickRange);            
% 		    set(gca,'ytick',ytickRange);
% 	end

% 	if ~exist('ldrType','var')
% 		saveas(gcf,['snr_merge_fig',num2str(ifig),'_skip',num2str(snrskipNum),'.bmp']);
% 	else
% 		saveas(gcf,['snr_merge_fig',num2str(ifig),'_',ldrType,'_skip',num2str(snrskipNum),'.bmp']);
% 	end 

% end

% close all;


% 
% 
% cntName = 'freTest-counterPhaseHalf.cnt';
% 
% [freBins,y]  = cntfft_bcl(cntName,[5 40]);
% plot(freBins,y);
% box off;
% linewidth    = 2;
% lineRange    = [0,0.1];
% % lineRange    = [0.2,0.3];
% lineColor    = [0 0 0];
% transparency = 0.5;
% 
% line([8 8],lineRange,'linestyle','--','linewidth',linewidth,'color',[lineColor transparency]);
% hold on;
% line([120/7 120/7],lineRange,'linestyle','--','linewidth',linewidth,'color',[lineColor transparency]);
% line([12 12],lineRange,'linestyle','--','linewidth',linewidth,'color',[lineColor transparency]);
% line([15 15],lineRange,'linestyle','--','linewidth',linewidth,'color',[lineColor transparency]);
% 
% % line([8 8]*2,lineRange,'linestyle','--','linewidth',linewidth,'color',[lineColor transparency]);
% % line([120/7 120/7]*2,lineRange,'linestyle','--','linewidth',linewidth,'color',[lineColor transparency]);
% % line([12 12]*2,lineRange,'linestyle','--','linewidth',linewidth,'color',[lineColor transparency]);
% % line([15 15]*2,lineRange,'linestyle','--','linewidth',linewidth,'color',[lineColor transparency]);
% 
% typeName = cntName([9:end-4]);
% title(typeName);
% saveas(gcf,[typeName,'.bmp']);
