
function plotNtopo_bcl(avgFile,topoNum,timeBin,amplitudeRange,viewPoint,figDir)

% written by Li Aisu
% currentPath should be avgFiles' path  

H = figure;
set(H, 'position', get(0,'ScreenSize'));

% avgFile={'T_easy_IOR_2days_ex1&9&12.avg','T_mid_IOR_2days_ex1&9&12.avg','T_hard_IOR_2days_ex1&9&12.avg','T_IOR_2days_ex1&9&12.avg'};

% avgFile={'T_hor_IOR_S16.avg','T_ver_IOR_S16.avg','T_IOR_S16.avg'};

% topoNum=10;



timeStart = timeBin(1);
avgFileNum = size(avgFile,2);

for i=1:topoNum

		subplot(avgFileNum,topoNum,i);
		[handle,Zi,contourLineVals] = topoplot_bcl(avgFile{1},'EEGChans',{},'baselinePeriod',200,'plotAreas',timeBin,'maplimitsDouble',amplitudeRange,'style','both','viewpoint',viewPoint);
		hold on;

		subplot(avgFileNum,topoNum,i+topoNum);
		[handle,Zi,contourLineVals] = topoplot_bcl(avgFile{2},'EEGChans',{},'baselinePeriod',200,'plotAreas',timeBin,'maplimitsDouble',amplitudeRange,'style','both','viewpoint',viewPoint);
		hold on;

		subplot(avgFileNum,topoNum,i+2*topoNum);
		[handle,Zi,contourLineVals] = topoplot_bcl(avgFile{3},'EEGChans',{},'baselinePeriod',200,'plotAreas',timeBin,'maplimitsDouble',amplitudeRange,'style','both','viewpoint',viewPoint);
		hold on;

		% subplot(4,topoNum,i+3*topoNum);
		% [handle,Zi,contourLineVals] = topoplot_bcl(avgFile{4},'EEGChans',{},'baselinePeriod',200,'plotAreas',timeBin,'maplimitsDouble',amplitudeRange,'style','both','viewpoint','top');
		% hold on;



        if  topoNum>1       
	        if timeBin(end)<200
				timeBin=timeBin+10;
			else
				timeBin=timeBin+20;
			end
		end

end


set(gcf,'color','w');

if topoNum>1; 
      print('-dpdf','-painters',fullfile(figDir,[num2str(timeStart),'ms_',num2str(topoNum),'topo.pdf']));
      % print('-dpdf','-painters',[figureName,'_fig4framePercent.pdf']);
else
	  print('-dpdf','-painters',fullfile(figDir,[num2str(timeBin(1)),'To',num2str(timeBin(end)),'topo.pdf']));
end

% close all;






% %%%%%%%%%%  get erpData  %%%%%%%%%%%    
% clear all;
% totsub   = [1:16];
% filepath = 'D:\expData\HorVer_IOR\Exp1_20161107\Data\rawData_shift\average';

% analysisParameter.polarDirection = 0;
% isJackknife                 = false;
% electrodeNames              = {'F1','FZ','F2','FC1','FCZ','FC2','C1','CZ','C2'}; % 'P7','P5','P3','P1','PZ','P2','P4','P6','P8','PO7','PO5','PO3','POZ','PO4','PO6','PO8'
% analysisParameter.TimeRange = [240,280];
% analysisParameter.PeakType  = 'average';

% analysisParameter.PeakCritertion  = 0.5;
% analysisParameter.avgFile       = {'hor_cue_S16','hor_uncue_S16','ver_cue_S16','ver_uncue_S16'};
% avgNameFormat='sub%';

% MergeMode = [];
% [data]    = erpDetection(totsub,filepath,isJackknife,electrodeNames,analysisParameter,MergeMode,[],avgNameFormat);



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%% plotErp %%%%%%%%%%%%

% files={'T_hor_cue_S16.avg','T_hor_uncue_S16.avg','T_ver_cue_S16.avg','T_ver_uncue_S16.avg','T_hor_IOR_S16.avg','T_ver_IOR_S16.avg'}; % ,'T_hor_IOR_S16.avg','T_ver_IOR_S16.avg'

% component=1; % 1 for P1, 2 for N1, 3 for Nd 

% switch component
% 	case 1
% 		chanNames={'P7','P5','P6','P8','PO7','PO5','PO6','PO8'};
% 		YRange=[-3,4];
% 		Yticks=[-1:1:1];
% 	case 2
% 		chanNames={'PO7','P7','P5','P3','P1','PZ','P2','P4','P6'};
% 		YRange=[-3,4];
% 		Yticks=[-1:1:1];
% 	case 3
% 		chanNames={'C1','CZ','C2','CP1','CPZ','CP2'};
% 		YRange=[-3,4];
% 		Yticks=[-1:1:1];			
% end

% xTimeBinMsec =[-200:2:400]; 
% Xticks=[-200:100:400];
% axisLineWidth=1;
% colorSet=[1 0 0 ; 0 0 1; 0 1 0; 0 0 0 ; 0.3 0.3 0.3 ; 0.8 0.8 0.8];
% [dataf,tData]=plotErp_bcl(files,chanNames,xTimeBinMsec,YRange,Xticks,Yticks,axisLineWidth,colorSet); 
% set(gcf,'color','w');

% switch component
% 	case 1
% 		print('-depsc','-painters','D:\expData\HorVer_IOR\Exp1_20161107\Figure\RawFigure\erp&IOR_P1.eps');
% 	case 2
% 		print('-depsc','-painters','D:\expData\HorVer_IOR\Exp1_20161107\Figure\RawFigure\erp&IOR_N1.eps');
% 	case 3
% 		print('-depsc','-painters','D:\expData\HorVer_IOR\Exp1_20161107\Figure\RawFigure\erp&IOR_Nd.eps');	

% end



		
% %%%%%%%%%%%%ssvep_2loc%%%%%%%%%%%%%%%%%%%%
% [handle,Zi,contourLineVals] = topoplot_bcl(eight,'EEGChans',e,'style','both','viewpoint','top');




% for i=1:size(chans,1)
%    e{1,i}=deblank(char(chans(i,:)));
% end

	
% subplot(3,1,1)
% 		[handle,Zi,contourLineVals] = topoplot_bcl(eight,'EEGChans',e,'maplimitsDouble',[-0.2,0.2],'style','both','viewpoint','top');
% 		% title('8hz');

% subplot(3,1,2)
% 		[handle,Zi,contourLineVals] = topoplot_bcl(sixteen,'EEGChans',e,'maplimitsDouble',[-0.2,0.2],'style','both','viewpoint','top');
% 		% title('16hz');
% subplot(3,1,3)
% 		[handle,Zi,contourLineVals] = topoplot_bcl(twenty,'EEGChans',e,'maplimitsDouble',[-0.2,0.2],'style','both','viewpoint','top');
% 		% title('20hz');



% expDesign.wsDesign ={'level','cue','e';3,2,size(P1eNo,2)};
% expDesign.InterRSE   = {'level*cue','cue*level'};
% cd('singleTimeSingleE');


% % expDesign.wsDesign ={'level','cue';3,2};
% % expDesign.InterRSE   = {'level*cue','cue*level'};
% % % [result] = runANOVA(fullfile(pwd,'rtacc'),rtacc,expDesign,1,'rt/acc','C:\Program Files (x86)\IBM\SPSS\Statistics\21');

% [result] = runANOVA(fullfile(pwd,[num2str(P1timeWin(1)),'To',num2str(P1timeWin(end)),'_',P1eName{:}]),meanP1,expDesign,1,[num2str(P1timeWin(1)),'To',num2str(P1timeWin(end)),'_',P1eName{:}],'C:\Program Files (x86)\IBM\SPSS\Statistics\21');


% spssRT=[];
% spssACC=[];
% for i=1:16;
%   spssRT= [spssRT;validRt([4*i-3:4*i],4)'];
%   spssACC=[spssACC;ACC([4*i-3:4*i],5)'];
% end


% horUncue_acc
% horCue_acc
% VerUncue_acc
% VerCue_acc