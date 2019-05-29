function analyBeh_SSVEP_8locDet
% written by Aisu Li
% Soochow University
% 2017-02-28

  clear all

  cd '/Users/las/Documents/MATLAB/preData_det_pre2';
  accRange=[0.75,1];
  accStick=[accRange(1):0.05:accRange(2)];  
  RtRange=[200,350; 100,250];
  RtStick=[RtRange(1,1):30:RtRange(1,2);RtRange(2,1):30:RtRange(2,2);];
  allMatFiles = dir(fullfile(pwd,'*.mat'));
  matNames = struct2cell(allMatFiles)';
  allSubData=[];
  allSubDataInfo={'subNum';'blockNo';'cueType/cue1';'cue2';'targetLoc';'hemi';'trigNum';'acc';'rt';'conditions';'UpperLimit of RT: mean+2.5SD';'LowerLimir of RT: 100'};
  eachConTrial=[4 4 4 4 4 4 4 4]'; % [8 8 4 12] [8 8 4 4 8]

  for iMatFile = 1:size(matNames,1)
      iMatData = [];
      iMatName = [];
      
      iMatName = matNames{iMatFile,1};
      
      iMat     = load(iMatName);
      
      iMatData(:,[1,2]) = repmat([str2num(iMat.SSVEP_IOR_8locDet.session),str2num(iMat.SSVEP_IOR_8locDet.num)],[size(iMat.SSVEP_IOR_8locDet.designMatrix,1),1]);
      iMatData(:,[3:9]) = [iMat.SSVEP_IOR_8locDet.designMatrix(:,[1,6,2:4]),...
                           iMat.SSVEP_IOR_8locDet.acc(:),...
                           iMat.SSVEP_IOR_8locDet.responseTimes(:)];

      % for iTrial=1:size(iMat.SSVEP_IOR_8locDet.designMatrix,1)

         % if iMatData(iTrial,5)<5 % isTarget
         
             % if mod(iMatData(iTrial,3),4)<2 % cues(2,1;3,4) separate 
 
             %     if ismember(iMatData(iTrial,5),iMatData(iTrial,[3,4]))
             %          iMatData(iTrial,10)= 2; % cue separate
             %     elseif iMatData(iTrial,5)==mean(iMatData(iTrial,[3,4]))
             %          iMatData(iTrial,10)= 3; % uncue intermediate 
             %     else 
             %          iMatData(iTrial,10)= 4; % uncue opposite
             %     end

             % else  % cues(1,3;2,4) adjacent
             %     if ismember(iMatData(iTrial,5),iMatData(iTrial,[3,4]))
             %          iMatData(iTrial,10)= 1; % cue adjacent
             %     else
             %          iMatData(iTrial,10)= 5; % uncue  oppsite    
             %     end
             % end

         % end % isTarget 

      % end % iTrial


       for iTrial = 1:size(iMat.SSVEP_IOR_8locDet.designMatrix,1)

              if  ismember(iMatData(iTrial,5),[2,3]) % targetLoc=2or3
     
                      if iMatData(iTrial,5) == iMatData(iTrial,3)         
                               iMatData(iTrial,10)=1; % adjacent cue2,3 
         
                      elseif iMatData(iTrial,5) == 5-iMatData(iTrial,3)                                    
                               iMatData(iTrial,10)=2; % adjacent uncue2,3 
         
                      elseif iMatData(iTrial,5) == iMatData(iTrial,4)         
                               iMatData(iTrial,10)=3; % separate cue2,3

                      elseif iMatData(iTrial,5) == mean(iMatData(iTrial,[3,4]))
                               iMatData(iTrial,10)=4; % separate Inter2,3 % -----makeSure     
         
                      end

              elseif ismember(iMatData(iTrial,5),[1,4]) 


                      if iMatData(iTrial,5) == iMatData(iTrial,4)         
                               iMatData(iTrial,10)=5; % adjacent cue1,4 
         
                      elseif iMatData(iTrial,5) == 5-iMatData(iTrial,4)                                    
                               iMatData(iTrial,10)=6; % adjacent uncue1,4 
         
                      elseif iMatData(iTrial,5) == iMatData(iTrial,3)         
                               iMatData(iTrial,10)=7; % separate cue1,4

                      else%if iMatData(iTrial,5) == 5-iMatData(iTrial,3) 
                               iMatData(iTrial,10)=8; % separate Inter1,4 % -----makeSure     
         
                      end
              end
      end     
      allSubData = [allSubData;iMatData];
	end% iMat 



  rowFilter = ~~allSubData(:,8) & allSubData(:,10)>0;  % correct nonCatch


   % cutoff Values
   [gm,gsd,gname,corrgcout] = grpstats(allSubData(rowFilter,9),{allSubData(rowFilter,1),allSubData(rowFilter,10)},{'mean','std','gname','numel'});
   gname                    = cellfun(@str2double,gname);

   [zscores]                     = shiftzs_BCL([corrgcout],0);   


   for iRow = 1:size(gname,1)
       
       conFilter = allSubData(:,1) == gname(iRow,1)  & allSubData(:,10) == gname(iRow,2) ;
       
       allSubData(conFilter,11) = gm(iRow) + zscores(iRow)*gsd(iRow); 
       allSubData(conFilter,12) = 100;
   end
   
   ACC = [gname,corrgcout];

   ACC(:,end+1) = ACC(:,end)./repmat(eachConTrial*numel(unique(allSubData(:,2))),numel(unique(allSubData(:,1))),1); % isTarget


  % -----------allSubMeanACC----------------/
    [allSubgmCount,allSubgsd,allSubgname,allSubcorrgcout]=grpstats(ACC(:,3),{ACC(:,2)},{'mean','std','gname','numel'});
    [allSubgm,allSubgsem,allSubgname,allSubcorrgcout]=grpstats(ACC(:,4),{ACC(:,2)},{'mean','sem','gname','numel'});
    allSubgname = cellfun(@str2double,allSubgname);
    ACC([size(ACC,1)+1:size(ACC,1)+numel(unique(gname(:,2)))],:)=[allSubgsem,allSubgname,allSubgmCount,allSubgm]; 
  % ----------------------------------------\



    % caculate validRt
    rowFilter             = ~~allSubData(:,8) & allSubData(:,10)>0 & allSubData(:,9)>allSubData(:,12) & allSubData(:,9)<=allSubData(:,11);
    
    [gm,gsem,gname,gcout] = grpstats(allSubData(rowFilter,9),{allSubData(rowFilter,1),allSubData(rowFilter,10)},{'mean','sem','gname','numel'});
    
    gname                 = cellfun(@str2double,gname);
    
    validRt               = [gname,gm,gsem,gcout];
    
    OutlierNum            = [gname, corrgcout - gcout]; 

   %-----------allSubMeanRT----------------/
    [allSubgmCount,allSubgsem,allSubgname,allSubcorrgcout] = grpstats(validRt(:,end),{validRt(:,2)},{'mean','sem','gname','numel'});
    [allSubgm,allSubgsem,allSubgname,allSubcorrgcout]      = grpstats(validRt(:,end-2),{validRt(:,2)},{'mean','sem','gname','numel'});
    allSubgname = cellfun(@str2double,allSubgname);
    validRt([size(validRt,1)+1:size(validRt,1)+numel(unique(gname(:,2)))],:) = [allSubgsem,allSubgname,allSubgm,allSubgsem,allSubgmCount]; 
   %---------------------------------------\

    nPlot = size(validRt,1)/size(allSubgname,1);

%-----------different displacement of nplot---------------/
      switch nPlot
           case 1 
            nRow=1;nColm=1;
           case 2
            nRow=1;nColm=2;
           case 3
            nRow=1;nColm=3;
           case 4
            nRow=2;nColm=2;     
           case 5
            nRow=2;nColm=3;       
           case 6
            nRow=2;nColm=3;  
           case 7
            nRow=3;nColm=3;       
           case 8
            nRow=3;nColm=3;       
           case 9
            nRow=3;nColm=3;      
           case 10
            nRow=3;nColm=4;       
           case 11
            nRow=3;nColm=4;       
           case 12
            nRow=3;nColm=4;      
           case 13
            nRow=4;nColm=4;       
           case 14
            nRow=4;nColm=4;       
           case 15
            nRow=4;nColm=4;      
           case 16
            nRow=4;nColm=4;   
           case 17
            nRow=4;nColm=5;   
           case 18
            nRow=4;nColm=5;               
           case 19
            nRow=4;nColm=5;   
           case 20
            nRow=4;nColm=5;   
           case 21
            nRow=5;nColm=5;   
           case 22
            nRow=5;nColm=5;   
           case 23
            nRow=5;nColm=5;      
           case 24
            nRow=5;nColm=5;                            
        end
%----------------------------------------------------------\
    
    barColor = [200 200 169; 252 157 154; 131 175 155;254 67 101 ];

    validRtLoc23 = validRt(validRt(:,2)<5,:);
    validRtLoc14 = validRt(validRt(:,2)>4,:);
    ACCLoc23     = ACC(ACC(:,2)<5,:);
    ACCLoc14     = ACC(ACC(:,2)>4,:);

%----------------plot all conditions in one figure---------------/
    figure(1);
    set (gcf,'Position',[50,50,1800,900], 'color','w');
    for iSub=1:nPlot
         
         iSubValidRtLoc23 = validRtLoc23([(iSub-1)*4+1:iSub*4],3);
         ierror23         = validRtLoc23([(iSub-1)*4+1:iSub*4],4);

         iSubValidRtLoc14 = validRtLoc14([(iSub-1)*4+1:iSub*4],3);
         ierror14         = validRtLoc14([(iSub-1)*4+1:iSub*4],4);

         subplot(nRow,nColm,iSub);         
         h  = bar([iSubValidRtLoc23';iSubValidRtLoc14'],'grouped');
         ch = get(h,'children');
         set(ch{1},'Facecolor',barColor(1,:));
         set(ch{2},'Facecolor',barColor(2,:));
         set(ch{3},'Facecolor',barColor(3,:));
         set(ch{4},'Facecolor',barColor(4,:));

         set(gca,'xTickLabel',{'Target23','Target14'});
           
         if iSub<nPlot
               title(['Sub',num2str(iSub),' RT']);
         else
              l = legend([ch{1} ch{2} ch{3} ch{4} ],'AjCue','AjUnc','SepCue','SepUnc');
              set(l, 'Position',[0.9 0.2 0.1 0.01]);;
              legend('boxoff');
              title('allSubMean');
         end
         
         if iSub ~=4
            ylim(RtRange(1,:));
            set(gca,'ytick',RtStick(1,:));
         else 
            ylim(RtRange(2,:));
            set(gca,'ytick',RtStick(2,:));
         end 
         hold on        
         errorbar([0.73:0.18:1.4], iSubValidRtLoc23, ierror23, 'k', 'linestyle', 'none', 'linewidth', 1);
         errorbar([0.73:0.18:1.4]+1, iSubValidRtLoc14, ierror14, 'k', 'linestyle', 'none', 'linewidth', 1);
         grid on ; box off ;  hold on

    end     
    saveas(gcf,['RT',num2str(nPlot-1),'_8con.bmp']);

%----------------------------plot target2,3 or 1,4-------------------------/

    figure(2);
    set (gcf,'Position',[50,50,1800,900], 'color','w');
    for iSub=1:nPlot

         
         iSubValidRtLoc23 = validRtLoc23([(iSub-1)*4+1:iSub*4],3);
         ierror23         = validRtLoc23([(iSub-1)*4+1:iSub*4],4);

         subplot(nRow,nColm,iSub);         
         h  = bar(iSubValidRtLoc23,0.5);
         ch = get(h,'children');
         set(ch,'FaceVertexCData',barColor)
         set(gca,'xTickLabel',{'AjCue','AjUnc','SepCue','SepUnc'});
           
         if iSub<nPlot
               title(['Sub',num2str(iSub),' RT']);
         else
              % l = legend([ch{1} ch{2} ch{3} ch{4} ],'cueAj','cueSep','uncueInter','uncueOps');
              % legend('boxoff');
              title('allSubMean');
         end

         if iSub ~=4
            ylim(RtRange(1,:));
            set(gca,'ytick',RtStick(1,:));
         else 
            ylim(RtRange(2,:));
            set(gca,'ytick',RtStick(2,:));
         end 

         hold on        
         errorbar([1:4], iSubValidRtLoc23, ierror23, 'k', 'linestyle', 'none', 'linewidth', 1);
         grid on ; box off ;  hold on
    end     
    saveas(gcf,['RT',num2str(nPlot-1),'_Loc23.bmp']);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(3);
    set (gcf,'Position',[50,50,1800,900], 'color','w');
    for iSub=1:nPlot

         iSubValidRtLoc14 = validRtLoc14([(iSub-1)*4+1:iSub*4],3);
         ierror14         = validRtLoc14([(iSub-1)*4+1:iSub*4],4);

         subplot(nRow,nColm,iSub);         
         h  = bar(iSubValidRtLoc14,0.5);
         ch = get(h,'children');
         set(ch,'FaceVertexCData',barColor)
         set(gca,'xTickLabel',{'AjCue','AjUnc','SepCue','SepUnc'});
           
         if iSub<nPlot
               title(['Sub',num2str(iSub),' RT']);
         else
              % l = legend([ch{1} ch{2} ch{3} ch{4} ],'cueAj','cueSep','uncueInter','uncueOps');
              % legend('boxoff');
              title('allSubMean');
         end
         
         if iSub ~=4
            ylim(RtRange(1,:));
            set(gca,'ytick',RtStick(1,:));
         else 
            ylim(RtRange(2,:));
            set(gca,'ytick',RtStick(2,:));
         end 

         hold on        
         errorbar([1:4], iSubValidRtLoc14, ierror14, 'k', 'linestyle', 'none', 'linewidth', 1);
         grid on ; box off ;  hold on
    end     
    saveas(gcf,['RT',num2str(nPlot-1),'_Loc14.bmp']);
%-------------------------------------------------------------------------------\


%%%%%%%%%%%%%%%%%%% ACC %%%%%%%%%%%%%%%
    figure(4);
    set (gcf,'Position',[50,50,1800,900], 'color','w');

    for iSub=1:nPlot
         
         iSubACCLoc23 = ACCLoc23([(iSub-1)*4+1:iSub*4],4);
         iSubACCLoc14 = ACCLoc14([(iSub-1)*4+1:iSub*4],4);

         subplot(nRow,nColm,iSub);         

         h  = bar([iSubACCLoc23';iSubACCLoc14'],'grouped');
         ch = get(h,'children');

         set(ch{1},'Facecolor',barColor(1,:));
         set(ch{2},'Facecolor',barColor(2,:));
         set(ch{3},'Facecolor',barColor(3,:));
         set(ch{4},'Facecolor',barColor(4,:));

         set(gca,'xTickLabel',{'Target23','Target14'});
           
         if iSub<nPlot
               title(['Sub',num2str(iSub),' ACC']);
         else
              l = legend([ch{1} ch{2} ch{3} ch{4} ],'AjCue','AjUnc','SepCue','SepUnc');
              set(l, 'Position',[0.9 0.2 0.1 0.01]);;
              legend('boxoff');
              title('allSubMean');
              hold on;
              errorbar([0.73:0.18:1.4], iSubACCLoc23, ACCLoc23(end-3:end,1), 'k', 'linestyle', 'none', 'linewidth', 1);
              errorbar([0.73:0.18:1.4]+1, iSubACCLoc14, ACCLoc14(end-3:end,1), 'k', 'linestyle', 'none', 'linewidth', 1);
         end
         ylim(accRange);
         set(gca,'ytick',accStick);     
         grid on ; box off ;  hold on
    end    

    saveas(gcf,['ACC',num2str(nPlot-1),'_8con.bmp']);

%---------------------------------------------------------
    figure(5);
    set (gcf,'Position',[50,50,1800,900], 'color','w');
    for iSub=1:nPlot

         iSubACCLoc23 = ACCLoc23([(iSub-1)*4+1:iSub*4],4);

         subplot(nRow,nColm,iSub);         
         h  = bar(iSubACCLoc23,0.5);
         ch = get(h,'children');
         set(ch,'FaceVertexCData',barColor)
         set(gca,'xTickLabel',{'AjCue','AjUnc','SepCue','SepUnc'});
           
         if iSub<nPlot
               title(['Sub',num2str(iSub),' ACC']);
         else
              % l = legend([ch{1} ch{2} ch{3} ch{4} ],'cueAj','cueSep','uncueInter','uncueOps');
              % legend('boxoff');
              title('allSubMean');
              hold on        
              errorbar([1:4], iSubACCLoc23,ACCLoc23(end-3:end,1) , 'k', 'linestyle', 'none', 'linewidth', 1);
         end

         ylim(accRange);
         set(gca,'ytick',accStick);
         grid on ; box off ;  hold on
    end     
    saveas(gcf,['ACC',num2str(nPlot-1),'_Loc23.bmp']);

%----------------------------------------------------------
    figure(6);
    set (gcf,'Position',[50,50,1800,900], 'color','w');
    for iSub=1:nPlot

         iSubACCLoc14 = ACCLoc14([(iSub-1)*4+1:iSub*4],4);

         subplot(nRow,nColm,iSub);         
         h  = bar(iSubACCLoc14,0.5);
         ch = get(h,'children');
         set(ch,'FaceVertexCData',barColor)
         set(gca,'xTickLabel',{'AjCue','AjUnc','SepCue','SepUnc'});
           
         if iSub<nPlot
               title(['Sub',num2str(iSub),' ACC']);
         else
              % l = legend([ch{1} ch{2} ch{3} ch{4} ],'cueAj','cueSep','uncueInter','uncueOps');
              % legend('boxoff');
              title('allSubMean');
              hold on        
              errorbar([1:4], iSubACCLoc14,ACCLoc14(end-3:end,1) , 'k', 'linestyle', 'none', 'linewidth', 1);
         end

         ylim(accRange);
         set(gca,'ytick',accStick);
         grid on ; box off ;  hold on
    end     
    saveas(gcf,['ACC',num2str(nPlot-1),'_Loc14.bmp']);

    cd('result');
    save([num2str(nPlot-1),'.mat'],'allSubData','ACC','validRt');

    
% %----------for spss---------------


nCondition=numel(unique(validRt(:,2)));
spssRT=[];
for isub=1:size(validRt,1)/nCondition-1
 isubRT=validRt([isub*nCondition-7:isub*nCondition-4],3);
 spssRT=[spssRT;isubRT'];
end

spssACC=[];
for isub=1:size(ACC,1)/nCondition-1
 isubACC=ACC([isub*nCondition-7:isub*nCondition-4],4);
 spssACC=[spssACC;isubACC'];
end


