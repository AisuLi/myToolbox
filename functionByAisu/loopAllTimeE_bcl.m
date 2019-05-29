
function loopAllTimeE_bcl(condition,subNum,matOutputString,wsDesign,pValue)


% written by Li Aisu 

% %------------------prepare all data-----------------------------
% condition={'easy_diff_2days','easy_same_2days','mid_diff_2days','mid_same_2days','hard_diff_2days','hard_same_2days'};
% condition={'easy_IOR_2days','mid_IOR_2days','hard_IOR_2days'};
% wsDesign ={'level','cue';3,2};

% condition ={'hor_cue_S16','hor_uncue_S16','ver_cue_S16','ver_uncue_S16'};
% wsDesign ={'displace','cue';2,2};
% subNum=1:16;
% matOutputString='APP_m1m2_adp';



if ~exist('matOutputString','var')||isempty(matOutputString)
    matOutputString = '';
end 

if ~exist('pValue','var')||isempty(pValue)
    pValue = 0.05;
end 


for iCondition=1:size(condition,2)
    
    cd(condition{iCondition});


             iFile=0;

			 for iSubFile = subNum
			 	 signal=[];
				 eachName = [num2str(iSubFile),'-',condition{iCondition},'.avg'];
				 iFile=iFile+1;
				 [signal,chan_names,variance, pnts, rate, xmin, xmax,nsweeps]=loadavg_bcl(eachName);
				 signal=signal(:,[1:66]);
				 allsignal([(iFile-1)*size(signal,1)+1:iFile*size(signal,1)],[(iCondition-1)*size(signal,2)+1:iCondition*size(signal,2)])=signal;
		     end

    cd('..');

end

%-----------get singleTimeE data------------
startime=xmin*rate; % ms
startPoint=numel([xmin*rate:startime]);

% endPoint=startPoint+39;

expDesign.wsDesign = wsDesign;
% expDesign.InterRSE   = {'cue*level'};
expDesign.runSE   = 0;
isPrintInfo=0;
SPSSfolder='C:\Program Files (x86)\IBM\SPSS\Statistics\21';



for iTimebinPoint=startPoint:size(signal,1)

	for ie=1:size(signal,2)

        singletimeEdata=[];
        result=[];
        singletimeEdata=allsignal([iTimebinPoint:size(signal,1):size(allsignal,1)],[ie:size(signal,2):size(allsignal,2)]);
		
		iTimebin=iTimebinPoint-201;
		outputFileString=[num2str(iTimebin),'_',num2str(ie),'e'];

	    % [result] = runANOVA(fullfile(pwd,outputFileString),singletimeEdata,expDesign,isPrintInfo,outputFileString,SPSSfolder);
	    [result] = ranova_bcl(singletimeEdata,expDesign);

%--------------------main effect of 3level------------------
        % if result.MauchlysTestofSphericitya{3,5}<0.05
        %    levelRsult(iTimebinPoint,ie)=str2num(result.TestsofWithinSubjectsEffects{3,7});
        % else
        %    levelRsult(iTimebinPoint,ie)=str2num(result.TestsofWithinSubjectsEffects{2,7});
        % end

        levelRsult(iTimebinPoint,ie)=result.withinSubjectEffect{2,6};

        if levelRsult(iTimebinPoint,ie)>pValue
        	levelSig(iTimebinPoint,ie)=NaN;
        else
        	levelSig(iTimebinPoint,ie)=levelRsult(iTimebinPoint,ie);
        end

%--------------------main effect of cueing------------------
        % cueRsult(iTimebinPoint,ie)=str2num(result.TestsofWithinSubjectsEffects{10,7});
          cueRsult(iTimebinPoint,ie)=result.withinSubjectEffect{4,6};
  
        if cueRsult(iTimebinPoint,ie)>pValue
        	cueSig(iTimebinPoint,ie)=NaN;
        else
        	cueSig(iTimebinPoint,ie)=cueRsult(iTimebinPoint,ie);
        end

%--------------------interaction of cueing*level------------------
        % levelcueRsult(iTimebinPoint,ie)=str2num(result.TestsofWithinSubjectsEffects{18,7});
          levelcueRsult(iTimebinPoint,ie)=result.withinSubjectEffect{6,6};  

        if levelcueRsult(iTimebinPoint,ie)>pValue
           levelcueSig(iTimebinPoint,ie)=NaN;
        else
           levelcueSig(iTimebinPoint,ie)=levelcueRsult(iTimebinPoint,ie);
        end
        
        spss.cueRsult      = cueRsult;
        spss.cueSig        = cueSig;
        spss.levelRsult    = levelRsult;       
        spss.levelSig      = levelSig; 
        spss.levelcueRsult = levelcueRsult;
        spss.levelcueSig   = levelcueSig;        
		save(fullfile(pwd,[num2str(startime),matOutputString,'.mat']),'spss');

    end

end




