
function [data ] = geteachSubAvgData_bcl(directory,condition,subNum,eType,eName,timeBin,baseCor)


% written by Li Aisu 

ampdata  = [];
eavgAmp  = [];
eNo   = getElectrodeNo_bcl(eType,eName);


for iCondition=1:size(condition,2)
	
	iFile = 0;
	titleString(1,size(eNo,2)*(iCondition-1)+1) = condition(iCondition);
	titleString(2,[size(eNo,2)*(iCondition-1)+1:size(eNo,2)*iCondition]) = eName;

    for iAvg = subNum

		iFile      = iFile + 1;
		eachAvgDir = fullfile(fullfile(directory,condition{iCondition}),[num2str(iAvg),'-',condition{iCondition},'.avg']);
        [signal,chan_names,variance, pnts, rate, xmin, xmax,nsweeps] = loadavg_bcl(eachAvgDir);
        baseNo   = numel(xmin*rate:0); 
        point    = baseNo+timeBin;
        if baseCor(1)==1 % do baseline  
                baseCorPoint = baseNo+baseCor(2:3);
                baseMinus = mean(signal([baseCorPoint(1):baseCorPoint(end)],:),1);
                
                for ichan = 1:size(signal,2)
                        signal(:,ichan) = signal(:,ichan) - baseMinus(1,ichan);
                end

        elseif  baseCor(1)==0 % notdo

        end
        eachConditionData(iFile,:) = mean(signal([point(1):point(end)],eNo),1);
    end

    ampdata = [ampdata eachConditionData];

    %%%%%  compute avg of needed electrodes 2018,01,25 %%%%
    eachConditionAvgAmp = mean(eachConditionData,2);
    eavgAmp = [eavgAmp  eachConditionAvgAmp];
end


data.Tamplitude  = ampdata;

data.eavgAmplitude = eavgAmp;

data.titleString = titleString;