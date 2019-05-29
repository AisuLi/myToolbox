
function [data] = geteachSubContinuousData_bcl(directory,condition,subNum,eType,eName,timeBin)


% written by Li Aisu 


% directory = 'D:\expData\difficultyIOR_20160415\Data\scanData\20170731\average';
% condition = {'easy_diff_2days','easy_same_2days','mid_diff_2days','mid_same_2days','hard_diff_2days','hard_same_2days'};
% subNum    = [2:8,10,11,13:18];
% eType ='haveM1';
% eName = {'PO7','PO5','PO6','PO8','P1','PZ','P2'};
% timeBin = [100,500];

ampdata  = [];

eNo   = getElectrodeNo_bcl(eType,eName);


iFile = 0;

for iAvg = subNum

        iFile = iFile + 1;

        for iCondition=1:size(condition,2)
            
            titleString(1,size(eNo,2)*(iCondition-1)+1) = condition(iCondition);
            titleString(2,[size(eNo,2)*(iCondition-1)+1:size(eNo,2)*iCondition]) = eName;

            eachAvgDir = fullfile(fullfile(directory,condition{iCondition}),[num2str(iAvg),'-',condition{iCondition},'.avg']);
            
            [signal,chan_names,variance, pnts, rate, xmin, xmax,nsweeps] = loadavg_bcl(eachAvgDir);
            baseNo   = numel(xmin*rate:0); 
            point    = baseNo+timeBin;
    
            ampdata(:,[size(eNo,2)*(iCondition-1)+1:size(eNo,2)*iCondition],iFile) = signal([point(1):point(end)],eNo);

        end
end


data.Tamplitude  = ampdata;

data.titleString = titleString;



% P1:100~200
% easyUncue = data.Tamplitude(:,[1:4],:);
% easyCue   = data.Tamplitude(:,[5:8],:);
% midUncue  = data.Tamplitude(:,[9:12],:);
% midCue    = data.Tamplitude(:,[13:16],:);
% hardUncue = data.Tamplitude(:,[17:20],:);
% hardCue   = data.Tamplitude(:,[21:24],:);

% P1_easyIOR = squeeze(mean(easyUncue-easyCue,2));
% P1_midIOR  = squeeze(mean(midUncue-midCue,2));
% P1_hardIOR = squeeze(mean(hardUncue-hardCue,2));


% N1:100~300
% easyUncue = data.Tamplitude(:,[1:6],:);
% easyCue   = data.Tamplitude(:,[7:12],:);
% midUncue  = data.Tamplitude(:,[13:18],:);
% midCue    = data.Tamplitude(:,[19:24],:);
% hardUncue = data.Tamplitude(:,[25:30],:);
% hardCue   = data.Tamplitude(:,[31:36],:);

% N1_easyIOR = squeeze(mean(easyUncue-easyCue,2));
% N1_midIOR  = squeeze(mean(midUncue-midCue,2));
% N1_hardIOR = squeeze(mean(hardUncue-hardCue,2));



% P2:200~400
% easyUncue = data.Tamplitude(:,[1:9],:);
% easyCue   = data.Tamplitude(:,[10:18],:);
% midUncue  = data.Tamplitude(:,[19:27],:);
% midCue    = data.Tamplitude(:,[28:36],:);
% hardUncue = data.Tamplitude(:,[37:45],:);
% hardCue   = data.Tamplitude(:,[46:54],:);

% easyIOR = easyUncue-easyCue,;
% midIOR  = midUncue-midCue;
% hardIOR = hardUncue-hardCue;