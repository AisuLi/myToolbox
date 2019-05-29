        
function [ selecData] = cntSelection_bcl(transData,cnt,isEpoched,moreBinNum,stimTrigger)

% select usefulData(generally from firstTrigger To lastTriger),moreBinNum means u can get moreBinNum before firstTrigger and after lastTriger  
% written by Li Aisu, 04-17-2017

if ~exist('stimTrigger','var')||isempty(stimTrigger)
		stimTrigger = [1:40];
end 


if ~exist('isEpoched','var')||isempty(isEpoched)
	isEpoched = false;
end 




allTrigger = [cnt.event(:).stimtype];

if 	~isEpoched % cnt
        
	triggerLoc = find(allTrigger>0); % block is 0
	epochTime  = [cnt.event(1,triggerLoc(1)).offset-moreBinNum:cnt.event(1,triggerLoc(end)).offset+moreBinNum];
	selecData  = transData(:,epochTime);

else          %  epoch

	cnt.event(~ismember(allTrigger,stimTrigger))=[]; 
	eventTimePoint = [cnt.event(:).offset];
	
	for ievent = 1:size(eventTimePoint,2)
		iepochTime = [eventTimePoint(ievent)-moreBinNum(1):eventTimePoint(ievent)+moreBinNum(2)];
		selecData(:,:,ievent) = cnt.data(:,iepochTime);
	end

end

