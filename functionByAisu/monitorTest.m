function [eventBins,foundMarkers,epochedData] = monitorTest(cntfilename,nframes,ifi)


if ~exist('nframes','var')||isempty(nframes)
	nframes = 4;
end

% if ~exist('polarity','var')||isempty(polarity)
% 	polarity = 1;
% end

if ~exist('ifi','var')||isempty(ifi)
	ifi = 1000/120;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              begin 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
smoothRangeBins = 8;

c         = loadcnt_bcl(cntfilename);

triggerCode = [c.event(:).stimtype];

c.event(~ismember(triggerCode,[100]))=[]; % filtered out the no stim triggers


eventBins = [c.event(:).offset];

tempData = c.data(eventBins(1):eventBins(end));

%----- get the polarity of the Data ---------/
if abs(max(tempData(:)))>abs(min(tempData(:)))
	polarity =  1;
else
	polarity = -1;
end
%--------------------------------------------\

data      = c.data*polarity;
stimDur   = round(ifi*nframes*c.header.rate/1000);



[Noused,filenameOnly] = fileparts(cntfilename);



threshold       = max(data(:))*0.4;

subThresholdLow = threshold*0.9;
subThresholdUp  = threshold*1.1;

blindWins       = round((nframes - 1)*ifi*c.header.rate/1000) + 15;


figure;
set(gcf,'Name',filenameOnly);

subplot(3,1,1);
plot(data,'r');

hold on;


bePlotXs = [eventBins;eventBins;nan(size(eventBins))];
bePlotYs = [min(data(:))*ones(size(eventBins));max(data(:))*ones(size(eventBins));nan(size(eventBins))];

line(bePlotXs(:),bePlotYs(:),'Color',[0 0 1]);



markerIdx = data>threshold;
markerIdx3 = data>subThresholdLow;
markerIdx2 = data>subThresholdUp;


startBin = c.event(1).offset - blindWins;
endBin = c.event(end).offset + 1000;



foundMarkers = [];

foundMarkersOff = [];


% for iBin = max(startBin,blindWins+1):endBin
iBin = max(startBin,blindWins+1);

iLoop = 1;

while iBin < endBin 
    
    isIncreased = false;
%     
%      if ~mod(iLoop,1000)
%      	fprintf('%d\n',iBin);
%      end
	if all(~markerIdx2(iBin - smoothRangeBins:iBin - 1))&&markerIdx(iBin)
		foundMarkers = [foundMarkers, iBin];

		iBin = iBin + blindWins; % skipping the blind wins ranges
        
        fprintf('markerOn : %d -->%d\n',foundMarkers(end),iBin);
        
        isIncreased = true;
	end

	if all(markerIdx3(iBin - smoothRangeBins:iBin - 1))&&~markerIdx(iBin)
		foundMarkersOff = [foundMarkersOff, iBin];

		iBin = iBin + round(blindWins/2); % skipping the blind wins ranges
        
        fprintf('markerOff: %d -->%d\n',foundMarkersOff(end),iBin);
        isIncreased = true;
	end

	if ~isIncreased
		iBin = iBin +1;
	end

	iLoop = iLoop +1;
end % while



% end


bePlotXs = [foundMarkers;foundMarkers;nan(size(foundMarkers))];
bePlotYs = [min(data(:))*ones(size(foundMarkers));max(data(:))*ones(size(foundMarkers));nan(size(foundMarkers))];

line(bePlotXs(:),bePlotYs(:),'Color',[0 1 0]);


bePlotXs = [foundMarkersOff;foundMarkersOff;nan(size(foundMarkersOff))];
bePlotYs = [min(data(:))*ones(size(foundMarkersOff));max(data(:))*ones(size(foundMarkersOff));nan(size(foundMarkersOff))];

line(bePlotXs(:),bePlotYs(:),'Color',[0.5 0.5 0.5]);


lh = legend('rawData','trigger','thresholdOn','thresholdOff');
set(lh,'box','off');

xlabel('time bins');

ylabel('voltages (uV)');

xlim([max(eventBins(1) - 1000,1),min(foundMarkersOff(end)+1000,numel(data))]);
hold off;


subplot(3,1,2);

plot([(foundMarkers - eventBins)',(foundMarkersOff - foundMarkers)']);

xlabel('testing nums');
ylabel('deviations (bins)');

lh = legend('threshold on - trigger','threshold off - on');
set(lh,'box','off');

ylim([min(mean(foundMarkers - eventBins),mean(foundMarkersOff - foundMarkers))*0.2,max(mean(foundMarkers - eventBins),mean(foundMarkersOff - foundMarkers))*1.2]);







subplot(3,1,3);

epochedTriggerBins = eventBins - foundMarkers;

beplotedX          = [epochedTriggerBins;epochedTriggerBins;nan(size(epochedTriggerBins))];

beplotedY          = repmat([-5;10;NaN],1,numel(epochedTriggerBins));

plot(beplotedX(:),beplotedY(:));


hold on;

baselineRange = round(stimDur*0.1);

epochedData = zeros(stimDur + baselineRange+1,numel(foundMarkers));


for iEvent = 1:numel(foundMarkers)
	epochedData(:,iEvent) = data((foundMarkers(iEvent)-baselineRange):(foundMarkers(iEvent)+stimDur)  );
end


beplotedX = -baselineRange:stimDur;


plot(beplotedX,epochedData);

xlabel('time bins');
ylabel('epoched data (uV)');

ylim([-10,50]);
xlim([min(epochedTriggerBins - 5),beplotedX(end)]);



hold off;
