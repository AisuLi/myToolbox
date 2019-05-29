function [eventBins,foundMarkers,epochedTriggerBins,beplotedX,epochedData,randColor] = monitorTest2(cntfilename,nframes,ifi,iLoc)


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



bePlotXs = [eventBins;eventBins;nan(size(eventBins))];
bePlotYs = [min(data(:))*ones(size(eventBins));max(data(:))*ones(size(eventBins));nan(size(eventBins))];

% line(bePlotXs(:),bePlotYs(:),'Color',[0 0 1]);

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


% line(bePlotXs(:),bePlotYs(:),'Color',[0 1 0]);

bePlotXs = [foundMarkersOff;foundMarkersOff;nan(size(foundMarkersOff))];
bePlotYs = [min(data(:))*ones(size(foundMarkersOff));max(data(:))*ones(size(foundMarkersOff));nan(size(foundMarkersOff))];


% line(bePlotXs(:),bePlotYs(:),'Color',[0.5 0.5 0.5]);


% lh = legend('rawData','trigger','thresholdOn','thresholdOff');
% set(lh,'box','off');

% xlabel('time bins');

% ylabel('voltages (uV)');

% xlim([max(eventBins(1) - 1000,1),min(foundMarkersOff(end)+1000,numel(data))]);
% hold off;




% subplot(3,1,2);

% plot([(foundMarkers - eventBins)',(foundMarkersOff - foundMarkers)']);

% xlabel('testing nums');
% ylabel('deviations (bins)');

% lh = legend('threshold on - trigger','threshold off - on');
% set(lh,'box','off');

% ylim([min(mean(foundMarkers - eventBins),mean(foundMarkersOff - foundMarkers))*0.2,max(mean(foundMarkers - eventBins),mean(foundMarkersOff - foundMarkers))*1.2]);


% randColor=rand(50,3);

randColor = [0.0606800204231480	0.902009465824431	0.962599013218848
0.0907508274678313	0.195184184542984	0.0875689210675951
0.805585641902073	0.897667932373422	0.255660876251761
0.546980919566268	0.821605771707393	0.412617762743869
0.895124242734885	0.665945381323852	0.703928253296993
0.186858761589189	0.656676024099073	0.747978484978250
0.832460573908841	0.535541779825160	0.0374367222841506
0.0626397734212979	0.148039406542768	0.781343841536334
0.600029557848887	0.615372649538667	0.316967910740484
0.0600721266080804	0.373914880657112	0.414421711994423
0.233532910627902	0.438732996175382	0.0879032229549459
0.219706353145267	0.599164947022138	0.349375779722459
0.259250380460217	0.877968834495135	0.129423561006250
0.599285161743039	0.115913419416761	0.304327202838603
0.694915398370391	0.985739540883822	0.946217354716618
0.922127964308747	0.857254469293949	0.770128099543164
0.655575466925251	0.441606147265157	0.821307287538090
0.633643357559537	0.831619780151143	0.114782531011690
0.401798139293280	0.292644842379715	0.203310651903552
0.357841996705088	0.510935901808601	0.592598263314838
0.453987599494908	0.751237166041510	0.0106912757564572
0.838933580444433	0.380326632265418	0.271813176584230
0.390013418605796	0.266266942314067	0.660901376287224
0.889675973309932	0.262603921217760	0.200749943817415
0.811448379649953	0.582624380387433	0.685722449374200
0.822251657461199	0.443120531084340	0.107285913105537
0.864895217749837	0.446478632441557	0.446994972411219
0.338744647324380	0.978698027425513	0.127039953137372
0.0632638852819618	0.791360945545610	0.519774582833318
0.266862545953493	0.211472445653468	0.256774164943503
0.0925699514987907	0.948570817980897	0.159930772426874
0.671309991409510	0.0585511461756274	0.421714971604656
0.944915448571485	0.0239480487550577	0.558692031556859
0.656308816159899	0.208484420420066	0.300660109775817
0.416360101765498	0.294285371021574	0.735579040615798
0.241606360052552	0.366028343025171	0.340174651164455
0.350991432832415	0.850064183114790	0.0446923134250257
0.699740414476989	0.941759228807027	0.408185331322768
0.488717014809074	0.0637994243281149	0.322181351065950
0.820372499255613	0.409576227270505	0.358314252710203
0.329963318321199	0.947448788352265	0.0542195484063197
0.0940674616885745	0.222860614992551	0.748507645380094
0.831761747601366	0.0475658698086285	0.0797376003902225
0.426677551738374	0.00722300874743898	0.292543038806435
0.256949691340334	0.219968027004846	0.416568365953175
0.605892642557118	0.673440894019482	0.699781703746440
0.577125304829104	0.131854989652384	0.229094331408131
0.961446714134265	0.704347394731617	0.798982598411030
0.0633784441506314	0.384885483429272	0.760793310432207
0.719434069219439	0.362163503410741	0.340752818464625];


transparent = 0.4;

subplot(3,3,iLoc);

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

% baseData = mean(epochedData([1:20],:),1);
% epochedData = epochedData-baseData;

for iTrial=1:numel(foundMarkers)
	H(iTrial)=plot(beplotedX,epochedData(:,iTrial));
	set(H(iTrial),'lineWidth',1,'color',[randColor(iTrial,:) transparent]);
	hold on;
end

box off;



meanEpochedData = mean(epochedData,2);

plot(beplotedX,meanEpochedData,'lineWidth',1,'color',[0.7 0 0]);

xlabel('time bins');
ylabel('epoched data (uV)');
ylim([-5,85]);
set(gca,'ytick',[0:10:80]);
xlim([-50,beplotedX(end)]);
set(gca,'xtick',[0:100:beplotedX(end)]);
hold on





% nframes=4;
% ifi=1000/120;

% cnt=dir('*.cnt');
% for i=1
% 	cntfilename=cnt(i).name;
% 	[eventBins,foundMarkers,epochedData] = monitorTest(cntfilename,nframes,ifi)
% end
