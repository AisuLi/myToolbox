function [EEG] = RESS_Yang_su(epochedFilename,fftType,targetFreqs,cond2use,caculateRange,isPlot)

%% RESS empirical example
% this script accompanies the manuscript titled 
%    "Rhythmic entrainment source separation: 
%     Optimizing analyses of neural responses to rhythmic sensory stimulation"
%   from Mike X Cohen and Rasa Gulbinaite

% 1) We encourage working through and adapting this script; 
%    it is not designed to be used as a function.
% 2) You will need the following additional files: 
%              filterFGx.m
%              RESS_EEG_data.mat
%              eeglab toolbox (only for topographical plotting, otherwise just comment out those lines)
% question? -> mikexcohen@gmail.com or rasa.gulbinaite@gmail.com

% clear

 % Rev. by Yang Zhang Mon Jan 22 16:50:42 2018
 %     Soochow University, China

if nargin < 2
    error('we need at least two prameters: opochedFilename and targetFreqs!');
end 

if ~exist('cond2use','var')||isempty(cond2use)
    cond2use = [];
end 


if ~exist('caculateRange','var')||isempty(caculateRange)
    caculateRange = [];
end 


if ~exist('isPlot','var')||isempty(isPlot)
    isPlot = 1;
end 

% a secret parameter for yang
isEOGincluded = 0;

% written by Yang Zhang Mon Dec 25 23:29:41 2017
%      Soochow University, China
%% specify parameters

% the experiment had 10/17 Hz flicker. You can also try harmonics.
% peakfreq1 = 10; % hz
% peakfreq2 = 17; % hz

% used for 'best-electrode' analyses
electrode1 = 'oz';
electrode2 = 'o2';



% parameters for RESS:
peakwidt  = .5; % FWHM at peak frequency
neighfreq = 1;  % distance of neighboring frequencies away from peak frequency, +/- in Hz
neighwidt = 1;  % FWHM of the neighboring frequencies

%% load in the data

% contains EEG data (in eeglab format) and leadfield information, which is
% used to create the anatomical estimations.
load scan65chansLF;
% load RESS_EEG_data

% IMPORTANT NOTE ABOUT ANATOMICAL LOCALIZATION!!! This leadfield is
% specific to this particular electrode montage. You cannot use this
% leadfield with your EEG data!! If you want to use the same code on your
% data, you will need to create a leadfield, e.g., via the Brainstorm
% toolbox.

%% condition markers

% first number is condition code; thereafter is a brief description:

% visual conditions only
% 11:  10 Hz outside 17 Hz inside
% 12:  17 Hz outside 10 Hz inside
% 13:  16 Hz outside 17 Hz inside
% 14:  17 Hz outside 16 Hz inside

% multisensory conditions. Visual stimulation was 15 Hz in all conditions
% 31:  left  hemifield / left  ear - 25 Hz
% 32:  right hemifield / right ear - 25 Hz
% 33:  left  hemifield / left  ear - 40 Hz
% 34:  right hemifield / right ear - 40 Hz

if ischar(epochedFilename)
    [filepath,filenameOnly] = fileparts(epochedFilename);

    if ismac 
        folderDiv = '/';
    else
        folderDiv = '\';
    end 

    EEG = pop_loadeeg([filenameOnly,'.eeg'], [filepath,folderDiv], 'all','all','all','all','auto');
else
    EEG = epochedFilename;
end

% filter out the eog channels added by yang 



dataChans = lower({EEG.chanlocs.labels});

% excluding eog chans from leadfield
ressChans  = lower({lf.chanlocs.labels});


if ~isEOGincluded
	eogChanIdx = ismember(ressChans,{'heo','veo'})|ismember(ressChans,{'heog','veog'});

	lf.Gain(eogChanIdx,:)   = [];
	lf.chanlocs(eogChanIdx) = [];
	ressChans(eogChanIdx)   = [];

	% excluding eog chans from the data
	eogChanIdx = ismember(dataChans,{'heo','veo'})|ismember(dataChans,{'heog','veog'});

	EEG.data(eogChanIdx,:,:) = [];
	EEG.chanlocs(eogChanIdx) = [];
    EEG.nbchan = EEG.nbchan - sum(eogChanIdx);
	dataChans(eogChanIdx)    = [];
end 

[isRessChanInData,sortedOrderbyRess] = ismember(ressChans,dataChans);

if any(diff(sortedOrderbyRess)~=1)

	fprintf('LD   data\n');

	for iChan = 1:min([numel(ressChans),numel(dataChans)])
		fprintf('%s  %s\n',ressChans{iChan},ressChans{dataChans});
	end 
	error('chan orders are not corresponded between the data and the leadfield!');
end 

% find which epochs contain which markers
epochvect = zeros(1,EEG.trials);

epochvect = [EEG.epoch(:).type];
% for i=1:EEG.trials
%     [~,t0] = min(abs(cell2mat(EEG.epoch(i).eventlatency)));
%     epochvect(i) = EEG.epoch(i).eventtype{t0};
% end

%% Start RESS for condition "1"

% FFT parameters
nfft = ceil( EEG.srate/.1 ); % .1 Hz resolution
% tidx = dsearchn(EEG.times',[500 10000]'); % we use data from .5-10 seconds, to avoid using the stimulus transient

% checking the caculateRange Parameter
if isempty(caculateRange)
    caculateRange = repmat([EEG.times(1) EEG.times(end)],numel(targetFreqs),1);
else
    switch size(caculateRange,1)
    case 1
        caculateRange = repmat(caculateRange,numel(targetFreqs),1);
    case  numel(targetFreqs)
        % do nothing
    otherwise
        error('size(caculateRange,1) should be of 1 or numel(targetFreqs)!!!');
    end 
end 


if isempty(caculateRange)

    tidx = zeros(numel(targetFreqs),2);

    for iFre = 1:numel(targetFreqs)
        tidx(iFre,:) = dsearchn(EEG.times',[EEG.times(1) EEG.times(end)]');
    end
else
    for iFre = 1:numel(targetFreqs)
        tidx(iFre,:) = dsearchn(EEG.times',[caculateRange(iFre,1) caculateRange(iFre,end)]');
    end
end 

if isempty(cond2use)
    cond2use = unique(epochvect);
end 


% extract EEG data
data  = EEG.data(:,:,ismember(epochvect,cond2use));

if strcmp(fftType,'avgFFT')
    data = mean(data,3);
end


hz    = linspace(0,EEG.srate,nfft);




for iFreq = 1:numel(targetFreqs) 

    dataX = mean(abs( fft(data(:,tidx(iFreq,1):tidx(iFreq,2),:),nfft,2)/diff(tidx(iFreq,:)) ).^2,3);

    ssvep(iFreq).targetFreq = targetFreqs(iFreq);
    ssvep(iFreq).range      = caculateRange(iFreq,:);
    % compute covariance matrix at peak frequency
    fdatAt = filterFGx(data,EEG.srate,targetFreqs(iFreq),peakwidt);
    fdatAt = reshape( fdatAt(:,tidx(iFreq,1):tidx(iFreq,2),:), EEG.nbchan,[] );
    fdatAt = bsxfun(@minus,fdatAt,mean(fdatAt,2));
    covAt  = (fdatAt*fdatAt')/diff(tidx(iFreq,:));

    ssvep(iFreq).fdatAt = fdatAt ;
    ssvep(iFreq).covAt  = covAt  ;
   
    % compute covariance matrix for lower neighbor
    fdatLo = filterFGx(data,EEG.srate,targetFreqs(iFreq)+neighfreq,neighwidt);
    fdatLo = reshape( fdatLo(:,tidx(iFreq,1):tidx(iFreq,2),:), EEG.nbchan,[] );
    fdatLo = bsxfun(@minus,fdatLo,mean(fdatLo,2));
    covLo  = (fdatLo*fdatLo')/diff(tidx(iFreq,:));

    ssvep(iFreq).fdatLo = fdatLo;
    ssvep(iFreq).covLo  = covLo ;
    % compute covariance matrix for upper neighbor
    fdatHi = filterFGx(data,EEG.srate,targetFreqs(iFreq)-neighfreq,neighwidt);
    fdatHi = reshape( fdatHi(:,tidx(iFreq,1):tidx(iFreq,2),:), EEG.nbchan,[] );
    fdatHi = bsxfun(@minus,fdatHi,mean(fdatHi,2));
    covHi  = (fdatHi*fdatHi')/diff(tidx(iFreq,:));

    ssvep(iFreq).fdatHi = fdatHi;
    ssvep(iFreq).covHi  = covHi ;
    % perform generalized eigendecomposition. This is the meat & potatos of RESS
    [evecs,evals] = eig(covAt,(covHi+covLo)/2);
    [~,comp2plot] = max(diag(evals)); % find maximum component
    evecs         = bsxfun(@rdivide,evecs,sqrt(sum(evecs.^2,1))); % normalize vectors (not really necessary, but OK)

    % extract components and force sign
    % maps = inv(evecs'); % get maps (this is fine for full-rank matrices)
    maps    = covAt * evecs / (evecs' * covAt * evecs); % this works either way
    [~,idx] = max(abs(maps(:,comp2plot))); % find biggest component
    maps    = maps * sign(maps(idx,comp2plot)); % force to positive sign

    ssvep(iFreq).maps = maps;
    % reconstruct RESS component time series
    ress_ts = zeros(EEG.pnts,size(data,3));
    for ti=1:size(data,3)
        ress_ts(:,ti) = evecs(:,comp2plot)'*squeeze(data(:,:,ti));
    end

    ssvep(iFreq).ress_ts = ress_ts;
    %% compute SNR spectrum

    ressx = mean(abs( fft(ress_ts(tidx(iFreq,1):tidx(iFreq,2),:),nfft,1)/diff(tidx(iFreq,:)) ).^2,2);
    elecx = dataX(strcmpi(electrode1,{EEG.chanlocs.labels}),:,:);


    [snrR,snrE] = deal(zeros(size(hz)));
    skipbins    =  5; % .5 Hz, hard-coded!
    numbins     = 20+skipbins; %  2 Hz, also hard-coded!

    % loop over frequencies and compute SNR
    for hzi=numbins+1:length(hz)-numbins-1
        numer = ressx(hzi);
        denom = mean( ressx([hzi-numbins:hzi-skipbins hzi+skipbins:hzi+numbins]) );
        snrR(hzi) = numer./denom;
        
        numer = elecx(hzi);
        denom = mean( elecx([hzi-numbins:hzi-skipbins hzi+skipbins:hzi+numbins]) );
        snrE(hzi) = numer./denom;
    end

    ssvep(iFreq).snrR = snrR;
    ssvep(iFreq).snrE = snrE;
    ssvep(iFreq).hz   = hz;
    %% estimate anatomical distribution

    beamsource = zeros(3,size(lf.Gain,3));
    for voxi=1:size(lf.Gain,3)
        for ori=1:3
            beamsource(ori,voxi) = corr( maps(:,comp2plot) , squeeze(lf.Gain(:,ori,voxi)) );
        end
    end
    forw_brain = sqrt(sum(beamsource.^2));


    ssvep(iFreq).beamsource = beamsource;
    ssvep(iFreq).forw_brain = forw_brain;


    if isPlot
        %% some plotting...

        figure(1), clf
        xlim = [3 45];

        subplot(numel(targetFreqs),4,(iFreq-1)*numel(targetFreqs)+1)
        plot(hz,snrR,'ro-','linew',1,'markersize',5,'markerface','w')
        hold on
        plot(hz,snrE,'ko-','linew',1,'markersize',5,'markerface','w')
        set(gca,'xlim',xlim)
        axis square
        xlabel('Frequency (Hz)'), ylabel('SNR')
        legend({'RESS';electrode1})

        subplot(numel(targetFreqs),4,(iFreq-1)*numel(targetFreqs)+2)
        map2plot = maps(:,comp2plot);
        topoplot(map2plot./max(map2plot),lf.chanlocs,'maplimits',[-.7 .7],'numcontour',0,'conv','on','electrodes','off','shading','interp');
        title([ 'RESS for ' num2str(targetFreqs(iFreq)) ' Hz' ])

        subplot(numel(targetFreqs),4,(iFreq-1)*numel(targetFreqs)+3)
        map2plot = dataX(:,dsearchn(hz',targetFreqs(iFreq)));
        topoplot(map2plot./max(map2plot),lf.chanlocs,'maplimits',[-.7 .7],'numcontour',0,'conv','on','electrodes','off','emarker2',{find(strcmpi({EEG.chanlocs.labels},electrode1)) 'o' 'w' 4},'shading','interp');
        title([ 'RESS for ' num2str(targetFreqs(iFreq)) ' Hz' ])
        title([ 'Electrode power at ' num2str(targetFreqs(iFreq)) ' Hz' ])



        subplot(numel(targetFreqs),4,(iFreq-1)*numel(targetFreqs)+4)
        % convert scalar colors to RGB
        cmap   = colormap;
        cid    = sum(bsxfun(@gt,forw_brain,linspace(0,1,length(cmap))'));
        newcol = cmap(cid,:);

        % threshold map at >median
        thresh = median(forw_brain) + std(forw_brain)*0;

        newcol(forw_brain<thresh,:) = repmat(.7,sum(forw_brain<thresh),3);

        h = patch(lf.p);
        set(h,'FaceVertexCData',newcol,'facecolor','flat','edgecolor','none')
        axis image, axis off, rotate3d on

        light, lightangle(-90,25), material shiny
        view([-30 30])
        % zoom(1.4) % computer-dependent

    end % isPlot

end % iFreq

% EEG.data = [];
EEG.ssvep = ssvep;








% %% now repeat for the other frequency/condition

% % compute covariance matrix at peak frequency
% fdatAt = filterFGx(data,EEG.srate,peakfreq2,peakwidt);
% fdatAt = reshape( fdatAt(:,tidx(iFreq,1):tidx(iFreq,2),:), EEG.nbchan,[] );
% fdatAt = bsxfun(@minus,fdatAt,mean(fdatAt,2));
% covAt  = (fdatAt*fdatAt')/diff(tidx(iFreq,:));

% % compute covariance matrix for lower neighbor
% fdatLo = filterFGx(data,EEG.srate,peakfreq2+neighfreq,neighwidt);
% fdatLo = reshape( fdatLo(:,tidx(iFreq,1):tidx(iFreq,2),:), EEG.nbchan,[] );
% fdatLo = bsxfun(@minus,fdatLo,mean(fdatLo,2));
% covLo  = (fdatLo*fdatLo')/diff(tidx(iFreq,:));

% % compute covariance matrix for upper neighbor
% fdatHi = filterFGx(data,EEG.srate,peakfreq2-neighfreq,neighwidt);
% fdatHi = reshape( fdatHi(:,tidx(iFreq,1):tidx(iFreq,2),:), EEG.nbchan,[] );
% fdatHi = bsxfun(@minus,fdatHi,mean(fdatHi,2));
% covHi  = (fdatHi*fdatHi')/diff(tidx(iFreq,:));

% % perform generalized eigendecomposition. This is the meat & potatos of RESS
% [evecs,evals] = eig(covAt,(covHi+covLo)/2);
% [~,comp2plot] = max(diag(evals)); % find maximum component
% evecs = bsxfun(@rdivide,evecs,sqrt(sum(evecs.^2,1))); % normalize vectors (not really necessary, but OK)

% % extract components and force sign
% % maps = inv(evecs'); % get maps (this is fine for full-rank matrices)
% maps = covAt * evecs / (evecs' * covAt * evecs); % this works either way
% [~,idx] = max(abs(maps(:,comp2plot))); % find biggest component
% maps = maps * sign(maps(idx,comp2plot)); % force to positive sign


% % reconstruct RESS component time series
% ress_ts2 = zeros(EEG.pnts,size(data,3));
% for ti=1:size(data,3)
%     ress_ts2(:,ti) = evecs(:,comp2plot)'*squeeze(data(:,:,ti));
% end

% %% compute SNR spectrum

% ressx = mean(abs( fft(ress_ts2(tidx(iFreq,1):tidx(iFreq,2),:),nfft,1)/diff(tidx(iFreq,:)) ).^2,2);
% elecx = dataX(strcmpi(electrode2,{EEG.chanlocs.labels}),:,:);


% [snrR,snrE] = deal(zeros(size(hz)));
% skipbins =  5; % .5 Hz, hard-coded!
% numbins  = 20+skipbins; %  2 Hz, also hard-coded!

% % loop over frequencies and compute SNR
% for hzi=numbins+1:length(hz)-numbins-1
%     numer = ressx(hzi);
%     denom = mean( ressx([hzi-numbins:hzi-skipbins hzi+skipbins:hzi+numbins]) );
%     snrR(hzi) = numer./denom;
    
%     numer = elecx(hzi);
%     denom = mean( elecx([hzi-numbins:hzi-skipbins hzi+skipbins:hzi+numbins]) );
%     snrE(hzi) = numer./denom;
% end

% %% estimate anatomical distribution

% beamsource = zeros(3,size(lf.Gain,3));
% for voxi=1:size(lf.Gain,3)
%     for ori=1:3
%         beamsource(ori,voxi) = corr( maps(:,comp2plot) , squeeze(lf.Gain(:,ori,voxi)) );
%     end
% end
% forw_brain = sqrt(sum(beamsource.^2));

% %% some plotting...

% subplot(245)
% plot(hz,snrR,'ro-','linew',1,'markersize',5,'markerface','w')
% hold on
% plot(hz,snrE,'ko-','linew',1,'markersize',5,'markerface','w')
% set(gca,'xlim',xlim)
% axis square
% xlabel('Frequency (Hz)'), ylabel('SNR')
% legend({'RESS';electrode2})

% subplot(246)
% map2plot = maps(:,comp2plot);
% topoplot(map2plot./max(map2plot),EEG.chanlocs,'maplimits',[-.7 .7],'numcontour',0,'conv','on','electrodes','off','shading','interp');
% title([ 'RESS for ' num2str(peakfreq2) ' Hz' ])

% subplot(247)
% map2plot = dataX(:,dsearchn(hz',peakfreq2));
% topoplot(map2plot./max(map2plot),EEG.chanlocs,'maplimits',[-.7 .7],'numcontour',0,'conv','on','electrodes','off','emarker2',{find(strcmpi({EEG.chanlocs.labels},electrode2)) 'o' 'w' 4},'shading','interp');
% title([ 'RESS for ' num2str(peakfreq2) ' Hz' ])
% title([ 'Electrode power at ' num2str(peakfreq2) ' Hz' ])



% subplot(248)
% % convert scalar colors to RGB
% cmap = colormap;
% cid = sum(bsxfun(@gt,forw_brain,linspace(0,1,length(cmap))'));
% newcol = cmap(cid,:);

% % threshold map at >median
% thresh = median(forw_brain) + std(forw_brain)*0;
% newcol(forw_brain<thresh,:) = repmat(.7,sum(forw_brain<thresh),3);

% h = patch(lf.p);
% set(h,'FaceVertexCData',newcol,'facecolor','flat','edgecolor','none')
% axis image, axis off, rotate3d on

% light, lightangle(-90,25), material shiny
% view([-30 30])
% zoom(1.4) % computer-dependent

