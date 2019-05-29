function  [amp hsns hcns] = AdaptiveRLS_foravg(filename,f0,chanNameList,isdetrend,orders,isAdaptiveRLS)

% Extract amplitudes of a specifical frequency via a RLS adaptive filter proposed by Tang and Norcia in 1995 
% Tang, Y., & Norcia, A. (1995). An adaptive filter for steady-state evoked responses. Electroencephalography and Clinical Neurophysiology, 96(3), 268â€?77.
% the data points should be much bigger than fs/(4*f0)
% 
% Written by Yang Zhang Tue Mar  1 17:37:28 2016
% Soochow University, China
% 
% 
% argins:
% 
% filename      [string]    filename of the avg file with full filepath , e.g., 'C:\uncued.avg';
% f0            [double]    a specifical frequency u want to extract from
% chanNameList  [cell of string]     channel names, e.g, {'FP1','FP2'}; Default is {'all'}
% isdetrend     [boolean]     if true will linear detrend the data





if ~exist('chanNameList','var')|isempty(chanNameList)
  chanNameList ={'all'};
end

if ~exist('isdetrend','var')|isempty(isdetrend)
  isdetrend = 0;
end


if ~exist('isAdaptiveRLS','var')
  isAdaptiveRLS = 1;
end

if ~exist('orders','var')|isempty(orders)
  orders = [1];
end
try
[signalData,chan_names,variance, pnts, rate, xmin, xmax,nsweeps]=loadavg_bcl(filename,chanNameList);

if isdetrend
  signalData = detrend(signalData);
end

if numel(orders) > 1
  fprintf('Will caculate more than 1 harmonics!\n');
  for iOrder = 1:numel(orders)
      fprintf('%d th harmonic: %3.1f Hz\n',orders(iOrder),f0*orders(iOrder));
  end

end
% orders = [1 2];% only caculates the main frequency instead of the harmonic frequencies

fs         = rate;

signalData = signalData'; % so that each row is corresponded to each chan


	%==============================
    %---- define the window length ------/
          for circles = 1:1000
              
              pin    = circles*(fs/f0);
              lambda = pin/(1+pin);
              
              if lambda>0.995
                  break
              end
          end
    %------------------------------------\
		  
		  display('======================  begin ================');

      sprintf('Window length: %d for target frequency : %5.2f/n',pin*1000/rate,f0);
		  display('    ');
      pin = round(pin);

      Hsn = zeros(pnts,numel(orders));
      Hcn = zeros(pnts,numel(orders));

	 for iSweep = 1:size(signalData,3)% define the cue type left vs. right
            for iChan = 1:size(signalData,1)
                
                y = signalData(iChan,:,iSweep);

        				if isAdaptiveRLS
                            for n = 1:pnts
                                for iOrder = 1:numel(orders)
                                    Hsn(n,iOrder) = sum((lambda.^(n-(1:n)).*y(1:n)).*sin(2*pi*orders(iOrder)*f0*(1:n)/fs))./...
                                        sum(lambda.^(n-(1:n)).*(sin(2*pi*orders(iOrder)*f0*(1:n)/fs)).^2);   %#ok<AGROW>
                                    
                                    Hcn(n,iOrder) = sum((lambda.^(n-(1:n)).*y(1:n)).*cos(2*pi*orders(iOrder)*f0*(1:n)/fs))./...
                                        sum(lambda.^(n-(1:n)).*(cos(2*pi*orders(iOrder)*f0*(1:n)/fs)).^2);   %#ok<AGROW>
                                end % iOrder
                            end% n
        				else 
                    for n = 1:pnts-pin+1
                      f        =abs(fft(y((1:pin)+n-1),512));
                      ampli(n) =f(round(f0/(fs/512))+1);	
                    end
        				end % isAdaptiveRLS
%                 CohHsn=conv(filters,Hsn);
%                 CohHcn=conv(filters,Hcn);
%                 CohHcn=CohHcn(2*pin+1:pnts);
%                 CohHsn=CohHsn(2*pin+1:pnts);
                if isAdaptiveRLS
                    % Hcn   =Hcn(pin+1:pnts);
                    % Hsn   =Hsn(pin+1:pnts);
                    ampli =(Hcn.^2+Hsn.^2).^0.5;  

                    ampli = mean(ampli,2); 
                end
					
                amp(iChan,:,iSweep) = ampli;  
                hcns(iChan,:,iSweep)= mean(Hcn,2);  
					      hsns(iChan,:,iSweep)= mean(Hsn,2);  

           end
     end
catch adaptiveRLS_ERROR;
    save adaptiveRLS_debug;
    
    rethrow(adaptiveRLS_ERROR);
    
end



        % for nfre=1:length(MarkerFre)
%              pnts=length(amp);
%              time=((pnts-1)*(1000/fs)-EpochRange(2))*-1:1000/fs:EpochRange(2);

%             for chanNum=1:size(amp,1)
%                 subplot(6,6,chanNum)
% %                 data=[];
%                 data=squeeze(amp(chanNum,:,:));
                
%                 chanidx(chanNum)=plot(time,data);
%                 hold on
%                 AmpliRange=[min(min(data))-abs(min(min(data)))*0.2,max(max(data))+abs(max(max(data)))*0.2];
%                 axis([time(1),time(length(time)),AmpliRange(1),AmpliRange(2)]);
%                 index(1)=fill([0,0,282,282],[AmpliRange(1),AmpliRange(2),AmpliRange(2),AmpliRange(1)],[0.6,0.6,0.6]);
%                 set(index(1),'EdgeAlpha',0,'FaceAlpha',0.4);
%                 index(2)=line([time(1),time(length(time))],[0,0],'LineWidth',1.5,'Color',[0 0 0]);
%                 index(3)=fill([1130,1130,282+1130,282+1130],[AmpliRange(1),AmpliRange(2),AmpliRange(2),AmpliRange(1)],[0,0,0.6]);
%                 set(index(3),'EdgeAlpha',0,'FaceAlpha',0.4);
%                 set(chanidx(chanNum),'XTick',[],'YTick',[]);
% %                 titleStr=(char(chan_names((chanNum-1)*2+nfre,:))); 
%                 title(name{chanNum});
%                 if chanNum==size(amp,1)
% 				le=legend('left','right');
%                 set(le,'')
%                 end
                
%             end
             
%              plot([1 1;2 2]);
              
%             printEvalStr=['print -dtiff -r600 Orig-sub-',name,'-',num2str(round(f0)),'hz'];
%             clf;
%             eval(printEvalStr);
%             close(figIndex);
        % end 
        
        
        
% duration=GetSecs-caculationStart;    
% disp(['The caculation procedure consume ',num2str(round(duration)),'secs  -_- ']);
    % cd (root);    
 
 
 

 
 
 
 
 
 
 
 
 
 
 
 
 

