function respIndexOrNameOrAddress = respDeviceConf(command,responseDeviceType,respIndexOrNameOrAddress)

global  osType_BCL 

%-- get OS info ----/
if isempty(osType_BCL)
    if IsWin
        osType_BCL = 1; % win
    elseif IsLinux
        osType_BCL = 2; % linux
    else
        osType_BCL = 3; % mac
    end
end
%-------------------\


if ~exist('respIndexOrNameOrAddress','var')
    respIndexOrNameOrAddress = [];
end 


% 
% if ~exist('restrictButtonIDinput','var')
%     restrictButtonIDinput = [];
% end 



switch lower(command)
    case 'open'
        
        switch responseDeviceType
            case 1 % keyboard
                respIndexOrNameOrAddress = -1;
            case 2 % Mouse;
                respIndexOrNameOrAddress = -1;
                % do nothing
            case 3 % gamepad
                switch osType_BCL
                    case 3 % mac 
                        %--- for PTB Methods ---/
                        respIndexOrNameOrAddress = GetGamepadIndices;

                        respIndexOrNameOrAddress = respIndexOrNameOrAddress(end);
                        %-----------------------\
                    case 2 % linux

                        if isempty(respIndexOrNameOrAddress)
                            % respIndexOrNameOrAddress = GetGamepadIndices('Microsoft X-Box One pad');

                            errorStr = ('For Linux joysticks, you should input either :');
                            errorStr = [errorStr,char(13),'1) the joystick name; '];
                            errorStr = [errorStr,char(13),'please run [noUsed, productNames] = GetGamepadIndices;'];
                            errorStr = [errorStr,char(13),'to get the joystick name (usually the last one)'];
                            errorStr = [errorStr,char(13),'Or:'];
                            errorStr = [errorStr,char(13),'2) the joystick index;'];
                            errorStr = [errorStr,char(13),'Please run findGamepadIndex_bcl to get the correct index!'];

                            error(errorStr);
                        end 

                        if ischar(respIndexOrNameOrAddress)
                            respIndexOrNameOrAddress = GetGamepadIndices(respIndexOrNameOrAddress);
                        end
                
                        [x,y,keyCode] = GetMouse([], respIndexOrNameOrAddress);
                        %-- there is a bug in SDL methods so we aborted it ---/
                        % if isempty(respIndexOrNameOrAddress)
                        %     %--- for PTB Methods ---/
                        %     respIndexOrNameOrAddress = 0;
                        %     %-----------------------\
                        % end
                        % joystickLinux('open',respIndexOrNameOrAddress,true);
                        % joystickLinux('query',respIndexOrNameOrAddress,1);   % after initializing, run it once
                        %-----------------------------------------------------\



                    case 1 % win
                        
                        if isempty(respIndexOrNameOrAddress)
                            respIndexOrNameOrAddress = 0;
                        end 
                        joystick_bcl(respIndexOrNameOrAddress);
                    otherwise
                        % do nothing
                end
                
            case 4 % cedrus's responsebox
                
                %------------- clear the cedrus drivers --------------/
                clear global ptb_cedrus_devices ptb_cedrus_drivertype;
                %-----------------------------------------------------\
                
                switch osType_BCL
                    case 2 % linux
                        if isempty(respIndexOrNameOrAddress)
                            respIndexOrNameOrAddress = '/dev/ttyUSB0';
                        end
                        %  Rev. by Yang Zhang Wed Apr  4 09:27:20 2018 Soochow University, China
                        % the 3rd and 4th parameter indicate lowbaudrate and doing time calibration respectively
                        respIndexOrNameOrAddress = CedrusResponseBox('Open',respIndexOrNameOrAddress,0,1);
                    case 1 % win
                        % if isempty(respIndexOrNameOrAddress)
                        %     respIndexOrNameOrAddress = '/dev/ttyUSB0';
                        % end
                        %----------- open the device --------------/
                        isSuccess = false;
                        
                        for iCom = 2:10  % suppose we have less than ten 
                            try
                                respIndexOrNameOrAddress = CedrusResponseBox('Open',['COM',num2str(iCom)]);
                                
                                isSuccess         = true;
                            catch
                                % do nothing
                            end
                            
                            if isSuccess
                                break;
                            end
                        end
                        %------------------------------------------\
                    case 3 % mac
                        if isempty(respIndexOrNameOrAddress)
                            respIndexOrNameOrAddress = '/dev/cu.usbserial-FTDI125ZX9';
                        end
                        respIndexOrNameOrAddress = CedrusResponseBox('Open',respIndexOrNameOrAddress);
                        % respIndexOrNameOrAddress = CedrusResponseBox('Open','/dev/ttyUSB0'); to be confirmed
                    otherwise
                        % do nothing
                end
            otherwise
                error('unsupported response device, responseDeviceType should be of [1 2 3 4]!');
        end
        
    case 'close'
        
        %---- close gamepad OR cedrus response box -----/
        switch responseDeviceType
            case 4
                CedrusResponseBox('CloseAll');
            case 3

                switch osType_BCL
                    case 2 % linux
                    if isempty(respIndexOrNameOrAddress)
                        respIndexOrNameOrAddress = 0;
                    end 
                    %--- for yang's methods ----/
                    % joystickLinux('close',respIndexOrNameOrAddress,true);
                    %---------------------------\
                otherwise

                end
            case 2
                % do nothing
            case 1
                % do nothing
            otherwise
                % do nothing
        end % switch
        %----------------------------------------------\
        
        respIndexOrNameOrAddress = [];
        %-- clear the osType_BCL ----/
        clear global osType_BCL;
        %----------------------------\
        
    otherwise
        error('command should be of "Open" or "Close"!');
end