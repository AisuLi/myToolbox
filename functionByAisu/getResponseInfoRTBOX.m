function [] = getResponseInfoRTBOX(respIndexOrNameOrAddress)

 % Rev. by Yang Zhang Wed Apr  4 10:42:41 2018
 %     Soochow University, China
Evs = [];
iEv = 0;

while 1
	
   cEvt = CedrusResponseBox('GetButtons', respIndexOrNameOrAddress);


   if isempty(cEvt)
   		break;
   end 

   iEv      = iEv +1;
   Evs(iEv) = cEvt;
end 


actions             = [Evs(:).action];
bottons             = [Evs(:).botton]; 
respTimes           = [Evs(:).ptbtime]; 

bottons(~actions)   = [];
respTimes(~actions) = [];

