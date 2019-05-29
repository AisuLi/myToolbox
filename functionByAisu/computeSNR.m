		
function [SNR] = computeSNR(Power,snrskipNum)  

if ~exist('snrskipNum','var')||isempty(snrskipNum)
		snrskipNum = 1;
end 


SNR = zeros(size(Power));

for ie  = 1:size(Power,2)       
	for iFrebin = snrskipNum+2:size(Power,1)-(snrskipNum+1)
		SNR(iFrebin,ie) = Power(iFrebin,ie)/mean(Power([iFrebin-(snrskipNum+1),iFrebin+(snrskipNum+1)],ie));
	end
end	