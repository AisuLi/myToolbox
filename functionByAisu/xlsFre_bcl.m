
function [value] = xlsFre_bcl(fileName,whichSheet,row,groupNum)

% written by Li Aisu 

if ~exist('whichSheet','var')||isempty(whichSheet)
	whichSheet = 1;
end 

if ~exist('groupNum','var')||isempty(groupNum)
	groupNum = 100;
end 


[NUM] = xlsread(fileName,whichSheet,row);

j=1;
for i = 1:numel(NUM)
   if NUM(i)>0
  	 fre(j) = NUM(i);
  	 j= j+1;
   end
end

fre(:) = roundn(fre(:),-1);
value = tabulate(fre);

hist(fre,groupNum,'EdgeColor','none');
xlim([0 30]);
box off;
set(gca,'xtick',[5 7 8.6 10 12 15 17 20 23 25 27.8]);        
