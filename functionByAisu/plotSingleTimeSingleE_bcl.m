
function plotSingleTimeSingleE_bcl(matFileDirectory,condition,conditionNum,figStr)

%% after loopAllTimeE_bcl,plotResult
%% written by Li Aisu 


load(matFileDirectory);

cue(:,[1:64]) = spss.cueSig(1:end-1,[1:64]);
% cue(:,[65,66,68,69,71,72])= nan;
% cue(:,[70,73])= spss.cueSig(1:end-1,[65:66]);


level(:,[1:64]) = spss.levelSig(1:end-1,[1:64]);
% level(:,[65,66,68,69,71,72])= nan;
% level(:,[70,73])= spss.levelSig(1:end-1,[65:66]);


levelcue(:,[1:64]) = spss.levelcueSig(1:end-1,[1:64]);
% levelcue(:,[65,66,68,69,71,72])= nan;
% levelcue(:,[70,73])= spss.levelcueSig(1:end-1,[65:66]);



set(gcf,'Position',[100,100,1800,800], 'color','w');

% condition={'cue','level','levelCue'};


if conditionNum==1
		H1=surf(cue);
elseif conditionNum==2
		H1=surf(level);
else
		H1=surf(levelcue);
end


set(H1,'LineStyle','none');
view(90,90);
grid off;
set(gca,'ytick',[0:50:1000]);
set(gca,'ytickLabel',[-200:50:800]);

title(condition{conditionNum}) ; % main effect of cueing interaction of cue*level
standard_eName={'FP1','FPZ','FP2','AF3','AF4',...
            'F7','F5','F3','F1','FZ','F2','F4','F6','F8',...
            'FT7','FC5','FC3','FC1','FCZ','FC2','FC4','FC6','FT8',...
            'T7','C5','C3','C1','CZ','C2','C4','C6','T8','M1',...
            'TP7','CP5','CP3','CP1','CPZ','CP2','CP4','CP6','TP8','M2',...
            'P7','P5','P3','P1','PZ','P2','P4','P6','P8',...
            'PO7','PO5','PO3','POZ','PO4','PO6','PO8',...
            'CB1','O1','OZ','O2','CB2','','','none','','','VEO','','','HEO'};
set(gca,'xtick',[1:3:64]);        
set(gca,'xTickLabel', standard_eName([1:3:64])); % [2,10,19,28,33,28,43,48,56,62]

hold on
saveas(gca,[condition{conditionNum},figStr,'.bmp']);
close all;

