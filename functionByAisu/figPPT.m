function [opt,ppt] = figPPT(pptTitleName)

    opt.fontSize      =18;
    opt.fontName      ='Calibri (Body)';
    opt.textAlignment ='ppAlignCenter';
    opt.picPosition   =[];
    

    ppt=PPT2007_zoe();% open a PPT

    
    ppt=ppt.addTitleSlide([pptTitleName,char(13),date]);
    
    
    %---- for single subject ----/
    opt.fontSize      =13;
    %----------------------------\
    
    % for i=1:10
    %     % ppt=ppt.addImageSlide(['added image ',num2str(i)],fullfile(pwd,'merge_power&snr_C2_m1m2.bmp'),opt);
    %     pause(1);
    % end
    
    % ppt.saveAs(fullfile(pwd,['makePPTbyMatlab.ppt']));
    
    % ppt.close;

    