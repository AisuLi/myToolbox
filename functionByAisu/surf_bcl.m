
function surf_bcl(data,colarbarAxis,xlimRange,xtick,xtickLabel)

	H1=surf(data);
	% set(H1,'LineStyle','none');
	shading interp;
	view(90,90);
	grid off;
	colorbar;
	caxis(colarbarAxis);
	xlim(xlimRange);

	if exist('xtick','var')
		set(gca,'xtick',xtick);
	end 

	if exist('xtickLabel','var')
		set(gca,'xtickLabel',xtickLabel);
	end 
