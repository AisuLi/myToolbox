classdef PPT2007_zoe
%Matlab class creating PowerPoint 2007 Presentations
% Always make sure to have opened PowerPoint before using this class.
% 
% Example:
% ppt = PPT2007(); %create a powerpoint presentation
% ppt = ppt.addTitleSlide('10 Plots of Random Data'); %adds a title slide
%
% for i = 1:10
%   a = rand(1,10); %a is some data
%   f = figure; %create a figure
%   plot(a); %plot your data
%   xlabel('Some x axis label'); %apply some formatting and decoration
%   ylabel('Some y axis label');
%   image_file_name = fullfile(pwd, ['image_' num2str(i) '.png']); 
%   saveas(f, image_file_name); %save the figure
%   close(f);
%   %add the figure to the powerpoint
%   ppt = ppt.addImageSlide(['Figure ' num2str(i)], image_file_name);
% end
% ppt.saveAs(fullfile(pwd,'presentation.ppt')); %save the presentation.
% 

% added opt to set the font size etc.
% added a opt to fixed the ration between image highth and width
% changed by YANG
%1-dec-2011

% Created by Alan Meeson. 
% Last Edited: 3rd December 2010

%	Rev. by Zy on 2012-06-11: now can added footnote as well
%

%layout types:
% 1 - title
% 2 - title and content
% 3 - section header
% 4 - Two Content
% 5 - Comparison
% 6 - Title only
% 7 - blank
% 8 - content with caption
% 9 - picture with caption
% 10 -title and vertical text
% 11 - vertical title and text

properties
app_handle
presentation
pt = 0.0352; %point size in cm.
newline = char(13); %the new line character.
end

methods
function obj = PPT2007_zoe()
	%creates a new presentation. You must have powerpoint
	%2007 open first.
	obj.app_handle = actxserver('PowerPoint.Application'); %open powerpoint
	set(obj.app_handle,'Visible','msoTrue');

	obj.presentation = obj.app_handle.Presentation.Add; %create presentation
end

function obj = addTitleSlide(obj, title_text, sub_title_text)
	%creates a title slide
	%   ppt = ppt.addTitleSlide(title_text, sub_title_text)
	%       title_text - the text to put in the title section.
	%       sub_title_text - the text to put in the sub title
	%       section.
	%
	% If you need to use a line break (start a new line) use the
	% ppt.newline field of this object. (aka char(13)).

	%create slide
	layout = obj.presentation.SlideMaster.CustomLayouts.Item(1); %title slide layout
	Slide = obj.presentation.Slides.AddSlide(obj.presentation.Slides.Count + 1,layout);

	%do title
	if exist('title_text', 'var') && ~isempty(title_text)
		Slide.Shapes.Item(1).TextFrame.TextRange.Text = title_text;
	end

	%do sub-title
	if exist('sub_title_text', 'var') && ~isempty(sub_title_text)
		Slide.Shapes.Item(2).TextFrame.TextRange.Text = sub_title_text;
	end
end

function obj = addSectionSlide(obj, title_text, main_text)
	%creates a section title slide
	%   ppt = ppt.addTitleSlide(title_text, main_text)
	%       title_text - the text to put in the title section.
	%       main_text - the text to put in the sub title
	%       section.
	%
	% If you need to use a line break (start a new line) use the
	% ppt.newline field of this object. (aka char(13)).

	%create slide
	layout = obj.presentation.SlideMaster.CustomLayouts.Item(3); %section header
	Slide = obj.presentation.Slides.AddSlide(obj.presentation.Slides.Count + 1,layout);

	%do title
	if exist('title_text', 'var') && ~isempty(title_text)
		tb = Slide.Shapes.Item(1);
		tb.TextFrame.TextRange.Text = title_text;
	end

	%do footer
	if exist('main_text', 'var') && ~isempty(main_text)
		tb2 = Slide.Shapes.Item(2);
		tb2.TextFrame.TextRange.Text = main_text;
	end
end

%         function obj = addImageOnlySlide(obj, title_text, image_file)
%             %addImageOnlySlide - synonymous with addImageSlide. 
%             %   This function is deprecated, and only included to avoid
%             %   breaking some of my old scripts. I may remove it in future versions.
%             %   Use addImageSlide instead.
%             
%             obj = obj.addImageSlide(title_text, image_file);
%         end

function obj = addImageSlide(obj, title_text, image_file,opt,foot_Text)
	%creates a slide consisting of only a title and an image
	%   ppt = ppt.addImageOnlySlide(title_text, image_file,opt)
	%       title_text - the text to put in the title section.
	%       image_file - the filename and path of the image to use
	if ~exist('opt','var')||isempty(opt)
		opt.fontSize=18;
		opt.fontName='Calibri (Body)';
		opt.textAlignment='ppAlignLeft';
		opt.picPosition=[];
	end

	disp(title_text);
	% 	disp(foot_Text);

	nlines     =numel(strfind(title_text,char(13)))+1;
	
	fontFactor =1.5;
	
	wMaxFactor =0.95;
	hMaxFactor =0.9;
	pageSize   =get(obj.presentation.PageSetup);    

% 	disp(pageSize)
	%create slide
	layout = obj.presentation.SlideMaster.CustomLayouts.Item(6); %title only
	Slide  = obj.presentation.Slides.AddSlide(obj.presentation.Slides.Count + 1,layout);

	%do image
	if exist('image_file', 'var') && ~isempty(image_file)
		%   if image size large than the optimal size, keep it
		%   orelse, scale it to fit the optimal size
		if isfield(opt,'picPosition')&&isempty(opt.picPosition)

			%----------------------------/
			imageData  =imread(image_file);     % read image to get imageSize info
			
			imageHight =size(imageData,1);
			imageWidth =size(imageData,2);

			clear imageData					 % free memory
			%----------------------------\

			if imageHight<pageSize.SlideHeight*hMaxFactor
				targetHight=imageHight;
				sacleFactorH=1;
			else
				targetHight=pageSize.SlideHeight*hMaxFactor-nlines*opt.fontSize*fontFactor;
				sacleFactorH=targetHight/imageHight;
			end

			if imageWidth<pageSize.SlideWidth*wMaxFactor
				targetWidth=imageWidth;
				sacleFactorW=1;
			else
				targetWidth=pageSize.SlideWidth*wMaxFactor;
				sacleFactorW=targetWidth/imageWidth;
			end
			%-------- to make sure the ration is fixed -----/
			scaleFactor =min([sacleFactorH,sacleFactorW]);
			
			targetHight =imageHight*scaleFactor;
			targetWidth =imageWidth*scaleFactor;
			
			targetLeft  =(pageSize.SlideWidth-targetWidth)/2;
			
			targetTop   = nlines*opt.fontSize*fontFactor+(1-hMaxFactor)*0.5*pageSize.SlideHeight+(pageSize.SlideHeight*hMaxFactor-nlines*opt.fontSize*fontFactor-targetHight)/2;
		else
			if any(opt.picPosition~=round(opt.picPosition))
				opt.picPosition=opt.picPosition/0.0352;
			end

			targetLeft  =opt.picPosition(1);
			targetTop   =opt.picPosition(2);
			
			targetWidth =opt.picPosition(3);
			targetHight =opt.picPosition(4);
		end
		%---------------------------------------------\

		Image1 = Slide.Shapes.AddPicture(image_file,'msoFalse','msoTrue',targetLeft,targetTop,targetWidth, targetHight);
		%         disp(targetLeft,targetTop,targetWidth, targetHight);
	end

	%do title
	if exist('title_text', 'var') && ~isempty(title_text)

		Slide.Shapes.Item(1).TextFrame.TextRange.Text = title_text;

		if isfield(opt,'fontSize')&&~isempty(opt.fontSize)
			Slide.Shapes.Item(1).TextFrame.TextRange.Characters.Font.Size=opt.fontSize;
		end

		if isfield(opt,'fontName')&&~isempty(opt.fontName)
			Slide.Shapes.Item(1).TextFrame.TextRange.Characters.Font.Name=opt.fontName;
		end

		if isfield(opt,'textAlignment')&&~isempty(opt.textAlignment)
			Slide.Shapes.Item(1).TextFrame.TextRange.ParagraphFormat.Alignment =opt.textAlignment;
		end
		Slide.Shapes.Item(1).Height=opt.fontSize*fontFactor;
	end

	%do footnotes
	if exist('foot_Text', 'var') && ~isempty(foot_Text)
		footSection = Slide.Shapes.Item(1);
		tb          = Slide.Shapes.AddTextbox('msoTextOrientationHorizontal', footSection.left, footSection.top,footSection.width,footSection.height);
		
		tb.Top      = pageSize.SlideHeight*(hMaxFactor+(1-hMaxFactor)/2-0.01);
		tb.Height   = opt.fontSize*fontFactor;
		
		tb.TextFrame.TextRange.Text = foot_Text;

		if isfield(opt,'fontSize')&&~isempty(opt.fontSize)
			tb.TextFrame.TextRange.Characters.Font.Size=opt.fontSize-2;
		end

		if isfield(opt,'fontName')&&~isempty(opt.fontName)
			tb.TextFrame.TextRange.Characters.Font.Name=opt.fontName;
		end

		if isfield(opt,'textAlignment')&&~isempty(opt.textAlignment)
			tb.TextFrame.TextRange.ParagraphFormat.Alignment =opt.textAlignment;
		end
	end
end

%         function obj = addSlide(obj, title_text, image_file, caption_text)
%             %addSlide - synonymous with addCaptionedPictureSlide.
%             %   This function is deprecated, and only included to avoid
%             %   breaking some of my old scripts. I may remove it in future versions.
%             %   Use addCaptionedPictureSlide instead.
%             obj = obj.addCaptionedPictureSlide(obj, title_text, image_file, caption_text);
%         end

function obj = addCaptionedPictureSlide(obj, title_text, image_file, footer_text)
	%creates a slide of title, image, and caption
	%   ppt = ppt.addTitleSlide(title_text, image_file, footer_text)
	%       title_text - the text to put in the title section.
	%       image_file - the filename and path of the image to use.
	%       footer_text - the text to put in the footer.
	%
	% If you need to use a line break (start a new line) use the
	% ppt.newline field of this object. (aka char(13)).

	%create slide
	layout = obj.presentation.SlideMaster.CustomLayouts.Item(9); %Picture with Caption
	Slide = obj.presentation.Slides.AddSlide(obj.presentation.Slides.Count + 1,layout);

	%do title
	if exist('title_text', 'var') && ~isempty(title_text)

		Slide.Shapes.Item(1).TextFrame.TextRange.Text = title_text;


	end

	%do image
	if exist('image_file', 'var') && ~isempty(image_file)
		li = Slide.Shapes.Item(2);
		Image1 = Slide.Shapes.AddPicture(image_file,'msoFalse','msoTrue',li.left, li.top, li.width, li.height);
	end

	%do footer
	if exist('footer_text', 'var') && ~isempty(footer_text)
		Slide.Shapes.Item(3).TextFrame.TextRange.Text = footer_text;
	end
end

function obj = addTwoImageSlide(obj, title_text, left_image, right_image)
	%creates a slide of title and two images
	%   ppt = ppt.addTitleSlide(title_text, left_image, right_image)
	%       title_text - the text to put in the title section.
	%       left_image - the filename and path of the left image to use.
	%       right_image - the filename and path of the right image to use.
	%
	% If you need to use a line break (start a new line) use the
	% ppt.newline field of this object. (aka char(13)).

	%create slide
	layout = obj.presentation.SlideMaster.CustomLayouts.Item(4); %Two Content
	Slide = obj.presentation.Slides.AddSlide(obj.presentation.Slides.Count + 1,layout);

	%do title
	if exist('title_text', 'var') && ~isempty(title_text)
		Slide.Shapes.Item(1).TextFrame.TextRange.Text = title_text;
	end

	%do images
	if exist('left_image', 'var') && ~isempty(left_image)
		li = Slide.Shapes.Item(2);
		Image1 = Slide.Shapes.AddPicture(left_image,'msoFalse','msoTrue',li.left, li.top, li.width, li.height);
	else
		Slide.Shapes.Item(2).TextFrame.TextRange.Text = ' ';
	end
	if exist('right_image', 'var') && ~isempty(right_image)
		li = Slide.Shapes.Item(3);
		Image1 = Slide.Shapes.AddPicture(right_image,'msoFalse','msoTrue',li.left, li.top, li.width, li.height);
	end
end

function obj = addTwoImageComparisonSlide(obj, title_text, left_image, right_image, left_text, right_text)
	%creates a slide of title, two images and captions for each image
	%   ppt = ppt.addTwoImageComparisonSlide(title_text, left_image, right_image, left_text, right_text)
	%       title_text - the text to put in the title section.
	%       left_image - the filename and path of the left image to use.
	%       right_image - the filename and path of the right image to use.
	%       left_text - text for the left image.
	%       right_text - text for the right image.
	%
	% If you need to use a line break (start a new line) use the
	% ppt.newline field of this object. (aka char(13)).

	%create slide
	layout = obj.presentation.SlideMaster.CustomLayouts.Item(5); %Comparison
	Slide = obj.presentation.Slides.AddSlide(obj.presentation.Slides.Count + 1,layout);

	%do title
	if exist('title_text', 'var') && ~isempty(title_text)
		Slide.Shapes.Item(1).TextFrame.TextRange.Text = title_text;
	end

	%do images
	if exist('left_image', 'var') && ~isempty(left_image)
		li = Slide.Shapes.Item(3);
		Image1 = Slide.Shapes.AddPicture(left_image,'msoFalse','msoTrue',li.left, li.top, li.width, li.height);
	end
	if exist('right_image', 'var') && ~isempty(right_image)
		li = Slide.Shapes.Item(5);
		Image1 = Slide.Shapes.AddPicture(right_image,'msoFalse','msoTrue',li.left, li.top, li.width, li.height);
	end

	%do text
	if exist('left_text', 'var') && ~isempty(left_text)
		Slide.Shapes.Item(2).TextFrame.TextRange.Text = left_text;
	end
	if exist('right_text', 'var') && ~isempty(right_text)
		Slide.Shapes.Item(4).TextFrame.TextRange.Text = right_text;
	end
end

function obj = addTwoImagePlusFooterSlide(obj, title_text, left_image, right_image, footer_text)
	%creates a slide of title, two images and a footer
	%   ppt = ppt.addTwoImagePlusFooterSlide(title_text, left_image, right_image, footer_text)
	%       title_text - the text to put in the title section.
	%       left_image - the filename and path of the left image to use.
	%       right_image - the filename and path of the right image to use.
	%       footer_text - the text to put in the footer.
	%
	% If you need to use a line break (start a new line) use the
	% ppt.newline field of this object. (aka char(13)).

	%create slide
	layout = obj.presentation.SlideMaster.CustomLayouts.Item(4); %Two Content
	Slide = obj.presentation.Slides.AddSlide(obj.presentation.Slides.Count + 1,layout);

	%do title
	if exist('title_text', 'var') && ~isempty(title_text)
		Slide.Shapes.Item(1).TextFrame.TextRange.Text = title_text;
	end

	%do images
	if exist('left_image', 'var') && ~isempty(left_image)
		li = Slide.Shapes.Item(2);
		Image1 = Slide.Shapes.AddPicture(left_image,'msoFalse','msoTrue',li.left, li.top, li.width, li.height * 0.7);
	end
	if exist('right_image', 'var') && ~isempty(right_image)
		li = Slide.Shapes.Item(3);
		Image1 = Slide.Shapes.AddPicture(right_image,'msoFalse','msoTrue',li.left, li.top, li.width, li.height * 0.7);
	end

	%do footer
	if exist('footer_text', 'var') && ~isempty(footer_text)
		li1 = Slide.Shapes.Item(2);
		li2 = Slide.Shapes.Item(3);
		tb = Slide.Shapes.AddTextbox('msoTextOrientationHorizontal', li1.left, li1.top + li1.height, li1.width + li2.width, (li1.height / 0.7) * 0.3);
		tb.TextFrame.TextRange.Text = footer_text;
	end

end

function obj = saveAs(obj, filename)
	%saves the presentation to the specified file
	obj.presentation.SaveAs(filename);
end

function obj = close(obj, noQuit)
	%closes the presentation. 
	%   close(noQuit)
	%       noQuit - if false function will attempt 
	%       to close powerpoint if no presentations remain open.

	if ~exist('noQuit', 'var') || (exist('noQuit', 'var') && isempty(noQuit))
		noQuit = false;
	end

	obj.presentation.Close;
	if ~noQuit && obj.app_handle.Presentations.Count <= 0
		obj.app_handle.Quit;
	end
	obj.app_handle.delete;
end
end

end

