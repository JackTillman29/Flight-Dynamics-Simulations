function newSlide = AddNewSlide(openPpt,varargin)

if(nargin > 1)
    % varargin{1}: slide object of the slide after which to insert a new
    % slide
    slideObj = varargin{1};
    slideCount = slideObj.SlideNumber;
else
    % insert new slide after the LAST slide in the ppt
    slideCount = get(openPpt.Slides,'Count');
end
currentSlide = int32(double(slideCount)+1);
newSlide = invoke(openPpt.Slides,'Add',currentSlide,11);

end