function [openPptObj,pptCOM] = OpenPowerpointPresentation(pptFilename)

% start powerpoint COM server, returns the handle to the COM server
% inteface
pptCOM = actxserver('PowerPoint.Application');

% make the powerpoint frame window visible
pptCOM.Visible = 1;

% open a presentation with the given filename
% openPptObj = invoke(pptHandle.Presentations,'Add',pptFilename);
% openPptObj = pptCOM.Presentations.Add(int32(1))
openPptObj = pptCOM.Presentations.Open(pptFilename);

% openPptObj.SaveAs(pptFilename);

% openPptObj.Slides.AddSlide(1,2);




end