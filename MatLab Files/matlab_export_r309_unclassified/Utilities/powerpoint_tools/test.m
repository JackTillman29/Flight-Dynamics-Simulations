close all; clear all; clc

pptTools = 'I:\MATLAB\Utilities\powerpoint_tools\';
addpath(pptTools)

[ppt,pptCOM] = OpenPowerpointPresentation([pptTools 'test.pptx']);
ppt.SaveAs('temp.pptx'); % create a local copy of template.

CL = ppt.Slides.Item(1).CustomLayout;

% ppt.CustomLayout
ppt.Slides.AddSlide(1,CL)
% ppt.Slides.AddSlide

% ppt.Slides.Item(1)



AddPicture_norm(ppt,[pptTools 'test.bmp'],[0.2 0.2 0.5 0.5])




