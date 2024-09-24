function AddPicture_norm(pptObj,filename,pos)

slidecount = pptObj.Slides.Count;
slide = pptObj.Slides.Item(slidecount);

slideHeight = pptObj.PageSetup.SlideHeight;
slideWidth  = pptObj.PageSetup.SlideWidth;

pos(1) = pos(1) * slideWidth;
pos(2) = pos(2) * slideHeight;
pos(3) = pos(3) * slideWidth;
pos(4) = pos(4) * slideHeight;

image = slide.Shapes.AddPicture(filename,0,1,1,1);
set(image,'Left',pos(1),'Top',pos(2),'Width',pos(3),'Height',pos(4));

end