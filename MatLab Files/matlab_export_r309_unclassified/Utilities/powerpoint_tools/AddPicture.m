function AddPicture(pptObj,filename,pos)

slidecount = this.Slides.Count;
slide = this.Slides.Item(slidecount);
image = slide.Shapes.AddPicture(filename,0,1,1,1);
set(image,'Left',pos(1),'Top',pos(2),'Width',pos(3),'Height',pos(4));

end