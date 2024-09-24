function varargout = png2mp4(varargin)

default_starting_path   = '';
default_imgdir          = '';
default_filename_method = 'consistent';
expected_filename_methods = {'consistent','inconsistent'};

p = inputParser;
p.CaseSensitive = logical(0);
passThrough = @(x) ischar(x);
p.addOptional('starting_dir',default_starting_path,passThrough)
p.addOptional('dir',default_imgdir,passThrough)
p.addOptional('method',default_filename_method, ...
    @(x) any(validatestring(x,expected_filename_methods)) );

parse(p,varargin{:});

starting_path = p.Results.starting_dir;
filename_method = p.Results.method;
imgdir = p.Results.dir;

% if(nargin > 0)
%     for k = 1:2:nargin
%         switch varargin{k}
%             case 'dir'
%                 starting_path = varargin{k+1};
%             case 'filename'
%                 filename_method = varargin{k+1};
%                 % can be ''
%         end
%     end
% else
%     starting_path = '';
% end




% assumes the filename is appended with an underscore followed by the frame
% number (e.g. "my_data_3.png" is the 3rd frame of "my_data")

% if each image has the appended frame number, but the first part of the
% filename varies, then pass in the argument

% get first file
if(isempty(imgdir) && strcmpi(filename_method,'consistent'))
    [png_file, imgdir, filterIdx] = uigetfile([starting_path '*.png']);
end

% get all files
switch lower(filename_method)
    case 'consistent'
        parts = strsplit(png_file,'_');
        for k = 1:length(parts)-1
            parts{k} = [parts{k} '_'];
        end
        png_base = strcat(parts{1:end-2});
        pngs = ls([imgdir png_base '*.png']);
    case 'inconsistent'
        % combine all imagesc in the directory
        pngs = ls([imgdir '*.png']);
end
        

npngs = length(pngs);
sorted_png = cell(1,npngs);

for k = 1:npngs
    
    this_png = strtrim(pngs(k,:));
    parts = strsplit(this_png,'_');
    frame_num = str2num(parts{end}(1:end-3));
    sorted_png{frame_num} = this_png;
    
end

vidobj = VideoWriter([imgdir 'combined.mp4'],'MPEG-4');
vidobj.FrameRate = 10;
vidobj.open;

% figure('Units','normalized','Position',[0.1888 0.0700 0.6250 0.7911]);
% hax = axes('Position',[0 0 1 1],'TickLength',[0 0]);
for k = 1:npngs
    disp(['frame ',num2str(k),' of ',num2str(npngs)])
    [A,map] = imread([imgdir sorted_png{k}]);
    
    vidobj.writeVideo(im2frame(A,map))
    
end
vidobj.close;




