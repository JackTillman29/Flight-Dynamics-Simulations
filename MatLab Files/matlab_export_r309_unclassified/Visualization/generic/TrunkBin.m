function varargout = TrunkBin(varargin)
% TRUNKBIN MATLAB code for TrunkBin.fig
%      TRUNKBIN, by itself, creates a new TRUNKBIN or raises the existing
%      singleton*.
%
%      H = TRUNKBIN returns the handle to a new TRUNKBIN or the handle to
%      the existing singleton*.
%
%      TRUNKBIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRUNKBIN.M with the given input arguments.
%
%      TRUNKBIN('Property','Value',...) creates a new TRUNKBIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TrunkBin_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TrunkBin_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TrunkBin

% Last Modified by GUIDE v2.5 10-Dec-2013 11:34:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TrunkBin_OpeningFcn, ...
                   'gui_OutputFcn',  @TrunkBin_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before TrunkBin is made visible.
function TrunkBin_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TrunkBin (see VARARGIN)

% Choose default command line output for TrunkBin
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TrunkBin wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TrunkBin_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btn_loadFile.
function btn_loadFile_Callback(hObject, eventdata, handles)
% hObject    handle to btn_loadFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile('*.*','Pick a FORTRAN Real Binary File');
handles.file = [pathname filename];
[nValuesPerWrite,totalLines] = GetParameters([pathname filename]);
handles.fileprops.nValuesPerWrite = nValuesPerWrite;
handles.fileprops.totalLines      = totalLines;
set(handles.st_loadedFile,'String',handles.file);
set(handles.st_stats,'String',{ ...
    ['Values per Row: ' num2str(nValuesPerWrite)], ...
    ['Total Rows: '     num2str(totalLines)] ...
    });
set(handles.sld_preview,'Max',totalLines-1,'Min',1,'Value',1,'SliderStep',[1/totalLines .1]);
set(handles.ed_start,'String',num2str(1));
set(handles.ed_stop,'String',num2str(totalLines));
set(handles.sld_preview,'Value',1);
guidata(hObject, handles);

function [nValuesPerWrite,totalLines] = GetParameters(f)
fid = fopen(f,'rb');
nBytesPerWrite  = fread(fid,1,'uint32');
nValuesPerWrite = nBytesPerWrite / 4.0;
fseek(fid,0,'eof');
totalLines = ftell(fid)/(8+nBytesPerWrite);
frewind(fid);
fclose(fid);
  
    


% --- Executes on slider movement.
function sld_preview_Callback(hObject, eventdata, handles)
% hObject    handle to sld_preview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val = get(hObject,'Value');
set(handles.st_prevLine,'String',num2str(round(val)));
dval = round(val);

fid = fopen(handles.file,'rb');
fseek(fid,dval*(4*handles.fileprops.nValuesPerWrite+8),'bof');
fread(fid,1,'uint32');
pdata = fread(fid,handles.fileprops.nValuesPerWrite,'single');
fread(fid,1,'uint32');
fclose(fid);
set(handles.ed_preview,'String',sprintf('%18.6e %18.6e %18.6e %18.6e %18.6e',pdata(1:5)));



% --- Executes during object creation, after setting all properties.
function sld_preview_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sld_preview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function ed_preview_Callback(hObject, eventdata, handles)
% hObject    handle to ed_preview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_preview as text
%        str2double(get(hObject,'String')) returns contents of ed_preview as a double


% --- Executes during object creation, after setting all properties.
function ed_preview_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_preview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_start_Callback(hObject, eventdata, handles)
% hObject    handle to ed_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_start as text
%        str2double(get(hObject,'String')) returns contents of ed_start as a double


% --- Executes during object creation, after setting all properties.
function ed_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_stop_Callback(hObject, eventdata, handles)
% hObject    handle to ed_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_stop as text
%        str2double(get(hObject,'String')) returns contents of ed_stop as a double


% --- Executes during object creation, after setting all properties.
function ed_stop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_makeNewFile.
function btn_makeNewFile_Callback(hObject, eventdata, handles)
% hObject    handle to btn_makeNewFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
startVal = str2num(get(handles.ed_start,'String'));
stopVal  = str2num(get(handles.ed_stop, 'String'));
nVals = stopVal - startVal + 1;

fid = fopen(handles.file,'rb');
fid_out = fopen([handles.file '_' num2str(startVal) '_to_' num2str(stopVal)],'wb');

fseek(fid,startVal*(4*handles.fileprops.nValuesPerWrite+8),'bof');
 h = waitbar(0,'Truncating File');
for k = 1 : nVals
    idata1 = fread(fid,1,'uint32');
    pdata  = fread(fid,handles.fileprops.nValuesPerWrite,'single');
    idata2 = fread(fid,1,'uint32');
    
    fwrite(fid_out,idata1,'uint32');
    fwrite(fid_out,pdata ,'single');
    fwrite(fid_out,idata2,'uint32');
    waitbar(k/nVals,h);
end
fclose(fid);
fclose(fid_out);
close(h);

 
 