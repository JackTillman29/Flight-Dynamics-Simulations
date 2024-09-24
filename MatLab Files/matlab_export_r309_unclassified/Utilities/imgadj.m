function varargout = imgadj(varargin)
% IMGADJ MATLAB code for imgadj.fig
%      IMGADJ, by itself, creates a new IMGADJ or raises the existing
%      singleton*.
%
%      H = IMGADJ returns the handle to a new IMGADJ or the handle to
%      the existing singleton*.
%
%      IMGADJ('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMGADJ.M with the given input arguments.
%
%      IMGADJ('Property','Value',...) creates a new IMGADJ or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before imgadj_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to imgadj_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help imgadj

% Last Modified by GUIDE v2.5 16-Oct-2014 08:00:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @imgadj_OpeningFcn, ...
                   'gui_OutputFcn',  @imgadj_OutputFcn, ...
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


% --- Executes just before imgadj is made visible.
function imgadj_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to imgadj (see VARARGIN)

% Choose default command line output for imgadj
handles.output = hObject;
handles.thisfigure = 1;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes imgadj wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = imgadj_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function e_slider_Callback(hObject, eventdata, handles)
% hObject    handle to e_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function e_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to e_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function e_min_Callback(hObject, eventdata, handles)
% hObject    handle to e_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of e_min as text
%        str2double(get(hObject,'String')) returns contents of e_min as a double
axes_handles = findobj('Type','axes','Parent',handles.thisfigure);
caxis(axes_handles(2),[ ...
    str2num(get(handles.e_min,'String')) ...
    str2num(get(handles.e_max,'String')) ...
    ]);

%set(handles.e_max,'String',num2str(caxis_out(2)));
%set(handles.e_min,'String',num2str(caxis_out(1)));



% --- Executes during object creation, after setting all properties.
function e_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to e_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function e_max_Callback(hObject, eventdata, handles)
% hObject    handle to e_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of e_max as text
%        str2double(get(hObject,'String')) returns contents of e_max as a double
axes_handles = findobj('Type','axes','Parent',handles.thisfigure);
caxis(axes_handles(2),[ ...
    str2num(get(handles.e_min,'String')) ...
    str2num(get(handles.e_max,'String')) ...
    ]);

% --- Executes during object creation, after setting all properties.
function e_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to e_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in e_rb_linear.
function e_rb_linear_Callback(hObject, eventdata, handles)
% hObject    handle to e_rb_linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of e_rb_linear


% --- Executes on button press in e_rb_log.
function e_rb_log_Callback(hObject, eventdata, handles)
% hObject    handle to e_rb_log (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of e_rb_log


% --- Executes on button press in e_btn_gcf.
function e_btn_gcf_Callback(hObject, eventdata, handles)
% hObject    handle to e_btn_gcf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected object is changed in e_scale_change.
function e_scale_change_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in e_scale_change 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected object is changed in e_colormap_change.
function e_colormap_change_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in e_colormap_change 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

axes_handles = findobj('Type','axes','Parent',handles.thisfigure);

if(get(handles.e_rb_c_jet,'Value'))
    colormap(axes_handles(1),jet);
elseif(get(handles.e_rb_c_bone,'Value'))
    colormap(axes_handles(1),bone);
elseif(get(handles.e_rb_c_cool,'Value'))
    colormap(axes_handles(1),cool);
elseif(get(handles.e_rb_c_gray,'Value'))
    colormap(axes_handles(1),gray);
elseif(get(handles.e_rb_c_hot,'Value'))
    colormap(axes_handles(1),hot);
elseif(get(handles.e_rb_c_hsv,'Value'))
    colormap(axes_handles(1),hsv);
elseif(get(handles.e_rb_c_green,'Value'))
    qcmap = [0*(0:255)' (0:255)' 0*(0:255)'];  % Green
    qcmap = qcmap ./ 255;
    colormap(axes_handles(1),qcmap);    
end



function e_fignum_Callback(hObject, eventdata, handles)
% hObject    handle to e_fignum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of e_fignum as text
%        str2double(get(hObject,'String')) returns contents of e_fignum as a double
handles.thisfigure = str2num(get(hObject,'String'));
axes_handles = findobj('Type','axes','Parent',handles.thisfigure);
caxis_out = caxis(axes_handles(2));

set(handles.e_max,'String',num2str(caxis_out(2)));
set(handles.e_min,'String',num2str(caxis_out(1)));

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function e_fignum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to e_fignum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
