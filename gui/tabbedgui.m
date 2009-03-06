function varargout = tabbedgui(varargin)
% TABBEDGUI M-file for tabbedgui.fig
%      TABBEDGUI, by itself, creates a new TABBEDGUI or raises the existing
%      singleton*.
%
%      H = TABBEDGUI returns the handle to a new TABBEDGUI or the handle to
%      the existing singleton*.
%
%      TABBEDGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TABBEDGUI.M with the given input arguments.
%
%      TABBEDGUI('Property','Value',...) creates a new TABBEDGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tabbedgui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tabbedgui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tabbedgui

% Last Modified by GUIDE v2.5 17-Feb-2008 02:18:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tabbedgui_OpeningFcn, ...
                   'gui_OutputFcn',  @tabbedgui_OutputFcn, ...
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


% --- Executes just before tabbedgui is made visible.
function tabbedgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tabbedgui (see VARARGIN)

% Choose default command line output for tabbedgui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tabbedgui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tabbedgui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in cluster.
function cluster_Callback(hObject, eventdata, handles)
% hObject    handle to cluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in correlate.
function correlate_Callback(hObject, eventdata, handles)
% hObject    handle to correlate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in window.
function window_Callback(hObject, eventdata, handles)
% hObject    handle to window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in general.
function general_Callback(hObject, eventdata, handles)
% hObject    handle to general (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in general.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to general (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pick.
function pick_Callback(hObject, eventdata, handles)
% hObject    handle to pick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in align.
function align_Callback(hObject, eventdata, handles)
% hObject    handle to align (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function general_CreateFcn(hObject, eventdata, handles)
% hObject    handle to general (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


