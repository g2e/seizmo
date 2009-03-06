function varargout = tabpanel(varargin)
% TABPANEL Example dialog with a tab panel.
%    S = TABPANEL Returns a structure with selected options from the GUI.

% Last Modified by GUIDE v2.0 20-May-2002 13:35:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tabpanel_OpeningFcn, ...
                   'gui_OutputFcn',  @tabpanel_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    varargout{1:nargout} = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before tabpanel is made visible.
function tabpanel_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TABPANEL (see VARARGIN)

% This is the tag that all tab push buttons share.  If you have multiple
% sets of tab push buttons, each group should have unique tag.
group_name = 'tab_group';

% This is a list of the UserData values used to link tab push buttons and
% the components on their linked panels.  To add a new tab panel to the group
%  Add the button using GUIDE
%  Assign the Tag based on the group name - in this case tab_group
%  Give the UserData a unique name - e.g. another_tab_panel
%  Add components to GUIDE for the new panel
%  Give the new components the same UserData as teh tab button
%  Add the new UserData name to the below cell array
%  Modify the get_outputs function so information on the new panel
%  is returned to the calling function.
panel_names = {'orientation','paper_type'};

% tabpanelfcn('makegroups',...) adds new fields to the handles structure,
% one for each panel name and another called 'group_name_all'.  These fields
% are used by the tabpanefcn when tab_group_handler is called.
handles = tabpanelfcn('make_groups',group_name, panel_names, handles, 1);

% Update handles structure
guidata(hObject, handles);

% Since this window blocks MATLAB execution with UIWAIT, it is simplist 
% to make it modal.  Set the window style to 'normal' for debugging.
% set(handles.figure1, 'WindowStyle', 'modal')

% UIWAIT makes tabpanel wait for user response (see UIRESUME)
uiwait(handles.figure1);

% Wait to close the figure in the output function so we can access the
% updated handles structure there.


% --- Outputs from this function are returned to the command line.
function varargout = tabpanel_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% close the figure
if ishandle(handles.figure1)
    delete(handles.figure1);
end


% --- Executes on button press in tab_group.
function tab_group_Callback(hObject, eventdata, handles)
% hObject    handle to tab_group (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Call the tab_group_handler.  This updates visiblity of components as needed to
% hide the components from the previous tab and show components on this tab.
% This also updates the last_tab field in the handles structure to keep track
% of which panel was hidden.
handles = tabpanelfcn('tab_group_handler',hObject, handles, get(hObject, 'Tag'));

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in orientation_group.
function orientation_group_Callback(hObject, eventdata, handles)
% hObject    handle to orientation_group (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% deselect other radio buttons in this group
% setdiff extracts a specific entry from a vector.
set(setdiff(handles.orientation_group, hObject), ...
    'Value', 0, ...
    'Enable', 'on')

% deactivate this control
set(hObject,'Enable','inactive')


% --- Executes during object creation, after setting all properties.
function pt_papertype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pt_papertype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% On the PC, listboxes typically have a white background.  On unix they
% typically use the system background color.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Error trapping - if a callback failed, the figure will still be there but
% no longer in a UIWAIT.  Make sure the window closes nicely if that is the
% case.
if ~strcmpi(get(hObject,'waitstatus'),'waiting')
    delete(hObject);
    return;
end

% get_outputs updates the output field in the handles structure
% based on the state of various elements in the GUI
handles = get_outputs(handles, 'Cancel');

% Update handles structure
guidata(hObject, handles);

% Use uiresume instead of delete so we can get the updated handles
% structure in the output function.
uiresume(hObject);


% --- Executes on button press in ok or cancel button.
function ok_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to ok_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get_outputs updates the output field in the handles structure
% based on the state of various elements in the GUI
handles = get_outputs(handles, get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);

uiresume(handles.figure1)


% --- helper function for close, ok, and cancel callbacks
function handles = get_outputs(handles, button_name)
switch(button_name)
 case 'OK'
  orientation_output = get(findobj(handles.figure1,'Tag','orientation_group',...
                                   'Value',1), 'String');
  contents = get(handles.pt_papertype,'String');
  paper_type_output = contents{get(handles.pt_papertype,'Value')};
 case 'Cancel'
  orientation_output = [];
  paper_type_output = [];
 otherwise
  error(['Unrecognized output option: ' button_name])
end
handles.output = struct('orientation', orientation_output, ...
                        'paper_type',  paper_type_output);


% --- GUI_MAINFCN A function to handle default GUIDE GUI creation and callback dispatch.
function varargout = gui_mainfcn(gui_State, varargin)
gui_StateFields =  {'gui_Name'
                    'gui_Singleton'
                    'gui_OpeningFcn'
                    'gui_OutputFcn'
                    'gui_LayoutFcn'
                    'gui_Callback'};
gui_Mfile = '';
for i=1:length(gui_StateFields)
    if ~isfield(gui_State, gui_StateFields{i})
        error('Could not find field %s in the gui_State struct in GUI M-file %s', gui_StateFields{i}, gui_Mfile);        
    elseif isequal(gui_StateFields{i}, 'gui_Name')
        gui_Mfile = [getfield(gui_State, gui_StateFields{i}), '.m'];
    end
end

numargin = length(varargin);

if numargin == 0
    % UNTITLED
    % create the GUI
    gui_Create = 1;
elseif numargin > 3 & ischar(varargin{1}) & ishandle(varargin{2})
    % UNTITLED('CALLBACK',hObject,eventData,handles,...)
    gui_Create = 0;
else
    % UNTITLED(...)
    % create the GUI and hand varargin to the openingfcn
    gui_Create = 1;
end

if gui_Create == 0
    varargin{1} = gui_State.gui_Callback;
    if nargout
        [varargout{1:nargout}] = feval(varargin{:});
    else
        feval(varargin{:});
    end
else
    if gui_State.gui_Singleton
        gui_SingletonOpt = 'reuse';
    else
        gui_SingletonOpt = 'new';
    end
    
    % Open fig file with stored settings.  Note: This executes all component
    % specific CreateFunctions with an empty HANDLES structure.
    
    % Do feval on layout code in m-file if it exists
    if ~isempty(gui_State.gui_LayoutFcn)
        gui_hFigure = feval(gui_State.gui_LayoutFcn, gui_SingletonOpt);
    else
        gui_hFigure = local_openfig(gui_State.gui_Name, gui_SingletonOpt);            
        % If the figure has InGUIInitialization it was not completely created
        % on the last pass.  Delete this handle and try again.
        if isappdata(gui_hFigure, 'InGUIInitialization')
            delete(gui_hFigure);
            gui_hFigure = local_openfig(gui_State.gui_Name, gui_SingletonOpt);            
        end
    end
    
    % Set flag to indicate starting GUI initialization
    setappdata(gui_hFigure,'InGUIInitialization',1);

    % Fetch GUIDE Application options
    gui_Options = getappdata(gui_hFigure,'GUIDEOptions');
    
    if ~isappdata(gui_hFigure,'GUIOnScreen')
        % Adjust background color
        if gui_Options.syscolorfig 
            set(gui_hFigure,'Color', get(0,'DefaultUicontrolBackgroundColor'));
        end

        % Generate HANDLES structure and store with GUIDATA
        guidata(gui_hFigure, guihandles(gui_hFigure));
    end
    
    % If user specified 'Visible','off' in p/v pairs, don't make the figure
    % visible.
    gui_MakeVisible = 1;
    for ind=1:2:length(varargin)
        if length(varargin) == ind
            break;
        end
        len1 = min(length('visible'),length(varargin{ind}));
        len2 = min(length('off'),length(varargin{ind+1}));
        if ischar(varargin{ind}) & ischar(varargin{ind+1}) & ...
                strncmpi(varargin{ind},'visible',len1) & len2 > 1
            if strncmpi(varargin{ind+1},'off',len2)
                gui_MakeVisible = 0;
            elseif strncmpi(varargin{ind+1},'on',len2)
                gui_MakeVisible = 1;
            end
        end
    end
    
    % Check for figure param value pairs
    for index=1:2:length(varargin)
        if length(varargin) == index
            break;
        end
        try, set(gui_hFigure, varargin{index}, varargin{index+1}), catch, break, end
    end

    % If handle visibility is set to 'callback', turn it on until finished
    % with OpeningFcn
    gui_HandleVisibility = get(gui_hFigure,'HandleVisibility');
    if strcmp(gui_HandleVisibility, 'callback')
        set(gui_hFigure,'HandleVisibility', 'on');
    end
    
    feval(gui_State.gui_OpeningFcn, gui_hFigure, [], guidata(gui_hFigure), varargin{:});
    
    if ishandle(gui_hFigure)
        % Update handle visibility
        set(gui_hFigure,'HandleVisibility', gui_HandleVisibility);
        
        % Make figure visible
        if gui_MakeVisible
            set(gui_hFigure, 'Visible', 'on')
            if gui_Options.singleton 
                setappdata(gui_hFigure,'GUIOnScreen', 1);
            end
        end

        % Done with GUI initialization
        rmappdata(gui_hFigure,'InGUIInitialization');
    end
    
    % If handle visibility is set to 'callback', turn it on until finished with
    % OutputFcn
    if ishandle(gui_hFigure)
        gui_HandleVisibility = get(gui_hFigure,'HandleVisibility');
        if strcmp(gui_HandleVisibility, 'callback')
            set(gui_hFigure,'HandleVisibility', 'on');
        end
        gui_Handles = guidata(gui_hFigure);
    else
        gui_Handles = [];
    end
    
    if nargout
        [varargout{1:nargout}] = feval(gui_State.gui_OutputFcn, gui_hFigure, [], gui_Handles);
    else
        feval(gui_State.gui_OutputFcn, gui_hFigure, [], gui_Handles);
    end
    
    if ishandle(gui_hFigure)
        set(gui_hFigure,'HandleVisibility', gui_HandleVisibility);
    end
end    

function gui_hFigure = local_openfig(name, singleton)
if nargin('openfig') == 3 
    gui_hFigure = openfig(name, singleton, 'auto');
else
    % OPENFIG did not accept 3rd input argument until R13,
    % toggle default figure visible to prevent the figure
    % from showing up too soon.
    gui_OldDefaultVisible = get(0,'defaultFigureVisible');
    set(0,'defaultFigureVisible','off');
    gui_hFigure = openfig(name, singleton);
    set(0,'defaultFigureVisible',gui_OldDefaultVisible);
end
