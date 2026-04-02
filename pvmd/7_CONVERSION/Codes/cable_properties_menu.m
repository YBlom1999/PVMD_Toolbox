function varargout = cable_properties_menu(type,varargin)
% Function related to the cable_properties_menu.fig
% Touch if you know what you're doing (I don't)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', [], ...
                   'gui_OutputFcn',  @STRING_CABLE_MENU_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if (nargin>1) && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if strcmp(type,'string')
    gui_State.gui_OpeningFcn = @STRING_CABLE_MENU_OpeningFcn;
elseif strcmp(type,'inverter')
    gui_State.gui_OpeningFcn = @INV_CABLE_MENU_OpeningFcn;
else
    disp('No type of cable selected')
end 

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

function STRING_CABLE_MENU_OpeningFcn(hObject, ~, handles, varargin)
handles.output = hObject;
handles.figure1.Name = 'String cable properties';
guidata(hObject, handles);

function INV_CABLE_MENU_OpeningFcn(hObject, ~, handles, varargin)
handles.output = hObject;
handles.figure1.Name = 'Inverter cable properties';
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = STRING_CABLE_MENU_OutputFcn(~, ~, handles) 
varargout{1} = handles.output;
uiwait

% --- Executes during object creation, after setting all properties.
function STRING_material_menu_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function STRING_crosssec_menu_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function STRING_length_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on selection change in STRING_material_menu.
function STRING_material_menu_Callback(hObject, ~, ~)
STRING_material = get(hObject,'Value');
assignin('base','STRING_material',STRING_material);

% --- Executes on selection change in STRING_crosssec_menu.
function STRING_crosssec_menu_Callback(hObject, ~, ~)
STRING_cross_sec = get(hObject,'Value');
assignin('base','STRING_cross_sec',STRING_cross_sec);

function STRING_length_Callback(hObject, ~, ~)
STRING_length = get(hObject,'string');
assignin('base','STRING_length',STRING_length);

% --- Executes on button press in confirm.
function confirm_Callback(~, ~, ~)
uiresume
closereq();

% --- Executes on button press in cancel.
function cancel_Callback(~, ~, ~)
closereq();

% --- Executes on button press in HELP.
function pushbutton4_Callback(~, ~, ~)
msgbox({'Please refer to the Use