function varargout = CoolDownGUI(varargin)
% COOLDOWNGUI MATLAB code for CoolDownGUI.fig
%      COOLDOWNGUI, by itself, creates a new COOLDOWNGUI or raises the existing
%      singleton*.
%
%      H = COOLDOWNGUI returns the handle to a new COOLDOWNGUI or the handle to
%      the existing singleton*.
%
%      COOLDOWNGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COOLDOWNGUI.M with the given input arguments.
%
%      COOLDOWNGUI('Property','Value',...) creates a new COOLDOWNGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CoolDownGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CoolDownGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CoolDownGUI

% Last Modified by GUIDE v2.5 03-Oct-2016 12:47:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CoolDownGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @CoolDownGUI_OutputFcn, ...
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

addpath('../requiredObjects/')
addpath('../requiredFunctions/')



% --- Executes just before CoolDownGUI is made visible.
function CoolDownGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CoolDownGUI (see VARARGIN)

% Choose default command line output for CoolDownGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
global isAborted
isAborted = false;

% UIWAIT makes CoolDownGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CoolDownGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
global isAborted
isAborted = false;


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
vars = evalin('base','who');
set(hObject,'String',vars)


% --- Executes on button press in updateButton.
function updateButton_Callback(hObject, eventdata, handles)
% hObject    handle to updateButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
vars = evalin('base','who');
set(handles.listbox1,'String',vars)
set(handles.text4, 'String', 'Updated Variables!')


% --- Executes on button press in runButton.
function runButton_Callback(hObject, eventdata, handles)
% hObject    handle to runButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ListBoxCheck = get(handles.listbox1,'String');
index_selected = get(handles.listbox1,'Value');
item_selected = ListBoxCheck{index_selected};

global isAborted

set(handles.text4, 'String', 'Running Algorithm')
adjMatrix = evalin('base', item_selected);

weightMatrix = CoolDownAlgorithm(adjMatrix);
assignin('base', 'weightMatrixObj', weightMatrix)
isAborted = false;


% --- Executes on button press in stopButton.
function stopButton_Callback(hObject, eventdata, handles)
% hObject    handle to stopButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global isAborted
isAborted = true;


