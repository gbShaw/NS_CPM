function varargout = NS_CPM_confirm(varargin)
% NS_CPM_CONFIRM MATLAB code for NS_CPM_confirm.fig
%      NS_CPM_CONFIRM, by itself, creates a new NS_CPM_CONFIRM or raises the existing
%      singleton*.
%
%      H = NS_CPM_CONFIRM returns the handle to a new NS_CPM_CONFIRM or the handle to
%      the existing singleton*.
%
%      NS_CPM_CONFIRM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NS_CPM_CONFIRM.M with the given input arguments.
%
%      NS_CPM_CONFIRM('Property','Value',...) creates a new NS_CPM_CONFIRM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NS_CPM_confirm_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NS_CPM_confirm_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NS_CPM_confirm

% Last Modified by GUIDE v2.5 04-Nov-2021 22:12:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NS_CPM_confirm_OpeningFcn, ...
                   'gui_OutputFcn',  @NS_CPM_confirm_OutputFcn, ...
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

%  guidata(hObject, handles);
 




% --- Executes just before NS_CPM_confirm is made visible.
function NS_CPM_confirm_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NS_CPM_confirm (see VARARGIN)

% Choose default command line output for NS_CPM_confirm
handles.output = hObject;

% Update handles structure

handles.Data = varargin{1};
handles.CPMcfg = varargin{2};
confirm = handles.Data.confirm;
set(handles.uitable1,'data',confirm)
 
set(handles.uitable1,'columnname',[{'FC filenames' 'scores' } ...
             repmat({'covariable'},1,size(confirm,2)-2 )  ])

%  handles.figure1.Position(1:2) = handles.CPMcfg.mainFigurePosition(1:2);   
%  figobj =
 position =  get(handles.figure1,'position');
 position(1:2) = handles.CPMcfg.mainFigurePosition(1:2);
 set(figobj,'Position',position);   
%  set(figboj,'Position',position);
 guidata(hObject, handles);



% UIWAIT makes NS_CPM_confirm wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = NS_CPM_confirm_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles

set(handles.figure1,'visible','off');
drawnow
 [ALLPrediction, ALLPerformance ] = run_NS_CPM(handles.Data,handles.CPMcfg);
 handles.ALLPerformance = ALLPerformance;
  guidata(hObject, handles);
assignin('base', 'ALLPerformance',ALLPerformance);
 delete(handles.figure1)


% --- Executes during object creation, after setting all properties.
function uitable1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
