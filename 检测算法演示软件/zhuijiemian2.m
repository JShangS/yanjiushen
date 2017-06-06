function varargout = zhuijiemian2(varargin)
% ZHUIJIEMIAN2 MATLAB code for zhuijiemian2.fig
%      ZHUIJIEMIAN2, by itself, creates a new ZHUIJIEMIAN2 or raises the existing
%      singleton*.
%
%      H = ZHUIJIEMIAN2 returns the handle to a new ZHUIJIEMIAN2 or the handle to
%      the existing singleton*.
%
%      ZHUIJIEMIAN2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ZHUIJIEMIAN2.M with the given input arguments.
%
%      ZHUIJIEMIAN2('Property','Value',...) creates a new ZHUIJIEMIAN2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before zhuijiemian2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to zhuijiemian2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help zhuijiemian2

% Last Modified by GUIDE v2.5 10-Jul-2016 20:45:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @zhuijiemian2_OpeningFcn, ...
                   'gui_OutputFcn',  @zhuijiemian2_OutputFcn, ...
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


% --- Executes just before zhuijiemian2 is made visible.
function zhuijiemian2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to zhuijiemian2 (see VARARGIN)

% Choose default command line output for zhuijiemian2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes zhuijiemian2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = zhuijiemian2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
MTD;


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Keystone;



% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
RFT;


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
TwoRFT;
