function varargout = RFT(varargin)
% RFT MATLAB code for RFT.fig
%      RFT, by itself, creates a new RFT or raises the existing
%      singleton*.
%
%      H = RFT returns the handle to a new RFT or the handle to
%      the existing singleton*.
%
%      RFT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RFT.M with the given input arguments.
%
%      RFT('Property','Value',...) creates a new RFT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RFT_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RFT_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% tishi_edit the above text to modify the response to help RFT

% Last Modified by GUIDE v2.5 29-Jun-2016 18:46:29

% Begin initialization code - DO NOT TISHI_EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RFT_OpeningFcn, ...
                   'gui_OutputFcn',  @RFT_OutputFcn, ...
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
% End initialization code - DO NOT TISHI_EDIT


% --- Executes just before RFT is made visible.
function RFT_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RFT (see VARARGIN)
set(handles.fc_edit,'String',1);
set(handles.B_edit,'String',1);
set(handles.Tp_edit,'String',128);
set(handles.Fs_edit,'String',1);
set(handles.PRF_edit,'String',1000);
set(handles.SNR_edit,'String',10);
set(handles.TCI_edit,'String',0.2);
set(handles.T1_edit,'String',0);
set(handles.T2_edit,'String',0);
set(handles.T3_edit,'String',1);
set(handles.V1_edit,'String',1);
set(handles.V2_edit,'String',2);
set(handles.V3_edit,'String',10);
set(handles.a1_edit,'String',0);
set(handles.a2_edit,'String',0);
set(handles.a3_edit,'String',0);
set(handles.R1_edit,'String',100);
set(handles.R2_edit,'String',20);
set(handles.R3_edit,'String',60);
set(handles.P_edit,'String',20);
set(handles.RFTcost_edit,'String',0);
set(handles.Tishi_edit,'String',['点击算法运行开始运算']);
C=3e8;
fc=str2num(get(handles.fc_edit,'String'))*1e9;
P=str2num(get(handles.P_edit,'String'));
PRF=str2num(get(handles.PRF_edit,'String'));
lamda=C/fc;
Tr=1/PRF
Vb=lamda/(2*Tr);
set(handles.Vb_edit,'String',num2str(Vb));
Fs=str2num(get(handles.Fs_edit,'String'))*1e6;
B=str2num(get(handles.B_edit,'String'))*1e6;
QuJian1=-P*Vb/round(Fs/B);
QuJian2=P*Vb/round(Fs/B);
set(handles.QuJian_edit,'String',[num2str(QuJian1),',',num2str(QuJian2)]);
% Choose default command line output for RFT
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes RFT wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RFT_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function fc_edit_Callback(hObject, eventdata, handles)
% hObject    handle to fc_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fc_edit as text
%        str2double(get(hObject,'String')) returns contents of fc_edit as a double


% --- Executes during object creation, after setting all properties.
function fc_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fc_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: Tishi_edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function B_edit_Callback(hObject, eventdata, handles)
% hObject    handle to B_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of B_edit as text
%        str2double(get(hObject,'String')) returns contents of B_edit as a double


% --- Executes during object creation, after setting all properties.
function B_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to B_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: Tishi_edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tp_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Tp_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tp_edit as text
%        str2double(get(hObject,'String')) returns contents of Tp_edit as a double


% --- Executes during object creation, after setting all properties.
function Tp_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tp_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: Tishi_edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Fs_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Fs_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Fs_edit as text
%        str2double(get(hObject,'String')) returns contents of Fs_edit as a double


% --- Executes during object creation, after setting all properties.
function Fs_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Fs_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: Tishi_edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PRF_edit_Callback(hObject, eventdata, handles)
% hObject    handle to PRF_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PRF_edit as text
%        str2double(get(hObject,'String')) returns contents of PRF_edit as a double


% --- Executes during object creation, after setting all properties.
function PRF_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PRF_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: Tishi_edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SNR_edit_Callback(hObject, eventdata, handles)
% hObject    handle to SNR_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SNR_edit as text
%        str2double(get(hObject,'String')) returns contents of SNR_edit as a double


% --- Executes during object creation, after setting all properties.
function SNR_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SNR_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: Tishi_edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TCI_edit_Callback(hObject, eventdata, handles)
% hObject    handle to TCI_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TCI_edit as text
%        str2double(get(hObject,'String')) returns contents of TCI_edit as a double


% --- Executes during object creation, after setting all properties.
function TCI_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TCI_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: Tishi_edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NTCI_edit_Callback(hObject, eventdata, handles)
% hObject    handle to NTCI_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NTCI_edit as text
%        str2double(get(hObject,'String')) returns contents of NTCI_edit as a double


% --- Executes during object creation, after setting all properties.
function NTCI_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NTCI_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: Tishi_edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function V1_edit_Callback(hObject, eventdata, handles)
% hObject    handle to V1_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of V1_edit as text
%        str2double(get(hObject,'String')) returns contents of V1_edit as a double


% --- Executes during object creation, after setting all properties.
function V1_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to V1_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: Tishi_edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function T1_edit_Callback(hObject, eventdata, handles)
% hObject    handle to T1_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T1_edit as text
%        str2double(get(hObject,'String')) returns contents of T1_edit as a double


% --- Executes during object creation, after setting all properties.
function T1_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T1_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: Tishi_edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function R1_edit_Callback(hObject, eventdata, handles)
% hObject    handle to R1_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R1_edit as text
%        str2double(get(hObject,'String')) returns contents of R1_edit as a double


% --- Executes during object creation, after setting all properties.
function R1_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R1_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: Tishi_edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function V2_edit_Callback(hObject, eventdata, handles)
% hObject    handle to V2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of V2_edit as text
%        str2double(get(hObject,'String')) returns contents of V2_edit as a double


% --- Executes during object creation, after setting all properties.
function V2_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to V2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: Tishi_edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function T2_edit_Callback(hObject, eventdata, handles)
% hObject    handle to T2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T2_edit as text
%        str2double(get(hObject,'String')) returns contents of T2_edit as a double


% --- Executes during object creation, after setting all properties.
function T2_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: Tishi_edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function R2_edit_Callback(hObject, eventdata, handles)
% hObject    handle to R2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R2_edit as text
%        str2double(get(hObject,'String')) returns contents of R2_edit as a double


% --- Executes during object creation, after setting all properties.
function R2_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: Tishi_edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function V3_edit_Callback(hObject, eventdata, handles)
% hObject    handle to V3_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of V3_edit as text
%        str2double(get(hObject,'String')) returns contents of V3_edit as a double


% --- Executes during object creation, after setting all properties.
function V3_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to V3_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: Tishi_edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function T3_edit_Callback(hObject, eventdata, handles)
% hObject    handle to T3_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T3_edit as text
%        str2double(get(hObject,'String')) returns contents of T3_edit as a double


% --- Executes during object creation, after setting all properties.
function T3_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T3_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: Tishi_edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function R3_edit_Callback(hObject, eventdata, handles)
% hObject    handle to R3_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R3_edit as text
%        str2double(get(hObject,'String')) returns contents of R3_edit as a double


% --- Executes during object creation, after setting all properties.
function R3_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R3_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: Tishi_edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function P_edit_Callback(hObject, eventdata, handles)
% hObject    handle to P_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of P_edit as text
%        str2double(get(hObject,'String')) returns contents of P_edit as a double


% --- Executes during object creation, after setting all properties.
function P_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to P_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: Tishi_edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MTD_radiobutton.
function MTD_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to MTD_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MTD_radiobutton


% --- Executes on button press in Keystone_radiobutton.
function Keystone_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to Keystone_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Keystone_radiobutton


% --- Executes on button press in RFT_radiobutton.
function RFT_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to RFT_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RFT_radiobutton



function MTDcost_edit_Callback(hObject, eventdata, handles)
% hObject    handle to MTDcost_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MTDcost_edit as text
%        str2double(get(hObject,'String')) returns contents of MTDcost_edit as a double


% --- Executes during object creation, after setting all properties.
function MTDcost_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MTDcost_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: Tishi_edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Run_pushbutton.
function Run_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Run_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.Tishi_edit,'String','正在运行请稍等....');
pause(0.0001);
fc=str2num(get(handles.fc_edit,'String'))*1e9;
B=str2num(get(handles.B_edit,'String'))*1e6;
Tao=str2num(get(handles.Tp_edit,'String'))*1e-6;
Fs=str2num(get(handles.Fs_edit,'String'))*1e6;
PRF=str2num(get(handles.PRF_edit,'String'));
SNR=str2num(get(handles.SNR_edit,'String'));
Tci=str2num(get(handles.TCI_edit,'String'));
Tr=1/PRF;
M=round(Tci*PRF);
A1=str2num(get(handles.T1_edit,'String'));
A2=str2num(get(handles.T2_edit,'String'));
A3=str2num(get(handles.T3_edit,'String'));
V1=str2num(get(handles.V1_edit,'String'))*340;
V2=str2num(get(handles.V2_edit,'String'))*340;
V3=str2num(get(handles.V3_edit,'String'))*340;
a1=str2num(get(handles.a1_edit,'String'));
a2=str2num(get(handles.a2_edit,'String'));
a3=str2num(get(handles.a3_edit,'String'));
C=3e8;
delt_R=C/2/Fs;
R1=round(Tao/Fs)-str2num(get(handles.R1_edit,'String'))*delt_R;
R2=round(Tao/Fs)-str2num(get(handles.R2_edit,'String'))*delt_R;
R3=round(Tao/Fs)-str2num(get(handles.R3_edit,'String'))*delt_R;
P=str2num(get(handles.P_edit,'String'));
t=-Tao/2:1/Fs:Tao/2-1/Fs;%脉冲时间
L=length(t);
mu=B/Tao;%条频率
lamda=C/fc;
Vb=lamda/(2*Tr);
set(handles.Vb_edit,'String',num2str(Vb));
pause(0.00000001);
QuJian1=-P*Vb/round(Fs/B);
QuJian2=P*Vb/round(Fs/B);
set(handles.QuJian_edit,'String',[num2str(QuJian1),',',num2str(QuJian2)]);
pause(0.00000001);
for i=1:M
    V1=V1+a1*(i-1)*Tr;
    fd1(i)=2*V1/lamda;
    delt_t1(i)=2*(R1+V1*Tr*(i-1))/C;%回拨延迟
    
    V2=V2+a2*(i-1)*Tr;
    fd2(i)=2*V2/lamda;
    delt_t2(i)=2*(R2+V2*Tr*(i-1))/C;%回拨延迟
    
    V3=V3+a3*(i-1)*Tr;
    fd3(i)=2*V3/lamda;
    delt_t3(i)=2*(R3+V3*Tr*(i-1))/C;%回拨延迟
end
%%
%%脉压系数
    ht_t=exp(-1j*2*pi*(mu/2*(t).^2));
    ht=conj((ht_t));
    ht_fft=(fft(ht));%fftshift
%%
for i=1:M
    echo1(i,:)=A1*exp(-1j*2*pi*(mu/2*(t+delt_t1(i)).^2)+-1j*2*pi*(fc)*(delt_t1(i)));
    echo2(i,:)=A2*exp(-1j*2*pi*(mu/2*(t+delt_t2(i)).^2)+-1j*2*pi*(fc)*(delt_t2(i)));
    echo3(i,:)=A3*exp(-1j*2*pi*(mu/2*(t+delt_t3(i)).^2)+-1j*2*pi*(fc)*(delt_t3(i)));
    echo(i,:)=echo1(i,:)+echo2(i,:)+echo3(i,:);
    echo(i,:)=awgn(echo(i,:),SNR);%%加噪声
    echo_fft(i,:)=(fft(echo(i,:)));%fftshift
    pc_result(i,:)=(ifft(echo_fft(i,:).*ht_fft));%ifftshift
    pc_result_fft(i,:)=(fft(pc_result(i,:)));%fftshift%快时间FFT
end
axes(handles.PC_figure);
imagesc ((abs(pc_result)))
title('脉压结果','FontSize',18)
xlabel('距离单元')
ylabel('脉冲单元')
axes(handles.PCmove_figure);
plot(abs(pc_result(1,:))/max(max(abs(pc_result(1,:)))),'r')
hold on
plot(abs(pc_result(M,:))/max(max(abs(pc_result(M,:)))),'b')
hold off
title('脉冲走动示意图','FontSize',18)
xlabel('距离单元')
ylabel('归一化幅度')
di1=['第1个脉冲'];
diM=['第',num2str(M),'个脉冲'];
legend(di1,diM);

%%%CZT_RFT%%%%%%%%%%%%%%%%%
fai=lamda*B/(L*C);
delta_V=lamda/(2*M*Tr);
%%系数生成
% V=Vb-delta_V:-delta_V:0;%%速度搜索
num_sou=M;
Sr=zeros(num_sou,L);
m=1:M;
m=m.';
a=1-fai.*(1:L);
Cn=exp(-1j*2*pi*a./num_sou);
Cn=Cn.';
kk = ( (-M+1):M ).';
kk2 = (kk .^ 2) ./ 2;
for i=1:length(Cn)
    ww(:,i)= Cn(i).^ (kk2); 
end
v=1 ./ ww;
fv = fft(v);
nfft = 2*M;%2^nextpow2(M+M-1);
L_i=(0:L-1);
Bu_mL=m*L_i;
[~,n]=size(pc_result_fft);
K=-P:P;
sp_final=zeros(M,L);
max_sp=0;
max_K=0;
    tic
    for Nk=1:length(K)
        set(handles.Tishi_edit,'String',[num2str(Nk/(2*P+1)*100),'%']);
        pause(0.0000001);
        %%乘以补偿因子补偿多普勒模糊
        BU=exp(-1j*2*pi*K(Nk).*Bu_mL.*lamda*B./(L*C));%%补偿因子
        x= pc_result_fft.*BU;
        nn = (0:(M-1))';
        y = x .* ww(M+nn,:);
        fy = fft( y, nfft);
        fy = fy .* fv;
        g  = ifft(fy);
        sp=g(M:2*M-1,:).* ww(M:2*M-1,:);%(M:2*M-1,:)
        sp=(ifft(sp,[],2));
      if (A1+A2+A3>1)
            bijiao_t=sp_final>sp;
            sp_final=sp_final.*bijiao_t+sp.*(ones(M,L)-bijiao_t);
      else
           max_t=max(max(abs(sp)));
%         比记录的最大峰值大则替换
            if max_t>max_sp
               max_sp=max_t;
               sp_final=sp;
               max_K=K(Nk);
            end
      end
    end
    time=toc;
    set(handles.RFTcost_edit,'String',num2str(time));
    set(handles.Tishi_edit,'String','运行结束');
    set(handles.Tishi_edit,'String','改变搜索模糊数可改变速度搜索范围');
    axes(handles.RFT_figure);
    V=(max_K-1)*Vb:delta_V:max_K*Vb-delta_V;
    [X,Y]=meshgrid((1:L),V);
    mesh((abs(sp_final))/max(max(abs(sp_final))));
    str_RFT=['SNR=',num2str(SNR),'时的RFT结果'];
    title(str_RFT,'FontSize',18)
    xlabel('距离单元')
    ylabel('速度单元')
    zlabel('归一化幅度')

%     V_all=Vb*(K(1)-1):delta_V:Vb*(K(end))-delta_V;
%     [Y_ALL,X_ALL]=meshgrid((1:L),V_all);
    
%     mesh(Y_ALL,X_ALL,(abs(sp_ALL)))


    
    


function Tishi_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Tishi_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tishi_edit as text
%        str2double(get(hObject,'String')) returns contents of Tishi_edit as a double


% --- Executes during object creation, after setting all properties.
function Tishi_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tishi_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: Tishi_edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  set(handles.Tishi_edit,'String','程序已停止，点击算法运行重新开始!');
  uiwait


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% set(gcf,'visible','off');
close(gcf)



function Vb_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Vb_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Vb_edit as text
%        str2double(get(hObject,'String')) returns contents of Vb_edit as a double


% --- Executes during object creation, after setting all properties.
function Vb_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Vb_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function QuJian_edit_Callback(hObject, eventdata, handles)
% hObject    handle to QuJian_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of QuJian_edit as text
%        str2double(get(hObject,'String')) returns contents of QuJian_edit as a double


% --- Executes during object creation, after setting all properties.
function QuJian_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to QuJian_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function a1_edit_Callback(hObject, eventdata, handles)
% hObject    handle to a1_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a1_edit as text
%        str2double(get(hObject,'String')) returns contents of a1_edit as a double


% --- Executes during object creation, after setting all properties.
function a1_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a1_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function a2_edit_Callback(hObject, eventdata, handles)
% hObject    handle to a2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a2_edit as text
%        str2double(get(hObject,'String')) returns contents of a2_edit as a double


% --- Executes during object creation, after setting all properties.
function a2_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function a3_edit_Callback(hObject, eventdata, handles)
% hObject    handle to a3_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a3_edit as text
%        str2double(get(hObject,'String')) returns contents of a3_edit as a double


% --- Executes during object creation, after setting all properties.
function a3_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a3_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function amin_edit_Callback(hObject, eventdata, handles)
% hObject    handle to amin_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of amin_edit as text
%        str2double(get(hObject,'String')) returns contents of amin_edit as a double


% --- Executes during object creation, after setting all properties.
function amin_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to amin_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function amax_edit_Callback(hObject, eventdata, handles)
% hObject    handle to amax_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of amax_edit as text
%        str2double(get(hObject,'String')) returns contents of amax_edit as a double


% --- Executes during object creation, after setting all properties.
function amax_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to amax_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Keystonecost_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Keystonecost_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Keystonecost_edit as text
%        str2double(get(hObject,'String')) returns contents of Keystonecost_edit as a double


% --- Executes during object creation, after setting all properties.
function Keystonecost_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Keystonecost_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function RFTcost_edit_Callback(hObject, eventdata, handles)
% hObject    handle to RFTcost_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RFTcost_edit as text
%        str2double(get(hObject,'String')) returns contents of RFTcost_edit as a double


% --- Executes during object creation, after setting all properties.
function RFTcost_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RFTcost_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TwoRFTcost_edit_Callback(hObject, eventdata, handles)
% hObject    handle to TwoRFTcost_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TwoRFTcost_edit as text
%        str2double(get(hObject,'String')) returns contents of TwoRFTcost_edit as a double


% --- Executes during object creation, after setting all properties.
function TwoRFTcost_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TwoRFTcost_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
