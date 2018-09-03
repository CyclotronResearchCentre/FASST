function varargout = crc_z3score_report(varargin)
% CRC_Z3SCORE_REPORT MATLAB code for crc_z3score_report.fig
%      CRC_Z3SCORE_REPORT, by itself, creates a new CRC_Z3SCORE_REPORT or raises the existing
%      singleton*.
%
%      H = CRC_Z3SCORE_REPORT returns the handle to a new CRC_Z3SCORE_REPORT or the handle to
%      the existing singleton*.
%
%      CRC_Z3SCORE_REPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CRC_Z3SCORE_REPORT.M with the given input arguments.
%
%      CRC_Z3SCORE_REPORT('Property','Value',...) creates a new CRC_Z3SCORE_REPORT or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before crc_z3score_report_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to crc_z3score_report_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help crc_z3score_report

% Last Modified by GUIDE v2.5 18-Jun-2018 17:57:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @crc_z3score_report_OpeningFcn, ...
                   'gui_OutputFcn',  @crc_z3score_report_OutputFcn, ...
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

% --- Executes just before crc_z3score_report is made visible.
function crc_z3score_report_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to crc_z3score_report (see VARARGIN)

% Choose default command line output for crc_z3score_report
handles.D = varargin{1};
handles.flags = varargin{2};

start = datetime([handles.D.info.date handles.D.info.hour]);
stop = datetime([handles.D.info.date handles.D.info.hour]) + seconds(nsamples(handles.D)/fsample(handles.D));

set(handles.start_text,'String',sprintf('Recording Started at %s (HH:MM:SS)',datestr(start)));
set(handles.stop_text,'String',sprintf('Recording Stopped at %s (HH:MM:SS)',datestr(stop)));

%check if structure exits handles.D.lightsOff
if(~isfieldRecursive(handles.D,'lightsOff')),
    handles.D.lightsOff = [];
    handles.D.lightsOff = datevec(start);
else
    start = datetime(handles.D.lightsOff);
end

%check if structure exits handles.D.lightsOff
if(~isfieldRecursive(handles.D,'lightsOn')),
    handles.D.lightsOn = [];
    handles.D.lightsOn = datevec(stop);
else
    stop = datetime(handles.D.lightsOn);
end

handles.flags.Dmeg{1} = save(handles.D);

set(handles.lightsOff,'String',datestr(start,'HH:MM:SS'));
set(handles.lightsOn,'String',datestr(stop,'HH:MM:SS'));
set(handles.scorer_selecter,'String',handles.D.CRC.score(2,:));
set(handles.scorer_selecter,'Value',handles.flags.currentscore);

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);



% --- Outputs from this function are returned to the command line.
function varargout = crc_z3score_report_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;






% --- Executes on selection change in scorer_selecter.
function scorer_selecter_Callback(hObject, eventdata, handles)
% hObject    handle to scorer_selecter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns scorer_selecter contents as cell array
%        contents{get(hObject,'Value')} returns selected item from scorer_selecter


% --- Executes during object creation, after setting all properties.
function scorer_selecter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scorer_selecter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function lightsOff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lightsOff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lightsOff_Callback(hObject, eventdata, handles)
% hObject    handle to lightsOff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lightsOff as text
%        str2double(get(hObject,'String')) returns contents of lightsOff as a double



function lightsOn_Callback(hObject, eventdata, handles)
% hObject    handle to lightsOn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lightsOn as text
%        str2double(get(hObject,'String')) returns contents of lightsOn as a double


% --- Executes during object creation, after setting all properties.
function lightsOn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lightsOn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in gen_report.
function gen_report_Callback(hObject, eventdata, handles)
% hObject    handle to gen_report (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%keyboard
start = datetime([handles.D.info.date handles.D.info.hour]);
stop = datetime([handles.D.info.date handles.D.info.hour]) + seconds(nsamples(handles.D)/fsample(handles.D));
try
    %By default it is assumed that lightsOff happened the same day
    %recording was started.
    lOff = sprintf('%s-%s',datestr(start,'yyyy-mm-dd'),get(handles.lightsOff,'String'));
    lOff_time = datetime(lOff,'InputFormat','yyyy-MM-dd-HH:mm:ss');
    %If the assumption is wrong, add one day
    if lOff_time < start,
        lOff_time = lOff_time + day(1);
    end
    % Again assumption is that the lightsOn happens the same day the
    % recording is stopped
    lOn = sprintf('%s-%s',datestr(stop,'yyyy-mm-dd'),get(handles.lightsOn,'String'));
    lOn_time = datetime(lOn,'InputFormat','yyyy-MM-dd-HH:mm:ss');
catch
   msgbox('Error reading lights on/off time.','Error','error');
   return;
end

if lOff_time < start,
    msgbox('Lights off cannot happen before start of record!','Error','error');
    return;
end

if lOff_time > stop,
    msgbox('Lights off cannot happen after stop of record!','Error','error');
    return;
end

if lOn_time > stop,
    msgbox('Lights on cannot happen after stop of record!','Error','error');
    return;
end

if lOn_time > stop,
    msgbox('Lights on cannot happen after stop of record!','Error','error');
    return;
end

handles.D.lightsOff = [];
handles.D.lightsOn = [];
handles.D.lightsOff = datevec(lOff_time);
handles.D.lightsOn = datevec(lOn_time);
handles.flags.Dmeg{1} = save(handles.D);
flags.currentscore = get(handles.scorer_selecter,'Value');
guidata(hObject, handles);
crc_z3score_genreport(handles.flags.Dmeg{1}, flags);



% --- Executes when user attempts to close sleep_report.
function sleep_report_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to sleep_report (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
crc_dis_main(handles.flags);
delete(hObject);
delete(hObject);
