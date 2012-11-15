function varargout = crc_SP_detection(varargin)
% CRC_SP_DETECTION M-file for crc_SP_detection.fig
%
% GUI to select SWS detection parameters
%
%      CRC_SP_DETECTION, by itself, creates a new CRC_SP_DETECTION or raises the existing
%      singleton*.
%
%      H = CRC_SP_DETECTION returns the handle to a new CRC_SP_DETECTION or the handle to
%      the existing singleton*.
%
%      CRC_SP_DETECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CRC_SP_DETECTION.M with the given input arguments.
%
%      CRC_SP_DETECTION('Property','Value',...) creates a new CRC_SP_DETECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before crc_SP_detection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to crc_SP_detection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%__________________________________________________________________________
% Copyright (C) 2010 Cyclotron Research Centre
%
% Written by J. Schrouff, 2010
% Cyclotron Research Centre, University of Liege, Belgium
% $Id:$

% Last Modified by GUIDE v2.5 01-Apr-2011 15:29:55


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @crc_SP_detection_OpeningFcn, ...
                   'gui_OutputFcn',  @crc_SP_detection_OutputFcn, ...
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


% --- Executes just before crc_SP_detection is made visible.
function crc_SP_detection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to crc_SP_detection (see VARARGIN)

% Choose default command line output for crc_SP_detection
handles.output = hObject;
crcdef = crc_get_defaults('sp');

%choose defaults
handles.highfc=crcdef.highfc;
handles.lowfc=crcdef.lowfc;
handles.review=0;
handles.fname=[];
handles.reref=0;
handles.scorer=1;
handles.analyse=1;
handles.Begpts=0;
handles.Endpts=1;
set(handles.edit_highpassfc,'String',num2str(handles.highfc))
set(handles.edit_lowpassfc,'String',num2str(handles.lowfc))
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes crc_SP_detection wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = crc_SP_detection_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function EEGfilename_Callback(hObject, eventdata, handles)
% hObject    handle to EEGfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EEGfilename as text
%        str2double(get(hObject,'String')) returns contents of EEGfilename as a double
handles.fname=get(hObject,'String');
% Update handles structure
guidata(hObject, handles);
try
    D = crc_eeg_load(handles.fname);
catch
    error('wrong file format or name')
end
if isfield(D,'CRC') && isfield(D.CRC,'score')
    sc=[];
    for i=1:size(D.CRC.score,2)
        sc=[sc; {D.CRC.score{2,i}}];
    end
    set(handles.scoringlist,'String',sc)
end
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function EEGfilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EEGfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browsefilename.
function browsefilename_Callback(hObject, eventdata, handles)
% hObject    handle to browsefilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.fname=spm_select(1, 'mat', 'Select EEG file');
set(handles.EEGfilename, 'string',handles.fname)
% Update handles structure
guidata(hObject, handles);
try
    D = crc_eeg_load(handles.fname);
catch
    error('wrong file format or name')
end
if isfield(D,'CRC') && isfield(D.CRC,'score')
    sc=[];
    for i=1:size(D.CRC.score,2)
        sc=[sc; {D.CRC.score{2,i}}];
    end
    set(handles.scoringlist,'String',sc)
end
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in radiobutton_rereferenced.
function radiobutton_rereferenced_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_rereferenced (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_rereferenced
handles.reref=get(hObject,'Value');
% Update handles structure
guidata(hObject, handles);

function edit_highpassfc_Callback(hObject, eventdata, handles)
% hObject    handle to edit_highpassfc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_highpassfc as text
%        str2double(get(hObject,'String')) returns contents of edit_highpassfc as a double
handles.highfc=str2double(get(hObject,'String'));
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_highpassfc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_highpassfc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_lowpassfc_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lowpassfc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lowpassfc as text
%        str2double(get(hObject,'String')) returns contents of edit_lowpassfc as a double
handles.lowfc=str2double(get(hObject,'String'));
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_lowpassfc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lowpassfc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton_review.
function radiobutton_review_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_review (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_review
handles.review=get(hObject,'Value');
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbutton_OK.
function pushbutton_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
crc_SP_detect(handles)


% --- Executes on button press in rb_score.
function rb_score_Callback(hObject, eventdata, handles)
% hObject    handle to rb_score (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_score


% --- Executes on button press in rd_wav.
function rd_wav_Callback(hObject, eventdata, handles)
% hObject    handle to rd_wav (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.wav=get(hObject,'Value');
% Hint: get(hObject,'Value') returns toggle state of rd_wav
guidata(hObject, handles);


% --- Executes on button press in cb_an_ft.
function cb_an_ft_Callback(hObject, eventdata, handles)
% hObject    handle to cb_an_ft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_an_ft
ft=get(hObject,'Value');
if ft==1
    set(handles.cb_an_all, 'Enable','off')
    set(handles.cb_an_3, 'Enable','off')
    set(handles.edit_scores, 'Enable','off')
    handles.analyse=1;
else
    set(handles.cb_an_all, 'Enable','on')
    set(handles.cb_an_3, 'Enable','on')
    set(handles.edit_scores, 'Enable','on')
end
% Update handles structure
guidata(hObject, handles);


function edit_starttime_Callback(hObject, eventdata, handles)
% hObject    handle to edit_starttime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_starttime as text
%        str2double(get(hObject,'String')) returns contents of edit_starttime as a double
set(handles.edit_starttime,'ForegroundColor','black')
sttime=get(hObject,'String');
if size(sttime,2)==8
    handles.start=[str2double(sttime(1:2)), str2double(sttime(4:5)),...
        str2double(sttime(7:8))];
    set(handles.edit_starttime,'ForegroundColor','blue')
elseif size(sttime,2)==6
    handles.start=[str2double(sttime(1:2)), str2double(sttime(3:4)),...
        str2double(sttime(5:6))];
    set(handles.edit_starttime,'ForegroundColor','blue')
else
    beep
    disp(['Error: enter time in the hh-mm-ss format'])
    set(handles.edit_starttime,'ForegroundColor','red')
    return
end

if ~isempty(handles.fname)
    D = crc_eeg_load(handles.fname);
    lgt = (nsamples(D)/fsample(D))/60^2;
    
    if isfield(D, 'info')
        if isfield(D.info, 'hour')
            hh=D.info.hour;
        else
            D.info=struct('date',[0 0 0],'hour',[0 0 0]);
            hh=D.info.hour;
            disp('Warning: time must be relative to recording start!')
        end
    else
        D.info = [];
        D.info = struct('date',[0 0 0],'hour',[0 0 0]);
        hh = D.info.hour;
        disp('Warning: time must be relative to recording start!')
    end
    
    save(D);
    datenumstart=datenum([0 0 0 hh]);
    datenumEnd=datenum([0 0 0 handles.start(1) handles.start(2) handles.start(3)]);

    datediff=datevec(datenumEnd-datenumstart);
    if datediff(4)>lgt
        datenumEnd=datenum([0 0 1 handles.start(1) handles.start(2) handles.start(3)]);
        datediff=datevec(datenumEnd-datenumstart);
    end
    Begpts=min(nsamples(D),max(1,round((datediff(4)*60^2+datediff(5)*60+datediff(6))*fsample(D))));
    if Begpts==nsamples(D)
        beep
        disp('Error: start time is not coherent with the file timing -- correct please')
        disp('Enter time in the hh-mm-ss format')
        set(handles.edit_starttime,'ForegroundColor','red')
    end
else
    beep
    disp('Select EEG file first!')
    set(handles.EEGfilename,'ForegroundColor','red')
end
handles.Begpts=Begpts;
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_starttime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_starttime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_stoptime_Callback(hObject, eventdata, handles)
% hObject    handle to edit_stoptime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_stoptime as text
%        str2double(get(hObject,'String')) returns contents of edit_stoptime as a double
set(handles.edit_stoptime,'ForegroundColor','black')
sptime=get(hObject,'String');
if size(sptime,2)==8
    handles.stop=[str2double(sptime(1:2)), str2double(sptime(4:5)),...
        str2double(sptime(7:8))];
    set(handles.edit_stoptime,'ForegroundColor','blue')
elseif size(sptime,2)==6
    handles.stop=[str2double(sptime(1:2)), str2double(sptime(3:4)),...
        str2double(sptime(5:6))];
    set(handles.edit_stoptime,'ForegroundColor','blue')
else
    beep
    disp(['Error: enter time in the hh-mm-ss format'])
    set(handles.edit_stoptime,'ForegroundColor','red')
    return
end

if ~isempty(handles.fname)
    Dmeg = crc_eeg_load(handles.fname);
    lgt = (nsamples(Dmeg)/fsample(Dmeg))/60^2;
    
    if isfield(Dmeg, 'info')
        if isfield(Dmeg.info, 'hour')
            hh = Dmeg.info.hour;
        else
            Dmeg.info = struct('date',[0 0 0],'hour',[0 0 0]);
            hh = Dmeg.info.hour;
            disp('Warning: time must be relative to recording start!')
        end
    else
        Dmeg.info = struct('date',[0 0 0],'hour',[0 0 0]);
        hh = Dmeg.info.hour;
        disp('Warning: time must be relative to recording start!')
    end
    
    datenumstart=datenum([0 0 0 hh]);
    datenumEnd=datenum([0 0 0 handles.stop(1) handles.stop(2) handles.stop(3)]);

    datediff=datevec(datenumEnd-datenumstart);
    if datediff(4)>lgt
        datenumEnd=datenum([0 0 1 handles.stop(1) handles.stop(2) handles.stop(3)]);
        datediff=datevec(datenumEnd-datenumstart);
    end
    Endpts=min(nsamples(Dmeg),max(1,round((datediff(4)*60^2+datediff(5)*60+datediff(6))*fsample(Dmeg))));
    if Endpts==nsamples(Dmeg)
        beep
        display('Warning: stop time larger than recording, considering size of file')
    end
else
    beep
    disp('Select EEG file first!')
    set(handles.EEGfilename,'ForegroundColor','red')
end
handles.Endpts=Endpts;

if handles.Begpts >= handles.Endpts
    beep
    disp('ERROR: The chosen starting point is invalid: it is placed after the chosen ending point');
    return
end
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_stoptime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_stoptime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cb_an_all.
function cb_an_all_Callback(hObject, eventdata, handles)
% hObject    handle to cb_an_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_an_all
all=get(hObject,'Value');
if all==1
    set(handles.cb_an_ft, 'Enable','off')
    set(handles.cb_an_3, 'Enable','off')
    set(handles.edit_scores, 'Enable','off')
    set(handles.edit_starttime, 'Enable','off')
    set(handles.edit_stoptime, 'Enable','off')
    handles.analyse=2;
else
    set(handles.cb_an_ft, 'Enable','on')
    set(handles.cb_an_3, 'Enable','on')
    set(handles.edit_scores, 'Enable','on')
    set(handles.edit_starttime, 'Enable','on')
    set(handles.edit_stoptime, 'Enable','on')
end
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in cb_an_3.
function cb_an_3_Callback(hObject, eventdata, handles)
% hObject    handle to cb_an_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_an_3
an3=get(hObject,'Value');
if an3==1
    set(handles.cb_an_all, 'Enable','off')
    set(handles.cb_an_ft, 'Enable','off')
    set(handles.edit_starttime, 'Enable','off')
    set(handles.edit_stoptime, 'Enable','off')
    set(handles.edit_scores, 'Enable','on')
    handles.analyse=3;
else
    set(handles.cb_an_all, 'Enable','on')
    set(handles.cb_an_ft, 'Enable','on')
    set(handles.edit_starttime, 'Enable','on')
    set(handles.edit_stoptime, 'Enable','on')
end
% Update handles structure
guidata(hObject, handles);



function edit_scores_Callback(hObject, eventdata, handles)
% hObject    handle to edit_scores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_scores as text
%        str2double(get(hObject,'String')) returns contents of edit_scores as a double
ss=get(hObject,'String');
% Update handles structure
guidata(hObject, handles);
st=zeros(1,size(ss,2));
for i=1:size(ss,2)
    try
        st(i)=str2double(ss(i));
    catch
        st(i)=0;
    end
end
handles.stagesp=st(~isnan(st));
if any(isnan(handles.stagesp)) || any(handles.stagesp==0)
    set(handles.edit_scores,'ForegroundColor','red')
    disp('Enter scores in the x x x format')
end
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_scores_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_scores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in scoringlist.
function scoringlist_Callback(hObject, eventdata, handles)
% hObject    handle to scoringlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns scoringlist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from scoringlist
handles.scorer=get(hObject,'Value');
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function scoringlist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scoringlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


