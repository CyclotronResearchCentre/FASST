function varargout = crc_sws_detection(varargin)
% CRC_SWS_DETECTION M-file for crc_sws_detection.fig
%
% GUI to select SWS detection parameters
%
%      CRC_SWS_DETECTION, by itself, creates a new CRC_SWS_DETECTION or raises the existing
%      singleton*.
%
%      H = CRC_SWS_DETECTION returns the handle to a new CRC_SWS_DETECTION or the handle to
%      the existing singleton*.
%
%      CRC_SWS_DETECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CRC_SWS_DETECTION.M with the given input arguments.
%
%      CRC_SWS_DETECTION('Property','Value',...) creates a new CRC_SWS_DETECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before crc_sws_detection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property
%      application
%      stop.  All inputs are passed to crc_sws_detection_OpeningFcn via varargin.
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

% Last Modified by GUIDE v2.5 13-Jan-2011 16:34:55


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @crc_sws_detection_OpeningFcn, ...
                   'gui_OutputFcn',  @crc_sws_detection_OutputFcn, ...
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


% --- Executes just before crc_sws_detection is made visible.
function crc_sws_detection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to crc_sws_detection (see VARARGIN)

% Choose default command line output for crc_sws_detection
handles.output = hObject;

%choose defaults
handles.analyse=1;
handles.highfc=0.2;
handles.lowfc=4;
handles.roisel=1;
handles.review=0;
handles.sensauto=1;
handles.sensfname=[];
handles.maps=1;
handles.fmri=0;
set(handles.edit_TR, 'Enable','off')
set(handles.edit_marker, 'Enable','off')
handles.Begpts=0;
handles.Endpts=1;
handles.fname=[];
handles.marker=[];
handles.TR=[];
handles.reref=0;
handles.scorer=1;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes crc_sws_detection wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = crc_sws_detection_OutputFcn(hObject, eventdata, handles) 
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
try
    D = crc_eeg_load(deblank(handles.fname));
    if isfield(D,'CRC') && isfield(D.CRC,'score')
        for i=1:size(D.CRC.score,2)
            sc{i}=D.CRC.score{2,i};
        end
        set(handles.scoringlist,'enable','on','String',sc);
    else
        set(handles.scoringlist,'enable','off');
    end
catch
    beep
    disp('Not a valid filename')
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
try
    D = crc_eeg_load(deblank(handles.fname));
    if isfield(D,'CRC') && isfield(D.CRC,'score')
        for i=1:size(D.CRC.score,2)
            sc{i}=D.CRC.score{2,i};
        end
        set(handles.scoringlist,'enable','on','String',sc);
    else
        set(handles.scoringlist,'enable','off');
    end
catch
    beep
    disp('Not a valid filename')
end
% Update handles structure
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
            DDmeg.info=struct('date',[0 0 0],'hour',[0 0 0]);
            hh = Dmeg.info.hour;
            disp('Warning: time must be relative to recording start!')
        end
    else
        Dmeg.info = struct('date',[0 0 0],'hour',[0 0 0]);
        hh = Dmeg.info.hour;
        disp('Warning: time must be relative to recording start!')
    end
    
    datenumstart =datenum([0 0 0 hh]);
    datenumEnd   =datenum([0 0 0 handles.stop(1) handles.stop(2) handles.stop(3)]);

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


% --- Executes on button press in cb_roi_auto.
function cb_roi_auto_Callback(hObject, eventdata, handles)
% hObject    handle to cb_roi_auto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_roi_auto
handles.roisel=get(hObject,'Value');
if handles.roisel==1
    set(handles.cb_roi_manual, 'Enable','off')
    set(handles.edit_numroi, 'Enable','off')
else
    set(handles.cb_roi_manual, 'Enable','on')
    set(handles.edit_numroi, 'Enable','on')
end
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in cb_roi_manual.
function cb_roi_manual_Callback(hObject, eventdata, handles)
% hObject    handle to cb_roi_manual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_roi_manual
roisel=get(hObject,'Value');
if roisel==1
    handles.roisel=0;
    set(handles.cb_roi_auto, 'Enable','off')
else
    set(handles.cb_roi_auto, 'Enable','on')
    set(handles.edit_numroi, 'Enable','on')
end
% Update handles structure
guidata(hObject, handles);



function edit_numroi_Callback(hObject, eventdata, handles)
% hObject    handle to edit_numroi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_numroi as text
%        str2double(get(hObject,'String')) returns contents of edit_numroi as a double
handles.numroi=str2double(get(hObject,'String'));
if handles.roisel==0
    if isfield (handles,'fname')
        D = crc_eeg_load(handles.fname);
    else
        beep
        disp('Select EEG file first!')
        return
    end
    handles.name_roi=cell(1,handles.numroi);
    for ir=1:handles.numroi
        handles.name_roi{ir}=spm_input(['Name of ROI ',num2str(ir)],'!+1','s');
    end
    close gcf
    varE = struct('file',[],'spmstruct',[],'Dmeg',[]);
    varE.spmstruct  = struct(D);
    varE.file       = handles.fname;
    varE.Dmeg       = D;
    handles.sel_roi=cell(1,handles.numroi);
    for iroi = 1:handles.numroi %regions of interest
        global gindex %handles
        h=dis_SWS_selchan(varE);
        waitfor(h,'BeingDeleted')
        handles.sel_ROI{iroi}=gindex;
    end
end
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_numroi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_numroi (see GCBO)
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

% --- Executes on button press in cb_sens_auto.
function cb_sens_auto_Callback(hObject, eventdata, handles)
% hObject    handle to cb_sens_auto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_sens_auto
handles.sensauto=get(hObject,'Value');
if handles.sensauto==1
    set(handles.cb_sens_file, 'Enable','off')
    set(handles.pushbutton_mapspot, 'Enable','off')
else
    set(handles.cb_sens_file, 'Enable','on')
    set(handles.pushbutton_mapspot, 'Enable','on')
end
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in cb_sens_file.
function cb_sens_file_Callback(hObject, eventdata, handles)
% hObject    handle to cb_sens_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_sens_file
sensfile=get(hObject,'Value');
if sensfile==1
    set(handles.cb_sens_auto, 'Enable','off')
    handles.sensauto=0;
    handles.sensfname=spm_select(1, 'any', 'Select electrodes positioning file','' ,pwd,'.*');
else
    set(handles.cb_sens_auto, 'Enable','on')
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_sensfile.
function pushbutton_sensfile_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_sensfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.sensfname=spm_select(1, 'any', 'Select electrodes positioning file','' ,pwd,'.*');
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_mapsdelays.
function pushbutton_mapsdelays_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_mapsdelays (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.maps=1;
set(handles.pushbutton_mapsdelays,'BackgroundColor',[186/255,212/255,244/255])
set(handles.pushbutton_mapspot,'BackgroundColor',[236/255,233/255,216/255])
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_mapspot.
function pushbutton_mapspot_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_mapspot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.maps=2;
set(handles.pushbutton_mapspot,'BackgroundColor',[186/255,212/255,244/255])
set(handles.pushbutton_mapsdelays,'BackgroundColor',[236/255,233/255,216/255])
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in radiobutton_fmri.
function radiobutton_fmri_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_fmri (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_fmri
handles.fmri=get(hObject,'Value');
if handles.fmri==1
    set(handles.edit_TR, 'Enable','on')
    set(handles.edit_marker, 'Enable','on')
else
    set(handles.edit_TR, 'Enable','off')
    set(handles.edit_marker, 'Enable','off')
end
% Update handles structure
guidata(hObject, handles);

function edit_TR_Callback(hObject, eventdata, handles)
% hObject    handle to edit_TR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_TR as text
%        str2double(get(hObject,'String')) returns contents of edit_TR as a double
handles.TR=str2double(get(hObject,'String'));
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_TR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_TR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_marker_Callback(hObject, eventdata, handles)
% hObject    handle to edit_marker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_marker as text
%        str2double(get(hObject,'String')) returns contents of edit_marker as a double
handles.marker=get(hObject,'String');
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_marker_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_marker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_OK.
function pushbutton_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
D=crc_SWS_detect(handles);



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
handles.stagesw=st(~isnan(st));
if any(isnan(handles.stagesw)) || any(handles.stagesw==0)
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


