function varargout = crc_z3score_score(varargin)
% CRC_Z3SCORE_SCORE MATLAB code for crc_z3score_score.fig
%      CRC_Z3SCORE_SCORE, by itself, creates a new CRC_Z3SCORE_SCORE or raises the existing
%      singleton*.
%
%      H = CRC_Z3SCORE_SCORE returns the handle to a new CRC_Z3SCORE_SCORE or the handle to
%      the existing singleton*.
%
%      CRC_Z3SCORE_SCORE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CRC_Z3SCORE_SCORE.M with the given input arguments.
%
%      CRC_Z3SCORE_SCORE('Property','Value',...) creates a new CRC_Z3SCORE_SCORE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before crc_z3score_score_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to crc_z3score_score_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help crc_z3score_score

% Last Modified by GUIDE v2.5 23-Dec-2018 18:01:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @crc_z3score_score_OpeningFcn, ...
                   'gui_OutputFcn',  @crc_z3score_score_OutputFcn, ...
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


% --- Executes just before crc_z3score_score is made visible.
function crc_z3score_score_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to crc_z3score_score (see VARARGIN)

% Choose default command line output for crc_z3score_score
handles.D = varargin{1};
handles.flags = varargin{2};

set(handles.C3,'String',handles.D.chanlabels);
set(handles.C4,'String',handles.D.chanlabels);
set(handles.EL,'String',handles.D.chanlabels);
set(handles.ER,'String',handles.D.chanlabels);
set(handles.EM,'String',['Not Available', handles.D.chanlabels]);

set(handles.C3_Ref,'String',['No Reref', handles.D.chanlabels]);
set(handles.C4_Ref,'String',['No Reref', handles.D.chanlabels]);
set(handles.EL_Ref,'String',['No Reref', handles.D.chanlabels]);
set(handles.ER_Ref,'String',['No Reref', handles.D.chanlabels]);
set(handles.EM_Ref,'String',['No Reref', handles.D.chanlabels]);

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes crc_z3score_score wait for user response (see UIRESUME)
% uiwait(handles.gui);


% --- Outputs from this function are returned to the command line.
function varargout = crc_z3score_score_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in C3.
function C3_Callback(hObject, eventdata, handles)
% hObject    handle to C3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns C3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from C3


% --- Executes during object creation, after setting all properties.
function C3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to C3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in C4.
function C4_Callback(hObject, eventdata, handles)
% hObject    handle to C4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns C4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from C4


% --- Executes during object creation, after setting all properties.
function C4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to C4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in EL.
function EL_Callback(hObject, eventdata, handles)
% hObject    handle to EL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns EL contents as cell array
%        contents{get(hObject,'Value')} returns selected item from EL


% --- Executes during object creation, after setting all properties.
function EL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ER.
function ER_Callback(hObject, eventdata, handles)
% hObject    handle to ER (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ER contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ER


% --- Executes during object creation, after setting all properties.
function ER_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ER (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

EM = get(handles.EM,'Value')-1;
EL = get(handles.EL,'Value');
ER = get(handles.ER,'Value');
C3 = get(handles.C3,'Value');
C4 = get(handles.C4,'Value');

EM_Ref = get(handles.EM_Ref,'Value')-1;
EL_Ref = get(handles.EL_Ref,'Value')-1;
ER_Ref = get(handles.ER_Ref,'Value')-1;
C3_Ref = get(handles.C3_Ref,'Value')-1;
C4_Ref = get(handles.C4_Ref,'Value')-1;

settings_path = fullfile(getuserdir,'/z3license.mat');
h = waitbar(0,'Please wait, authenticating user...');

if fsample(handles.D) < 200 && EM ~= 0,
    close(h);
    errordlg(' Error message: Sampling rate too low','Minimum sampling rate must be 200Hz or more for Z3Score V2.');
    return
end

if fsample(handles.D) < 100 && EM == 0,
    close(h);
    errordlg(' Error message: Sampling rate too low','Minimum sampling rate must be 100Hz or more for Z3Score V1.');
    return
end

if exist(settings_path, 'file') == 2,
    load(settings_path);
    serverURL = settings.serverURL;
    email = settings.email;
    key = settings.key;
    if(isfield(settings,'isLocalServer'))
        isLocalServer = settings.isLocalServer;
    else
        settings.isLocalServer = 0;
        isLocalServer = settings.isLocalServer;
    end
else
    close(h);
    errordlg('No license found. Please go to settings and add your license information.','License missing');
    return
end



try
    response.status = 0;
    response.message = 'Cannot communicate with local server.';
    if(isLocalServer),
        client = connectZ3Score(serverURL);
        res = client.execute('ping', '');
        if(isa(res,'java.util.HashMap'))
            response.status = res.get('status');
            response.message = 'Local server is up and running.';
        end
    else
        response = loadjson(urlreadpost([serverURL '/check'],...
                                        {'email',email,'key',key}));
    end
catch
    close(h);
    errordlg('Error connecting server. Please check server address and internet connection.','Error Connecting Server');
    return
end

if response.status == 0,
    close(h);
    errordlg([' Error message: ' response.message],'License check failed.');
    return
end
waitbar(0.333,h,'Converting data to compressed feature set (CFS)...');
p = dcblock(0.2/fsample(handles.D));
b = [1 -1];                         % set up differentiator
a = [1 -p];                         % set up integrator

if EL_Ref == 0
    EL_Data = handles.D(EL,:,1);
else
    EL_Data = handles.D(EL,:,1) - handles.D(EL_Ref,:,1);
end

if ER_Ref == 0
    ER_Data = handles.D(ER,:,1);
else
    ER_Data = handles.D(ER,:,1) - handles.D(ER_Ref,:,1);
end

if C3_Ref == 0
    C3_Data = handles.D(C3,:,1);
else
    C3_Data = handles.D(C3,:,1) - handles.D(C3_Ref,:,1);
end

if C4_Ref == 0
    C4_Data = handles.D(C4,:,1);
else
    C4_Data = handles.D(C4,:,1) - handles.D(C4_Ref,:,1);
end

if EM == 0
    waitbar(0.5,h,'No EMG channel, using Z3Score V1...');
    EEGData = filter(b,a,[C3_Data; C4_Data; EL_Data; ER_Data]')';
    [stream, status, message, ~] = streamCFS(EEGData, fsample(handles.D));
elseif EM_Ref == 0
    waitbar(0.5,h,'Using Z3Score V2 (NEO)...');
    EM_Data = handles.D(EM,:,1);
    [stream, status, message, ~] = streamCFS_V2(C3_Data', C4_Data', EL_Data', ER_Data', EM_Data', fsample(handles.D), fsample(handles.D), fsample(handles.D));
else
    waitbar(0.5,h,'Using Z3Score V2 (NEO)...');
    EM_Data = handles.D(EM,:,1) - handles.D(EM_Ref,:,1);
    [stream, status, message, ~] = streamCFS_V2(C3_Data', C4_Data', EL_Data', ER_Data', EM_Data', fsample(handles.D), fsample(handles.D), fsample(handles.D));
end

waitbar(0.6,h,'Doing quality check...');
if status,
    close(h);
    errordlg(message,'Quality check failed');
    return;
end
    

waitbar(0.666,h,'Now uploading data and waiting for results...');

try
    if(isLocalServer),
        client = connectZ3Score(serverURL);
        res = client.execute('score', base64encode(stream));
        response.status = res.get('status');
        response.message = loadjson(res.get('score'));
        response.artifact = loadjson(res.get('artifact'));
    else
        
        response = loadjson(urlreadpost([serverURL '/score'], ...
            {'email',email,'key',key,'file',stream}));
    end
catch
    close(h);
    errordlg('Error connecting server. Please check server address and internet connection.','Error Connecting Server');
    return;
end

if response.status == 0,
    close(h);
    errordlg(['Error message:' response.message],'Error scoring data.'); 
    return;
end
waitbar(0.9,h,'Saving scores...');
scores = response.message;

%check if structure exits
if(~isfieldRecursive(handles.D,'CRC','score')),
    n = 1;
    handles.D.CRC.score = {};
else
    n = find(ismember(handles.D.CRC.score(2,:),'Oracle'),1,'first');
    if(isempty(n))
        n = size(handles.D.CRC.score,2)+1;
    end
end

%save relative confidence 
handles.D.relConfidence = [];
handles.D.relConfidence = scores(:,2);

handles.D.CRC.score{1,n} = scores(:,1)';
handles.D.CRC.score{2,n} = 'Oracle';
handles.D.CRC.score{3,n} = 30;
handles.D.CRC.score{4,n} = [0.005,length(scores(:,1))*30];
handles.D.CRC.score{5,n} = [];
handles.D.CRC.score{6,n} = [];
handles.D.CRC.score{7,n} = [];

%Save artifacts if available
if isfield(response,'artifact'),
    artifacts = response.artifact;
   
    n_art = find(ismember(handles.D.CRC.score(2,:),'Oracle_Artifact'),1,'first');
    if(isempty(n_art))
        n_art = size(handles.D.CRC.score,2)+1;
    end
    
    handles.D.CRC.score{1,n_art} = artifacts';
    handles.D.CRC.score{2,n_art} = 'Oracle_Artifact';
    handles.D.CRC.score{3,n_art} = 5;
    handles.D.CRC.score{4,n_art} = [0.005,length(artifacts)*5];
    handles.D.CRC.score{5,n_art} = [];
    handles.D.CRC.score{6,n_art} = [];
    handles.D.CRC.score{7,n_art} = [];
    
end

handles.flags.Dmeg{1} = save(handles.D);
waitbar(1,h,'Done...');
close(h);
h = msgbox('Data scored successfully. Sleep is scored by Oracle and Sleep Artifacts are scored by Oracle_Artifact.','Success');
waitfor(h);
handles.flags.showOracle = 1;
guidata(hObject, handles);
close(handles.gui);




% --- Executes when user attempts to close gui.
function gui_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
crc_dis_main(handles.flags);
delete(hObject);


% --- Executes during object creation, after setting all properties.
function C3_Ref_CreateFcn(hObject, eventdata, handles)
% hObject    handle to C3_Ref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in C4_Ref.
function C4_Ref_Callback(hObject, eventdata, handles)
% hObject    handle to C4_Ref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns C4_Ref contents as cell array
%        contents{get(hObject,'Value')} returns selected item from C4_Ref


% --- Executes during object creation, after setting all properties.
function C4_Ref_CreateFcn(hObject, eventdata, handles)
% hObject    handle to C4_Ref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in EL_Ref.
function EL_Ref_Callback(hObject, eventdata, handles)
% hObject    handle to EL_Ref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns EL_Ref contents as cell array
%        contents{get(hObject,'Value')} returns selected item from EL_Ref


% --- Executes during object creation, after setting all properties.
function EL_Ref_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EL_Ref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ER_Ref.
function ER_Ref_Callback(hObject, eventdata, handles)
% hObject    handle to ER_Ref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ER_Ref contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ER_Ref


% --- Executes during object creation, after setting all properties.
function ER_Ref_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ER_Ref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in C3_Ref.
function C3_Ref_Callback(hObject, eventdata, handles)
% hObject    handle to C3_Ref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns C3_Ref contents as cell array
%        contents{get(hObject,'Value')} returns selected item from C3_Ref


% --- Executes on selection change in EM.
function EM_Callback(hObject, eventdata, handles)
% hObject    handle to EM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns EM contents as cell array
%        contents{get(hObject,'Value')} returns selected item from EM


% --- Executes during object creation, after setting all properties.
function EM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in EM_Ref.
function EM_Ref_Callback(hObject, eventdata, handles)
% hObject    handle to EM_Ref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns EM_Ref contents as cell array
%        contents{get(hObject,'Value')} returns selected item from EM_Ref


% --- Executes during object creation, after setting all properties.
function EM_Ref_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EM_Ref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
