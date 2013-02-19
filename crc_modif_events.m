function varargout = crc_modif_events(varargin)
% MODIF_EVENTS MATLAB code for crc_modif_events.fig
%      MODIF_EVENTS, by itself, creates a new MODIF_EVENTS or raises the existing
%      singleton*.
%
%      H = MODIF_EVENTS returns the handle to a new MODIF_EVENTS or the handle to
%      the existing singleton*.
%
%      MODIF_EVENTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MODIF_EVENTS.M with the given input arguments.
%
%      MODIF_EVENTS('Property','Value',...) creates a new MODIF_EVENTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before crc_modif_events_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to crc_modif_events_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help crc_modif_events

% Last Modified by GUIDE v2.5 05-Jun-2012 15:37:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @crc_modif_events_OpeningFcn, ...
                   'gui_OutputFcn',  @crc_modif_events_OutputFcn, ...
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


% --- Executes just before crc_modif_events is made visible.
function crc_modif_events_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to crc_modif_events (see VARARGIN)

% Choose default command line output for crc_modif_events
handles.output = hObject;

%load flag from crc_main
handles.Dmeg = varargin{1}.Dmeg;
handles.type = varargin{1}.type;
handles.file = varargin{1}.file;
handles.index = varargin{1}.index;
handles.currentscore = varargin{1}.user; 

%pannel color available
pannelcolor{1} = 'yellow';
pannelcolor{2} = 'magenta';
pannelcolor{3} = 'ciel';
pannelcolor{4} = 'red';
pannelcolor{5} = 'green';

handles.color = pannelcolor;

set(handles.colorevent,'String',pannelcolor,'Value',4)

%type available
if ~isempty(handles.type)
    set(handles.typeofevent,'String',handles.type(:,1),'Value',size(handles.type(:,1),1))
end
%initialization 
set(handles.colorevent,'enable','on','visible','on');
if ~isempty(handles.type)
    set(handles.name,'String',handles.type(end,1))
else 
    set(handles.name,'String','name')
end

set(gcf,'CloseRequestFcn',@my_closefcn)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes crc_modif_events wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function my_closefcn(hObject, handles)

handles = guidata(hObject);

flags.Dmeg      =   handles.Dmeg;
flags.file      =   handles.file;
flags.index     =   handles.index;
flags.type      =   handles.type;
flags.user      =   handles.currentscore;
flags.scoresleep = 1;

if isfield(handles,'delmap')
    flags.delmap=handles.delmap;
end

crc_dis_main(flags);

try
    close(handles.figure1)
end

% --- Outputs from this function are returned to the command line.
function varargout = crc_modif_events_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in typeofevent.
function typeofevent_Callback(hObject, eventdata, handles)
% hObject    handle to typeofevent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns typeofevent contents as cell array
%        contents{get(hObject,'Value')} returns selected item from typeofevent

if ~isempty(handles.type)
contents = cellstr(get(hObject,'String'));
name = contents{get(hObject,'Value')};
val = get(hObject,'Value');
color = handles.type(val,2);
[col intcol intpan] = intersect(color,handles.color);

set(handles.colorevent,'Value',intpan)
set(handles.name,'String',name)
end
% --- Executes during object creation, after setting all properties.
function typeofevent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to typeofevent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function name_Callback(hObject, eventdata, handles)
% hObject    handle to name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of name as text
%        str2double(get(hObject,'String')) returns contents of name as a double


% --- Executes during object creation, after setting all properties.
function name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in savevent.
function savevent_Callback(hObject, eventdata, handles)
% hObject    handle to savevent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

name = get(handles.name,'String');

if ischar(name)
    name = {name};
end

intcol  = get(handles.colorevent,'Value');
col = handles.color(intcol);
intevent = get(handles.typeofevent,'Value');

%change all the events whose name's changing
ev = events(handles.Dmeg{1});
tev = {ev(:).type};
itev = find(strcmpi(tev,handles.type(intevent,1)));

if iscell(name)
    namev = char(name);
else 
    namev = name;
end

for iitev = 1 : length(itev)
    ev(itev(iitev)).type = namev;
end
handles.Dmeg{1} = events(handles.Dmeg{1},1,ev);
save(handles.Dmeg{1});

handles.type(intevent,1) = name;
handles.type(intevent,2) = col;

set(handles.typeofevent,'String',handles.type(:,1));

handles.Dmeg{1}.CRC.Event = handles.type;

save(handles.Dmeg{1});

guidata(hObject,handles)

% --- Executes on button press in suppressevent.
function suppressevent_Callback(hObject, eventdata, handles)
% hObject    handle to suppressevent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

intevent = get(handles.typeofevent,'Value');
name = handles.type(intevent,1);

%to prevent error if the event selected is the last one
if intevent == size(handles.type,1) 
    new  =  intevent-1;
else 
    new  = intevent;
end

ButtonName = questdlg(['Are you sure do you want to suppress the event "',char(name),'" ?'], ...
                     'Suppression of event', ...
                     'YES', 'NO', 'NO');
switch ButtonName,
    case 'YES'
     change(1 : intevent-1,:) = handles.type(1 : intevent-1,:) ;
     change(intevent : size(handles.type,1)-1,:) = handles.type(intevent+1 : size(handles.type,1),:);
    case 'NO'
        
end
handles.type = change;
set(handles.typeofevent,'String',handles.type(:,1),'Value',new)

%change in events
D = handles.Dmeg{1};
ev = events(D);
nev = 1;
while nev <= size(ev,2)
    if strcmpi(ev(nev).type,name)
        stock1 = struct([]);
        stock2 = struct([]);
        stock1 = ev(1:nev-1);
        stock2 = ev(nev + 1 : size(ev,2));
        ev = [stock1 stock2];
    else 
        nev = nev + 1;
    end
end

%save
D = events(D,1,ev);
save(D);
handles.Dmeg{1} = D;

guidata(hObject,handles)


% --- Executes on selection change in colorevent.
function colorevent_Callback(hObject, eventdata, handles)
% hObject    handle to colorevent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns colorevent contents as cell array
%        contents{get(hObject,'Value')} returns selected item from colorevent


% --- Executes during object creation, after setting all properties.
function colorevent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to colorevent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function return_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

flags.type  = handles.type;
flags.Dmeg  = handles.Dmeg;
flags.file  = handles.file;
flags.index = handles.index;
flags.user  = handles.currentscore;
flags.scoresleep = 1;

crc_dis_main(flags)

try
    close(handles.figure1)
end
