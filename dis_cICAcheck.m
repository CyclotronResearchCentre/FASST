function varargout = dis_cICAcheck(varargin)
% DIS_CICACHECK M-file for dis_cICAcheck.fig
%      DIS_CICACHECK, by itself, creates a new DIS_CICACHECK or raises the existing
%      singleton*.
%
%      H = DIS_CICACHECK returns the handle to a new DIS_CICACHECK or the handle to
%      the existing singleton*.
%
%      DIS_CICACHECK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DIS_CICACHECK.M with the given input arguments.
%
%      DIS_CICACHECK('Property','Value',...) creates a new DIS_CICACHECK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dis_main_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dis_cICAcheck_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%__________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by Y. Leclercq & C. Phillips, 2008.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

% Edit the frqabv text to modify the response to help dis_cICAcheck

% Last Modified by GUIDE v2.5 14-May-2010 16:52:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dis_cICAcheck_OpeningFcn, ...
                   'gui_OutputFcn',  @dis_cICAcheck_OutputFcn, ...
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


% --- Executes just before dis_cICAcheck is made visible.
function dis_cICAcheck_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dis_cICAcheck (see VARARGIN)

crcdef = crc_get_defaults('one');
set(0,'CurrentFigure',handles.figure1);
warning('off')

% Choose default command line output for dis_cICAcheck
handles.output = hObject;
handles.file = varargin{1}.file;
handles.figz = 0;
handles.index = varargin{1}.index;
handles.Dmeg = varargin{1}.Dmeg{1};

if isfield(handles.Dmeg,'info')
    if isfield(handles.Dmeg.info,'hour')
        handles.offset = handles.Dmeg.info.hour(1) * 60^2 + ...
                         handles.Dmeg.info.hour(2)*60 + ...
                         handles.Dmeg.info.hour(3);
    end
end

handles.winsize = crcdef.winsize;
handles.scale   = crcdef.scale;

% Define filter

handles.filter.other = crcdef.filtEEG;
handles.filter.EMG = [crcdef.filtEMG(1) ...
          min(crcdef.filtEMG(2),round(fsample(handles.Dmeg)/2*.9))];
handles.filter.EOG = crcdef.filtEOG;

set(handles.slider1,...
    'Max', nsamples(handles.Dmeg)/fsample(handles.Dmeg)-handles.winsize,...
    'Value',1/fsample(handles.Dmeg),...
    'Min',1/fsample(handles.Dmeg))

try
    set(handles.slider1,'SliderStep',[handles.winsize/(nsamples(handles.Dmeg)/ ...
                                      fsample(handles.Dmeg)-handles.winsize) 0.1])
end

Nchdisp=str2double(get(handles.NbreChan,'String'));
set(handles.NbreChan,'String', num2str(min(Nchdisp, length(handles.index))))
Nchdisp=str2double(get(handles.NbreChan,'String'));

set(handles.Chanslider,...
    'Min',1,...
    'Max',length(handles.index)-Nchdisp+1,...
    'Value',1,...
    'SliderStep', [min(1/(length(handles.index)-Nchdisp),1) min(1/(length(handles.index)-Nchdisp),1)]);

set(handles.totaltime,'String',['/ ' ...
    num2str(round(nsamples(handles.Dmeg)/fsample(handles.Dmeg)))]);
set(handles.currenttime,'String',num2str(round(handles.winsize/2)));
set(handles.filename,'String',fname(handles.Dmeg));
set(handles.frqabv,'String',num2str(handles.filter.other(2)));
set(handles.upemg,'String',num2str(handles.filter.EMG(2)));

% Popmenu setting
popmenustring = chanlabels(handles.Dmeg) ;

if sum(strcmp(chanlabels(handles.Dmeg),'REF2'))
    popmenustring = [popmenustring 'MEAN OF REF'];
end

if and(sum(strcmp(chanlabels(handles.Dmeg),'M1')), ...
       sum(strcmp(chanlabels(handles.Dmeg),'M2')))
    popmenustring = [popmenustring 'M1-M2'];
end

popmenustring = [popmenustring 'REF1'];

set(handles.EEGpopmenu,...
    'String',popmenustring,...
    'Value',length(popmenustring))

set(handles.otherpopmenu,...
    'String',popmenustring,...
    'Value',length(popmenustring))

popmenustring = [popmenustring 'BIPOLAR'];

set(handles.EOGpopmenu,...
    'String',popmenustring,...
    'Value',length(popmenustring))

set(handles.EMGpopmenu,...
    'String',popmenustring,...
    'Value',length(popmenustring))

set(handles.popcorrmx,'String',1:length(handles.Dmeg.CRC.cICA{3}))

load CRC_electrodes.mat;
handles.names     = names;
handles.crc_types = crc_types;
handles.cmx       = 1:length(handles.Dmeg.CRC.cICA{3});

mainplot(handles)

xtick = get(handles.axes1,'XTick');

if isfield(handles.Dmeg,'info')
    if isfield(handles.Dmeg.info,'hour')
        xtick = mod(xtick + handles.offset,24*60^2);
    end
end

[time string] = crc_time_converts(xtick);
set(handles.axes1,'XTickLabel',string)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dis_cICAcheck wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = dis_cICAcheck_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

mainplot(handles)

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
handles.winsize = min(str2double(get(hObject,'String')), ...
        round(nsamples(handles.Dmeg)/(2*fsample(handles.Dmeg))));
set(handles.edit1,'String',num2str(handles.winsize));

set(handles.slider1,'Max',nsamples(handles.Dmeg)/ ...
                                fsample(handles.Dmeg)-handles.winsize)
set(handles.slider1,'Min',1/fsample(handles.Dmeg))
set(handles.slider1,'SliderStep', ...
                    [handles.winsize/(nsamples(handles.Dmeg)/ ...
                     fsample(handles.Dmeg)-handles.winsize) 0.1])

mainplot(handles)

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
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

if not(str2double(get(hObject,'String'))>0)&not(str2double(get(hObject,'String'))<0)
    set(handles.edit2,'String',num2str(handles.scale));
else
    handles.scale=str2double(get(hObject,'String'));
end

% Update handles structure
guidata(hObject, handles);

mainplot(handles)

i=length(handles.index);
set(handles.axes1,'YTick',[handles.scale/2:handles.scale/2:i*handles.scale+handles.scale/2]);

ylabels=[num2str(round(handles.scale/2))];
for j=1:length(handles.index)
    ylabels=[ylabels chanlabels(handles.Dmeg,handles.index(j))];
    ylabels=[ylabels num2str(round(handles.scale/2))];
end
set(handles.axes1,'YTickLabel',ylabels);

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

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

flags.index = handles.index;
flags.Dmeg  = handles.Dmeg;
flags.file  = handles.file;
dis_selchan(flags);

delete(handles.figure1)
try
    close(handles.figz)
end

function currenttime_Callback(hObject, eventdata, handles)
% hObject    handle to currenttime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of currenttime as text
%        str2double(get(hObject,'String')) returns contents of currenttime as a double

slidval = str2num(get(hObject,'String'));
slidval = max(slidval - handles.winsize/2,1/fsample(handles.Dmeg));
slidval = min(slidval,nsamples(handles.Dmeg)/ ...
                        fsample(handles.Dmeg)-handles.winsize);
set(handles.slider1,'Value',slidval)

mainplot(handles)

% --- Executes during object creation, after setting all properties.
function currenttime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to currenttime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dis_selchan;
delete(handles.figure1);
try
    close(handles.figz)
end

% --- Executes on button press in push_wfile.
function push_wfile_Callback(hObject, eventdata, handles)
% hObject    handle to push_wfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val=get(handles.popcorrmx,'Value');

handles.cmx(:) = val;
guidata(hObject, handles);

function frqabv_Callback(hObject, eventdata, handles)
% hObject    handle to frqabv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frqabv as text
%        str2double(get(hObject,'String')) returns contents of frqabv as a double

frbe=str2double(get(hObject,'String'));

if not(isnan(frbe))
    if frbe>fsample(handles.Dmeg)/2
        handles.filter.other(2)=fsample(handles.Dmeg)/2;
        set(hObject,'String',num2str(fsample(handles.Dmeg)/2));
    else
        handles.filter.other(2) = frbe;
    end

else
    beep
    set(hObject,'String',num2str(handles.filter.other(2)))
    return
end

guidata(hObject, handles);
mainplot(handles)
    
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function frqabv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frqabv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function frqblw_Callback(hObject, eventdata, handles)
% hObject    handle to frqblw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frqblw as text
%        str2double(get(hObject,'String')) returns contents of frqblw as a double

frbe=str2double(get(hObject,'String'));

if not(isnan(frbe))
    if frbe<0.0001
        handles.filter.other(1)=0.0001;
        set(hObject,'String',num2str(0.0001));
    else
        handles.filter.other(1) = frbe;
    end
else
    beep
    set(hObject,'String',num2str(handles.frqbelow))
    return
end

guidata(hObject, handles);

mainplot(handles)

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function frqblw_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frqblw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in checkfil.
function checkfil_Callback(hObject, eventdata, handles)
% hObject    handle to checkfil (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkfil
mainplot(handles)

% --- Executes on button press in EMGfil.
function EMGfil_Callback(hObject, eventdata, handles)
% hObject    handle to EMGfil (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of EMGfil
mainplot(handles)

function upemg_Callback(hObject, eventdata, handles)
% hObject    handle to upemg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of upemg as text
%        str2double(get(hObject,'String')) returns contents of upemg as a double

frbe=str2double(get(hObject,'String'));

if not(isnan(frbe))
    if frbe>fsample(handles.Dmeg)/2*.9
        handles.filter.EMG(2) = fsample(handles.Dmeg)/2*.9;
        set(hObject,'String',num2str(fsample(handles.Dmeg)/2*.9));
    else
        handles.filter.EMG(2) = frbe;
    end
else
    beep
    set(hObject,'String',num2str(handles.filter.EMG(2)))
    return
end

guidata(hObject, handles);

mainplot(handles)

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function upemg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to upemg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function downemg_Callback(hObject, eventdata, handles)
% hObject    handle to downemg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of downemg as text
%        str2double(get(hObject,'String')) returns contents of downemg as a double
frbe=str2double(get(hObject,'String'));

if not(isnan(frbe))
    if frbe<0.0001
        handles.filter.EMG(1)=0.0001;
        set(hObject,'String',num2str(0.0001));
    else
        handles.filter.EMG(1) = frbe;
    end
else
    beep
    set(hObject,'String',num2str(handles.filter.EMG(1)))
    return
end

guidata(hObject, handles);

mainplot(handles)

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function downemg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to downemg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in EOGpopmenu.
function EOGpopmenu_Callback(hObject, eventdata, handles)
% hObject    handle to EOGpopmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns EOGpopmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from EOGpopmenu

mainplot(handles)

% --- Executes during object creation, after setting all properties.
function EOGpopmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EOGpopmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in EMGpopmenu.
function EMGpopmenu_Callback(hObject, eventdata, handles)
% hObject    handle to EMGpopmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns EMGpopmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from EMGpopmenu

mainplot(handles)

% --- Executes during object creation, after setting all properties.
function EMGpopmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EMGpopmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in EEGpopmenu.
function EEGpopmenu_Callback(hObject, eventdata, handles)
% hObject    handle to EEGpopmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns EEGpopmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from EEGpopmenu
mainplot(handles)

% --- Executes during object creation, after setting all properties.
function EEGpopmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EEGpopmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in otherpopmenu.
function otherpopmenu_Callback(hObject, eventdata, handles)
% hObject    handle to otherpopmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns otherpopmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from otherpopmenu
mainplot(handles)

% --- Executes during object creation, after setting all properties.
function otherpopmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to otherpopmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --------------------------------------------------------------------
function Pwr_spect_Callback(hObject, eventdata, handles)
% hObject    handle to Pwr_spect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function Cmp_Pwr_Sp_Callback(hObject, eventdata, handles)
% hObject    handle to Cmp_Pwr_Sp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hand=get(handles.axes1,'Children');
Mouse=get(handles.axes1,'CurrentPoint');

Chan=ceil((Mouse(1,2)-handles.scale/2)/handles.scale);
Chan=min(max(1,Chan),length(handles.index));

if handles.figz~=0
    z=handles.figz;
else
    z=figure;
end

figure(z)
axs=get(handles.figz,'Children');

cleargraph(handles.figz)

fs=fsample(handles.Dmeg);

slidval = get(handles.slider1,'Value');

tdeb=round(slidval*fsample(handles.Dmeg));
temps=tdeb:1:min(tdeb+(fsample(handles.Dmeg)*handles.winsize), ...
                 nsamples(handles.Dmeg));
toshow=temps;
temps=temps/fsample(handles.Dmeg);
index=handles.index;

%%
mxnb=get(handles.popcorrmx,'Value');

badindex=sort(handles.Dmeg.CRC.cICA{4}); % gives the badindex
selected_chan = setdiff(1:nchannels(handles.Dmeg), badindex);

plotdat=handles.Dmeg(:,toshow); % load current data;
cormx=handles.Dmeg.CRC.cICA{3}{mxnb};

plotdat(selected_chan,:)=cormx*plotdat(selected_chan,:);

%%
for i=Chan
    hold on
    [dumb1,dumb2,index2] = ...
        intersect(upper(chanlabels(handles.Dmeg,index(i))),handles.names);
    if abs(handles.crc_types(index2))>1
        if handles.crc_types(index2)>0
            [dumb1,index1,dumb2] = ...
                intersect(upper(chanlabels(handles.Dmeg)), ...
                          upper(handles.names(handles.crc_types(index2))));
            bipolar = plotdat(index(i),:) - ...
                        handles.Dmeg(index1,toshow);
            X=bipolar;
            Col=[0.2 0.9 0.5];
        else
            range=max(plotdat(index(i),:)) - ...
                    min(plotdat(index(i),:));
            X=(handles.scale)*plotdat(index(i),:)/range;
            Col=[1 0 0];
        end
    else
        X=plotdat(index(i),:);
        Col=[0 0 1];
    end
    
    X=filterforspect(handles,X,[0.001 fsample(handles.Dmeg)/2]);
    
    [P,F] = pwelch(X,[],[],[],fs);
    P=log(P);    
    plot(F,P,'Color',Col)

end
    
grid on

titre=['Power on ' chanlabels(handles.Dmeg,handles.index(Chan))];
title(titre)
ylabel('Log of power')
xlabel('Frequency in Hz')
xlim([0 20])

handles.figz=z;

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function Export_Callback(hObject, eventdata, handles)
% hObject    handle to Export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.figz~=0
    z=handles.figz;
else
    z=figure;
end

figure(z)

cleargraph(z)

mainplot(handles)

ax = get(z, 'Children');

i=length(handles.index);

set(ax,'YTick',[handles.scale/2:handles.scale/2:i*handles.scale+handles.scale/2]);

ylabels=[num2str(round(handles.scale/2))];
for j=1:length(handles.index)
    ylabels=[ylabels chanlabels(handles.Dmeg,handles.index(j))];
    ylabels=[ylabels num2str(round(handles.scale/2))];
end
set(ax,'YTickLabel',ylabels);

xtick = get(ax,'XTick');
if isfield(handles.Dmeg.info,'hour')
    xtick = mod(xtick + handles.offset,24*60^2);
end
[time string] = crc_time_converts(xtick);
set(ax,'XTickLabel',string)

handles.figz=z;

% Update handles structure
guidata(hObject, handles);

return

% --- Executes on button press in EOGfil.
function EOGfil_Callback(hObject, eventdata, handles)
% hObject    handle to EOGfil (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of EOGfil
mainplot(handles)

function downeog_Callback(hObject, eventdata, handles)
% hObject    handle to downeog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of downeog as text
%        str2double(get(hObject,'String')) returns contents of downeog as a double

frbe=str2double(get(hObject,'String'));

if not(isnan(frbe))
    if frbe<0.0001
        handles.filter.EOG(1)=0.0001;
        set(hObject,'String',num2str(0.0001));
    else
        handles.filter.EOG(1) = frbe;
    end
else
    beep
    set(hObject,'String',num2str(handles.filter.EOG(1)))
    return
end

guidata(hObject, handles);

mainplot(handles)

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function downeog_CreateFcn(hObject, eventdata, handles)
% hObject    handle to downeog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function upeog_Callback(hObject, eventdata, handles)
% hObject    handle to upeog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of upeog as text
%        str2double(get(hObject,'String')) returns contents of upeog as a double
frbe=str2double(get(hObject,'String'));

if not(isnan(frbe))
    if frbe>fsample(handles.Dmeg)/2
        handles.filter.EOG(2)=fsample(handles.Dmeg)/2;
        set(hObject,'String',num2str(fsample(handles.Dmeg)/2));
    else
        handles.filter.EOG(2) = frbe;
    end
else
    beep
    set(hObject,'String',num2str(handles.filter.EOG(2)))
    return
end

guidata(hObject, handles);

mainplot(handles)

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function upeog_CreateFcn(hObject, eventdata, handles)
% hObject    handle to upeog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function Chanslider_Callback(hObject, eventdata, handles)
% hObject    handle to Chanslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

slidval = get(handles.Chanslider,'Value');
set(handles.Chanslider,'Value',slidval-rem(slidval,1))
mainplot(handles);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Chanslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Chanslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function NbreChan_Callback(hObject, eventdata, handles)
% hObject    handle to NbreChan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NbreChan as text
%        str2double(get(hObject,'String')) returns contents of NbreChan as a double
Nchdisp=str2double(get(handles.NbreChan,'String'));
if ~isnan(Nchdisp)
    set(handles.NbreChan,'String', num2str(min(Nchdisp, length(handles.index))))
    Nchdisp=str2double(get(handles.NbreChan,'String'));

    set(handles.Chanslider,...
        'Min',1,...
        'Max',length(handles.index)-Nchdisp+1,...
        'Value',1,...
        'SliderStep', [min(1/(length(handles.index)-Nchdisp),1) min(1/(length(handles.index)-Nchdisp),1)]);
    mainplot(handles);
else
    beep
    set(handles.NbreChan,'String', num2str(10))
end

% --- Executes during object creation, after setting all properties.
function NbreChan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NbreChan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popcorrmx.
function popcorrmx_Callback(hObject, eventdata, handles)
% hObject    handle to popcorrmx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popcorrmx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popcorrmx

slidval = get(handles.slider1,'Value');
mxnb=get(handles.popcorrmx,'Value');
mxtarget=floor(slidval/(handles.Dmeg.CRC.cICA{1})+1); % give the matrix number to be used
handles.cmx(mxtarget)=mxnb;

guidata(hObject, handles);

mainplot(handles);

% --- Executes during object creation, after setting all properties.
function popcorrmx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popcorrmx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in push_cICA.
function push_cICA_Callback(hObject, eventdata, handles)
% hObject    handle to push_cICA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Di = handles.Dmeg;
ch_Mi = handles.cmx;
crc_par_cICAsav(Di,ch_Mi);
 
close(handles.figure1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mainplot(handles)

slidval = get(handles.slider1,'Value');

tdeb=round(slidval*fsample(handles.Dmeg));
temps=tdeb:1:min(tdeb+(fsample(handles.Dmeg)*handles.winsize), ...
                    nsamples(handles.Dmeg));
toshow=temps;
temps=temps/fsample(handles.Dmeg);

Chanslidval = get(handles.Chanslider,'Value');
slidpos=Chanslidval-rem(Chanslidval,1);

NbreChandisp=str2double(get(handles.NbreChan,'String'));

index=handles.index(slidpos:1:slidpos+NbreChandisp-1);
cleargraph(handles.figure1)

%%
mxtarget=floor(slidval/(handles.Dmeg.CRC.cICA{1})+1); % give the matrix number to be used
set(handles.popcorrmx,'Value',handles.cmx(min(mxtarget,length(handles.cmx))));



mxnb=get(handles.popcorrmx,'Value');

badindex=sort(handles.Dmeg.CRC.cICA{4}); % gives the badindex
selected_chan = setdiff(1:nchannels(handles.Dmeg), badindex);

plotdat=handles.Dmeg(:,toshow); % load current data;
cormx=handles.Dmeg.CRC.cICA{3}{mxnb};

plotdat(selected_chan,:)=cormx*plotdat(selected_chan,:);

%%

for i=1:length(index)
    hold on
    [dumb1,dumb2,index2] = ...
        intersect(upper(chanlabels(handles.Dmeg,index(i))),handles.names);
    
    crc_type = handles.crc_types(index2);
    
    if length(crc_type)==0
        crc_type = -1;
    end
    
%    switch crc_type
        
%        case -1
         testcond = chanlabels(handles.Dmeg,index(i));
            if  strfind(testcond{:},'EMG')
                
                contents = get(handles.EMGpopmenu,'String');
                selchan=upper(contents{get(handles.EMGpopmenu,'Value')});
                
                switch  selchan
                
                    case    'REF1'
                            
                        plt=plot(temps,handles.scale*i+plotdat(index(i),:));

                    case    'MEAN OF REF'
                        
                        basedata = plotdat(index(i),:);
                        ref2idx = find(strcmp(chanlabels(handles.Dmeg),'REF2'));
                        scnddata = plotdat(index(i),:) - ...
                                    plotdat(ref2idx,:);
                        
                        toplotdat = mean([basedata ; scnddata]);
                        
                        plt=plot(temps,handles.scale*i+toplotdat);
                    
                    case    'M1-M2'
                        
                        basedata = plotdat(index(i),:);
                        M1idx = find(strcmp(chanlabels(handles.Dmeg),'M1'));
                        M2idx = find(strcmp(chanlabels(handles.Dmeg),'M2'));
                        meanM = mean([plotdat(M1idx,:) ; ...
                                        plotdat(M2idx,:)]);
                        
                        toplotdat = basedata - meanM;
                        
                        plt=plot(temps,handles.scale*i+toplotdat);

                    case    'BIPOLAR'
                        
                        if handles.crc_types(index2)>0
                            [dumb1,index3,dumb2] = ...
                                intersect(upper(chanlabels(handles.Dmeg)), ...
                                      upper(handles.names(handles.crc_types(index2))));
                        else
                            index3 = [];
                        end
                        if ~isempty(index3)
                            bipolar = plotdat(index(i),:) - ...
                                        plotdat(index3,:);
                            plt=plot(temps,handles.scale*i+bipolar,'Color',[0.2 0.9 0.5]);
                        else
                            plt=plot(temps,handles.scale*i+plotdat(index(i),:), ...
                                        'Color',[0.2 0.5 0.9]);
                        end

                    otherwise

                        [dumb1,index3,dumb2] = ...
                            intersect(upper(chanlabels(handles.Dmeg)),selchan);
                        basedata = plotdat(index(i),:);
                        toplotdat = plotdat(index(i),:) - ...
                                    plotdat(index3,:);
                        
                        plt=plot(temps,handles.scale*i+toplotdat);
                        
                end
                
                if get(handles.EMGfil,'Value')
                    filterlowhigh(plt,i,handles,handles.filter.EMG)
                end

            elseif  strfind(testcond{:},'EOG')

                contents = get(handles.EOGpopmenu,'String');
                selchan=upper(contents{get(handles.EOGpopmenu,'Value')});
                
                switch  selchan
                
                    case    'REF1'

                        plt=plot(temps,handles.scale*i + ...
                                    plotdat(index(i),:));

                    case    'MEAN OF REF'
                        
                        basedata = plotdat(index(i),:);
                        ref2idx = find(strcmp(chanlabels(handles.Dmeg),'REF2'));
                        scnddata = plotdat(index(i),:) - ...
                                    plotdat(ref2idx,:);
                        
                        toplotdat = mean([basedata ; scnddata]);
                        
                        plt=plot(temps,handles.scale*i+toplotdat);
            
                    case    'M1-M2'

                        basedata = plotdat(index(i),:);
                        M1idx = find(strcmp(chanlabels(handles.Dmeg),'M1'));
                        M2idx = find(strcmp(chanlabels(handles.Dmeg),'M2'));
                        meanM = mean([plotdat(M1idx,:) ; ...
                                        plotdat(M2idx,:)]);

                        toplotdat = basedata - meanM;

                        plt=plot(temps,handles.scale*i+toplotdat);

                        
                    case    'BIPOLAR'

                        if handles.crc_types(index2)>0
                            [dumb1,index3,dumb2] = ...
                                intersect(upper(chanlabels(handles.Dmeg)), ...
                                      upper(handles.names(handles.crc_types(index2))));
                        else
                            index3 = [];
                        end
                        if ~isempty(index3)
                            bipolar = plotdat(index(i),:) - ...
                                        plotdat(index3,:);
                            plt=plot(temps,handles.scale*i+bipolar,'Color',[0.2 0.9 0.5]);
                        else
                            plt=plot(temps,handles.scale*i+plotdat(index(i),:), ...
                                        'Color',[0.2 0.5 0.9]);
                        end

                    otherwise

                        [dumb1,index3,dumb2] = ...
                            intersect(upper(chanlabels(handles.Dmeg)),selchan);
                        basedata = plotdat(index(i),:);
                        toplotdat = plotdat(index(i),:) - ...
                                    plotdat(index3,:);
                        
                        plt=plot(temps,handles.scale*i+toplotdat);
                        
                end
                if get(handles.EOGfil,'Value')
                    filterlowhigh(plt,i,handles,handles.filter.EOG)
                end

            elseif  strcmp(upper(chantype(handles.Dmeg,index(i))),'EEG')                

                contents = get(handles.EEGpopmenu,'String');
                selchan=upper(contents{get(handles.EEGpopmenu,'Value')});
                
                switch  selchan
                
                    case    'REF1'

                        plt=plot(temps,handles.scale*i+plotdat(index(i),:));

                    case    'MEAN OF REF'
                        
                        basedata = plotdat(index(i),:);
                        ref2idx = find(strcmp(chanlabels(handles.Dmeg),'REF2'));
                        scnddata = plotdat(index(i),:) - ...
                                    plotdat(ref2idx,:);
                        
                        toplotdat = mean([basedata ; scnddata]);
                        
                        plt=plot(temps,handles.scale*i+toplotdat);

                    case    'M1-M2'
                        
                        basedata = plotdat(index(i),:);
                        M1idx = find(strcmp(chanlabels(handles.Dmeg),'M1'));
                        M2idx = find(strcmp(chanlabels(handles.Dmeg),'M2'));
                        meanM = mean([plotdat(M1idx,:) ; ...
                                      plotdat(M2idx,:)]);
                        
                        toplotdat = basedata - meanM;
                        
                        plt=plot(temps,handles.scale*i+toplotdat);

                        
                    case    'BIPOLAR'

                        if handles.crc_types(index2)>0
                            [dumb1,index3,dumb2] = ...
                                intersect(upper(chanlabels(handles.Dmeg)), ...
                                        upper(handles.names(handles.crc_types(index2))));
                        else
                            index3 = [];
                        end
                        if ~isempty(index3)
                            bipolar = plotdat(index(i),:) - ...
                                        plotdat(index3,:);
                            plt=plot(temps,handles.scale*i+bipolar,'Color',[0.2 0.9 0.5]);
                        else
                            plt=plot(temps,handles.scale*i+plotdat(index(i),:), ...
                                        'Color',[0.2 0.5 0.9]);                            
                        end

                    otherwise

                        [dumb1,index3,dumb2] = ...
                            intersect(upper(chanlabels(handles.Dmeg)),selchan);
                        basedata = plotdat(index(i),:);
                        toplotdat = plotdat(index(i),:) - ...
                                    plotdat(index3,:);
                        
                        plt=plot(temps,handles.scale*i+toplotdat);
                        
                end
                if get(handles.checkfil,'Value')
                    filterlowhigh(plt,i,handles)
                end
                
            elseif strfind(testcond{:},'ECG')|strfind(testcond{:},'EKG')

                contents = get(handles.otherpopmenu,'String');
                selchan=upper(contents{get(handles.otherpopmenu,'Value')});
                
                switch  selchan
                
                    case    'REF1'
                            
                        toplotdat = plotdat(index(i),:);
                        
                        
                        range=max(toplotdat)-min(toplotdat);
                        if range == 0;range=1;end;
                        plt=plot(temps,handles.scale*i+handles.scale*toplotdat/range,'r');

                    case    'MEAN OF REF'
                        
                        basedata = plotdat(index(i),:);
                        ref2idx = find(strcmp(chanlabels(handles.Dmeg),'REF2'));
                        scnddata = plotdat(index(i),:) - ...
                                    plotdat(ref2idx,:);
                        
                        toplotdat = mean([basedata ; scnddata]);

                        range=max(toplotdat)-min(toplotdat);
                        if range == 0;range=1;end;
                        plt=plot(temps,handles.scale*i+handles.scale*toplotdat/range,'r');

                    case    'M1-M2'
                        
                        basedata = plotdat(index(i),:);
                        M1idx = find(strcmp(chanlabels(handles.Dmeg),'M1'));
                        M2idx = find(strcmp(chanlabels(handles.Dmeg),'M2'));
                        meanM = mean([plotdat(M1idx,:) ; ...
                                        plotdat(M2idx,:)]);
                        
                        toplotdat = basedata - meanM;
                        
                        range=max(toplotdat)-min(toplotdat);
                        if range == 0;range=1;end;
                        plt=plot(temps,handles.scale*i+handles.scale*toplotdat/range,'r');
                        
                    case    'BIPOLAR'

                        if handles.crc_types(index2)>0
                            [dumb1,index3,dumb2] = ...
                                intersect(upper(chanlabels(handles.Dmeg)), ...
                                        upper(handles.names(handles.crc_types(index2))));
                        else
                            index3 = [];
                        end
                        if ~isempty(index3)
                            bipolar = plotdat(index(i),:) - ...
                                        plotdat(index3,:);
                        	plt=plot(temps,handles.scale*i+bipolar,'Color',[0.2 0.9 0.5]);
                        else
                            plt=plot(temps,handles.scale*i+handles.scale*toplotdat/range, ...
                                        'Color',[0.2 0.5 0.9]);
                        end

                    otherwise

                        [dumb1,index3,dumb2] = ...
                            intersect(upper(chanlabels(handles.Dmeg)),selchan);
                        basedata = plotdat(index(i),:);
                        toplotdat = plotdat(index(i),:) - ...
                                    plotdat(index3,:);
                        
                        range=max(toplotdat)-min(toplotdat);
                        if range == 0;range=1;end;
                        plt=plot(temps,handles.scale*i+handles.scale*toplotdat/range,'r');
                        
                end
                if get(handles.checkfil,'Value')
                    filterlowhigh(plt,i,handles)
                end

            else

                contents = get(handles.EEGpopmenu,'String');
                selchan=upper(contents{get(handles.EEGpopmenu,'Value')});
                
                switch  selchan
                
                    case    'REF1'

                        plt=plot(temps,handles.scale*i+plotdat(index(i),:));

                    case    'MEAN OF REF'
                        
                        basedata = plotdat(index(i),:);
                        ref2idx = find(strcmp(chanlabels(handles.Dmeg),'REF2'));
                        scnddata = plotdat(index(i),:) - ...
                                    plotdat(ref2idx,:);
                        
                        toplotdat = mean([basedata ; scnddata]);
                        
                        plt=plot(temps,handles.scale*i+toplotdat);

                    case    'M1-M2'
                        
                        basedata = plotdat(index(i),:);
                        M1idx = find(strcmp(chanlabels(handles.Dmeg),'M1'));
                        M2idx = find(strcmp(chanlabels(handles.Dmeg),'M2'));
                        meanM = mean([plotdat(M1idx,:) ; ...
                                      plotdat(M2idx,:)]);
                        
                        toplotdat = basedata - meanM;
                        
                        plt=plot(temps,handles.scale*i+toplotdat);

                        
                    case    'BIPOLAR'

                        if handles.crc_types(index2)>0
                            [dumb1,index3,dumb2] = ...
                                intersect(upper(chanlabels(handles.Dmeg)), ...
                                        upper(handles.names(handles.crc_types(index2))));
                        else
                            index3 = [];
                        end
                        if ~isempty(index3)
                            bipolar = plotdat(index(i),:) - ...
                                        plotdat(index3,:);
                            plt=plot(temps,handles.scale*i+bipolar,'Color',[0.2 0.9 0.5]);
                        else
                            plt=plot(temps,handles.scale*i+plotdat(index(i),:), ...
                                        'Color',[0.2 0.5 0.9]);                            
                        end

                    otherwise

                        [dumb1,index3,dumb2] = ...
                            intersect(upper(chanlabels(handles.Dmeg)),selchan);
                        basedata = plotdat(index(i),:);
                        toplotdat = plotdat(index(i),:) - ...
                                    plotdat(index3,:);
                        
                        plt=plot(temps,handles.scale*i+toplotdat);
                        
                end
                if get(handles.checkfil,'Value')
                    filterlowhigh(plt,i,handles)
                end
                            
            end

            
            
end
            
ylim([0 handles.scale*(i+1)])

lg=length(index);
set(handles.axes1,'YTick', ...
    [handles.scale/2:handles.scale/2:lg*handles.scale+handles.scale/2]);

ylabels=[num2str(round(handles.scale/2))];
for j=1:length(index)
    ylabels=[ylabels chanlabels(handles.Dmeg,index(j))];
    ylabels=[ylabels num2str(round(handles.scale/2))];
end
set(handles.axes1,'YTickLabel',ylabels);


xlim([slidval slidval+handles.winsize])
xtick = get(handles.axes1,'XTick');
if isfield(handles.Dmeg,'info')
    if isfield(handles.Dmeg.info,'hour')
        xtick = mod(xtick + handles.offset,24*60^2);
    end
end
[time string] = crc_time_converts(xtick);
set(handles.axes1,'XTickLabel',string)

set(handles.currenttime,'String', ...
    num2str(round(tdeb/fsample(handles.Dmeg)+handles.winsize/2)));

grid on

% Display trigger
ev = events(handles.Dmeg);
if isempty(ev)
    ev=struct('time', -500,'value', -2000);
end

try
    [int indextoshow indextrig] = intersect(toshow,round([ev(:).time]*fsample(handles.Dmeg)));
catch
    [int indextoshow indextrig] = intersect(toshow,round([ev{:}.time]*fsample(handles.Dmeg)));
end

%Check if some events happens at the same time
if ~isempty(indextrig)
    if indextrig(1)~=1
        if and(ev(indextrig(1)).time == ev(indextrig(1)-1).time,length(indextrig)~=length(indextrig(1):indextrig(end)))

            trou=diff(indextrig);
            trouve=find(trou==2)+1;

            int = sort([int(1) int int(trouve)]);
            %    indextoshow = sort([indextoshow indextoshow(trouve)]);
            indextrig = sort([(indextrig(1)-1) indextrig (indextrig(trouve)-1)]);

        elseif ev(indextrig(1)).time == ev(indextrig(1)-1).time
            int = sort([int(1) int]);
            indextrig = sort([(indextrig(1)-1) indextrig]);

        end
    end

    if length(indextrig)~=length(indextrig(1):indextrig(end))

        trou=diff(indextrig);
        trouve=find(trou==2)+1;

        int = sort([int int(trouve)]);
        indextrig = sort([indextrig (indextrig(trouve)-1)]);
    end

end
int = int/fsample(handles.Dmeg);
Nev_dis = length(int);

if Nev_dis % do all this if there are triggers to be displayed !
    try
        tmp_val = ev(indextrig(1)).value;
    catch
        tmp_val = ev{indextrig(1)}.value;
    end
    if ischar(tmp_val)
        use_numv = 0;
    else
        use_numv = 0;
    end
    
    try
        if use_numv
            etype = [ev(indextrig).value];
        else
            etype = strvcat(num2str([ev(indextrig).value]'));
        end
    catch
        if use_numv
            etype = [ev{indextrig}.value];
        else
            etype = strvcat(num2str([ev{indextrig}.value]'));
        end
    end

    plot(int,0.5*ones(1,length(int))*handles.scale/50*NbreChandisp, ...
                'k^','LineWidth',2.5)

    for jj = 1:Nev_dis
        if use_numv
            if etype(jj)<0
                msg = ['R' num2str(abs(etype(jj)))];
            elseif etype(jj)>0 && etype(jj)<1000
                msg = ['S' num2str(etype(jj))];
            else
                msg = ['K' num2str(etype(jj)-1000)];
            end
        else
            msg = etype(jj,:);
        end
        lgmsg = length(msg);
        text(int(jj)-(0.4*lgmsg/(lgmsg+1))*handles.winsize/20, ...
                3*handles.scale/50*NbreChandisp,msg);
    end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function filterlowhigh(plt,ii,handles,frqcut)

if nargin<4
    flc = handles.filter.other(1)/fsample(handles.Dmeg);
    fhc = handles.filter.other(2)/fsample(handles.Dmeg);
else
	flc = frqcut(1)/fsample(handles.Dmeg);
    fhc = frqcut(2)/fsample(handles.Dmeg);
end

k = .7; % cut-off value

alphal = (1-k*cos(2*pi*flc)-sqrt(2*k*(1-cos(2*pi*flc))-k^2*sin(2*pi*flc)^2))/(1-k);
alphah = (1-k*cos(2*pi*fhc)-sqrt(2*k*(1-cos(2*pi*fhc))-k^2*sin(2*pi*fhc)^2))/(1-k);

    
X = get(plt,'YData');
   
% Apply low pass filter
Y = filtfilt(1-alphal,[1 -alphal],X);
	
% Apply high pass filter
Y = filtfilt(1-alphah,[1 -alphah],X-Y);

set(plt,'YData',Y+(ii)*handles.scale - mean(Y))%,'Color',Col{ii})

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Y=filterforspect(handles,X,frqcut)

	flc = frqcut(1)/fsample(handles.Dmeg);
    fhc = frqcut(2)/fsample(handles.Dmeg);

k = .7; % cut-off value

alphal = (1-k*cos(2*pi*flc)-sqrt(2*k*(1-cos(2*pi*flc))-k^2*sin(2*pi*flc)^2))/(1-k);
alphah = (1-k*cos(2*pi*fhc)-sqrt(2*k*(1-cos(2*pi*fhc))-k^2*sin(2*pi*fhc)^2))/(1-k);
 
% Apply low pass filter
Y = filtfilt(1-alphal,[1 -alphal],X);
	
% Apply high pass filter
Y = filtfilt(1-alphah,[1 -alphah],X-Y);

return
%%%%%%%%%%%%%%%%%%%%%%%

function cleargraph(figure)

A=get(figure,'Children');

idx=find(strcmp(get(A,'Type'),'axes')==1);

delete(get(A(idx),'Children'))



