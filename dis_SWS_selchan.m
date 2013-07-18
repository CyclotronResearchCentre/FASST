function varargout = dis_SWS_selchan(varargin)
% dis_SWS_selchan M-file for dis_SWS_selchan.fig
%
% GUI to select electrodes in order to build ROI on which the signal will
% be averaged before SWS detection
%
%      dis_SWS_selchan, by itself, creates a new dis_SWS_selchan or raises the existing
%      singleton*.
%
%      H = dis_SWS_selchan returns the handle to a new dis_SWS_selchan or the handle to
%      the existing singleton*.
%
%      dis_SWS_selchan('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in dis_SWS_selchan.M with the given input arguments.
%
%      dis_SWS_selchan('Property','Value',...) creates a new dis_SWS_selchan or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dis_SWS_selchan_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dis_SWS_selchan_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%__________________________________________________________________

% Written by J. Schrouff, Y. Leclercq & C. Phillips, 2008.
% Cyclotron Research Centre, University of Liège, Belgium
% $Id$

% Last Modified by GUIDE v2.5 02-Oct-2009 11:56:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dis_SWS_selchan_OpeningFcn, ...
                   'gui_OutputFcn',  @dis_SWS_selchan_OutputFcn, ...
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


% --- Executes just before dis_SWS_selchan is made visible.
function dis_SWS_selchan_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dis_SWS_selchan (see VARARGIN)
% Choose default command line output for dis_SWS_selchan
handles.output = hObject;

load CRC_electrodes.mat;
% handles.names=names;
for i=1:size(names,2)
    handles.names{i}=upper([names{i}]);
end
% handles.pos=pos;
handles.pos=pos';
% handles.types=types;
handles.crc_types=crc_types;

if length(varargin)<1
    Dmeg = crc_eeg_load;
    file = fullfile(D.path,D.fname);
    
    handles.file           = file;
    handles.spmstruct      = struct(Dmeg);
    handles.Dmeg           = Dmeg;
    set(handles.Deselect,'String',upper(chanlabels(handles.Dmeg)));
else
    handles.file      = varargin{1}.file;
    handles.spmstruct = varargin{1}.spmstruct;
    handles.Dmeg      = varargin{1}.Dmeg;
    set(handles.Deselect,'String',upper(chanlabels(handles.Dmeg)));
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dis_SWS_selchan wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = dis_SWS_selchan_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
% varargout{2}=sel_ROI;

% --- Executes on selection change in Select.
function Select_Callback(hObject, eventdata, handles)
% hObject    handle to Select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Select contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Select
contents = get(hObject,'String');
if length(contents)==0
else
%Remove the "activated" item from the list "Available Channels"
    [dumb1,dumb2,index]=intersect(contents{get(hObject,'Value')},contents);
    temp=[contents(1:index-1) ; contents(index+1:length(contents))];
    set(handles.Select,'String',temp);
    
    [dumb1,dumb2,index2]=intersect(temp,handles.names);
    
    idxred=index2(find(handles.crc_types(index2)<-1));
    idxblue=index2(find(handles.crc_types(index2)>-2));

    xred=handles.pos(1,idxred);
    yred=handles.pos(2,idxred);
    
    xblu=handles.pos(1,idxblue);
    yblu=handles.pos(2,idxblue);
    
    cleargraph(handles)
    hold on
    plot(xred,yred,'r+'), plot(xblu,yblu,'b+')
    hold off
    if and(length(xblu)==0,length(xred)==0)
        cleargraph(handles)
    end
    
   
    set(handles.Localizer,'XTick',[]);
    set(handles.Localizer,'YTick',[]);

    %Add the "activated" in the list "Selected Channels"
    if length(get(handles.Deselect,'String'))==0
        temp={contents{get(hObject,'Value')}};
    else
        temp=[contents{get(hObject,'Value')} ; get(handles.Deselect,'String')];
    end
    set(handles.Deselect,'String',temp);

%Prevent crashing if the first/last item of the list is selected.
set(handles.Select,'Value',max(index-1,1));
set(handles.Deselect,'Value',1);
end

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Select_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Deselect.
function Deselect_Callback(hObject, eventdata, handles)
% hObject    handle to Deselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = get(hObject,'String');
if length(contents)==0
else
%Remove the "activated" item from the list "Available Channels"
    [dumb1,dumb2,index]=intersect(contents{get(hObject,'Value')},contents);
    temp=[contents(1:index-1) ; contents(index+1:length(contents))];
    set(handles.Deselect,'String',temp);

%Add the "activated" in the list "Selected Channels"
    if length(get(handles.Select,'String'))==0
        temp={contents{get(hObject,'Value')}};
    else
        temp=[contents{get(hObject,'Value')} ; get(handles.Select,'String')];
    end
    set(handles.Select,'String',temp);

%Prevent crashing if the first/last item of the list is selected.
    set(handles.Deselect,'Value',max(index-1,1));
    set(handles.Select,'Value',1);

    [dumb1,dumb2,index]=intersect(temp,handles.names);

    idxred=index(find(handles.crc_types(index)<-1));
    idxblue=index(find(handles.crc_types(index)>-2));

    xred=handles.pos(1,idxred);
    yred=handles.pos(2,idxred);
    
    xblu=handles.pos(1,idxblue);
    yblu=handles.pos(2,idxblue);
    
    cleargraph(handles)
    hold on
    plot(xred,yred,'r+'), plot(xblu,yblu,'b+')
    hold off
    if and(length(xblu)==0,length(xred)==0)
        cleargraph(handles)
    end

    xlim([0 1])
    ylim([0 1])

    set(handles.Localizer,'XTick',[]);
    set(handles.Localizer,'YTick',[]);

end
% Update handles structure
guidata(hObject, handles);


% Hints: contents = get(hObject,'String') returns Deselect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Deselect


% --- Executes during object creation, after setting all properties.
function Deselect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Deselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in pushbuttonok.
function pushbuttonok_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 [dumb1,index,dumb2]=intersect(upper(chanlabels(handles.Dmeg)),get(handles.Select,'String'));
try
    flags.index = sortch(handles.spmstruct,index);
catch
    flags.index = fliplr(sort(index))  ;  
end
flags.spmstruct = handles.spmstruct;
flags.file = handles.file;
handles.index = flags.index;


% Update handles structure
guidata(hObject, handles);
global gindex 
gindex = handles.index;
delete(handles.figure1)


% --- Executes on button press in desall.
function desall_Callback(hObject, eventdata, handles)
% hObject    handle to desall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.Deselect,'String',upper(chanlabels(handles.Dmeg)));
set(handles.Select,'String',cell(0));

clmo all
xlim([0 1])
ylim([0 1])
set(handles.Localizer,'XTick',[]);
set(handles.Localizer,'YTick',[]);
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbuttonload.
function pushbuttonload_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonload (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

indfile=spm_select(1,'mat','Select ROI configuration','',pwd);
load(indfile)
index=ind;
set(handles.Select,'String',upper(chanlabels(handles.Dmeg,index)));
diff=setdiff(upper(chanlabels(handles.Dmeg)),upper(chanlabels(handles.Dmeg,index)));
set(handles.Deselect,'String',diff);

[dumb1,dumb2,index2]=intersect(upper(chanlabels(handles.Dmeg,index)),handles.names);

idxred=index2(find(handles.crc_types(index2)<-1));
idxblue=index2(find(handles.crc_types(index2)>-2));

xred=handles.pos(1,idxred);
yred=handles.pos(2,idxred);

xblu=handles.pos(1,idxblue);
yblu=handles.pos(2,idxblue);

hold on
plot(xred,yred,'r+'), plot(xblu,yblu,'b+')
hold off
if and(length(xblu)==0,length(xred)==0)
    clmo all
end


xlim([0 1])
ylim([0 1])

% Update handles structure
guidata(hObject, handles);
    

% --- Executes on button press in pushbuttonsave.
function pushbuttonsave_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonsave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[dumb1,index,dumb2]=intersect(upper(chanlabels(handles.Dmeg)),get(handles.Select,'String'));
try
    handles.index = sortch(handles.spmstruct,index);
catch
    handles.index = fliplr(sort(index))  ;  
end
ind=handles.index;
uisave('ind')
% Update handles structure
guidata(hObject, handles);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SORT CHANNEL FX
function idx = sortch(D,index)

% Find the EOG channels
A = strfind(upper(chanlabels(D,index)),'EOG');
EOGchan=[];
for i=1:length(A);
	if A{i}>0
    	EOGchan=[EOGchan index(i)];
	end
end

% Find the ECG channels
A=strfind(upper(chanlabels(D,index)),'ECG');
ECGchan=[];
for i=1:length(A);
	if A{i}>0
        ECGchan=[ECGchan index(i)];
	end
end

% Find the EMG channels
A=strfind(upper(chanlabels(D,index)),'EMG');
EMGchan=[];
for i=1:length(A);
	if A{i}>0
        EMGchan=[EMGchan index(i)];
	end
end

A=strfind(upper(chantype(D)),'EEG');
EEGchan=[];
for i=1:length(A);
	if A{i}>0
        EEGchan=[EEGchan index(i)];
	end
end


allbad = [EOGchan ECGchan EMGchan EEGchan];
other = setdiff(index,allbad);

otherknown=[];
othernotknown=[];
for ff = other
    if d.channels.order(ff)==0
        othernotknown = [othernotknown ff];
    else
        otherknown = [otherknown ff];
    end
end
other = [otherknown othernotknown];


allbad=[EOGchan ECGchan EMGchan other];
eeg=setdiff(index,allbad);

AFrontal = intersect(find(strncmp(chanlabels(handles.Dmeg),'AF',2) ==1),eeg);
Frontal = intersect(find(strncmp(chanlabels(handles.Dmeg),'F',1) ==1),eeg);
Coronal = intersect(find(strncmp(chanlabels(handles.Dmeg),'C',1) ==1),eeg);
Temporal = intersect(find(strncmp(chanlabels(handles.Dmeg),'T',1) ==1),eeg);
Parietal= intersect(find(strncmp(chanlabels(handles.Dmeg),'P',1) ==1),eeg);
Occipital= intersect(find(strncmp(chanlabels(handles.Dmeg),'O',1) ==1),eeg);

neweeg = [Occipital Parietal Temporal Coronal Frontal AFrontal];
eeg2=setdiff(eeg,neweeg);

eeg = [eeg2 neweeg];

idx=[otherknown othernotknown ECGchan EMGchan EOGchan eeg]; 
%%%%%%%%%%%%%%%
function cleargraph(handles)

A=get(handles.figure1,'Children');

idx=find(strcmp(get(A,'Type'),'axes')==1);

delete(get(A(idx),'Children'))

