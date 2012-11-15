function varargout = crc_par_selchan(varargin)
% DIS_SELCHAN M-file for dis_selchan.fig
%      DIS_SELCHAN, by itself, creates a new DIS_SELCHAN or raises the existing
%      singleton*.
%
%      H = DIS_SELCHAN returns the handle to a new DIS_SELCHAN or the handle to
%      the existing singleton*.
%
%      DIS_SELCHAN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DIS_SELCHAN.M with the given input arguments.
%
%      DIS_SELCHAN('Property','Value',...) creates a new DIS_SELCHAN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dis_selchan_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dis_selchan_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help dis_selchan

% Last Modified by GUIDE v2.5 21-Jan-2009 11:52:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @crc_par_selchan_OpeningFcn, ...
    'gui_OutputFcn',  @crc_par_selchan_OutputFcn, ...
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
end

% --- Executes just before dis_selchan is made visible.
function crc_par_selchan_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dis_selchan (see VARARGIN)

% Choose default command line output for dis_selchan
handles.output = hObject;

set(0,'CurrentFigure',handles.figure1);

if length(varargin)<1
	Dmeg = crc_eeg_load;
else
    Dmeg = crc_eeg_load(deblank(varargin{1}{1}));
end
file = fullfile(Dmeg.path,Dmeg.fname);
set(handles.Filename,'String',Dmeg.fname);

handles.file = file;
handles.Dmeg = Dmeg;

MEEGchan = meegchannels(Dmeg);
BADchan = badchannels(Dmeg);

select = setdiff(MEEGchan,BADchan);
deselect = setdiff(1:Dmeg.nchannels,select);
set(handles.Deselect,'String',upper(chanlabels(handles.Dmeg,select)));
set(handles.Select,'String',upper(chanlabels(handles.Dmeg,deselect)));

guidata(hObject, handles);

end
% Update handles structure

% UIWAIT makes dis_selchan wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = crc_par_selchan_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

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
end

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

end
% Update handles structure
guidata(hObject, handles);

end
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
end
% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[dumb1,index,dumb2]=intersect(upper(chanlabels(handles.Dmeg)),get(handles.Select,'String'));
handles.output = index;

set(handles.figure1,'UserData',1);
guidata(hObject, handles);

end