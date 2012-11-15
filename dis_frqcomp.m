function varargout = dis_frqcomp(varargin)
% DIS_FRQCOMP M-file for dis_frqcomp.fig
%      DIS_FRQCOMP, by itself, creates a new DIS_FRQCOMP or raises the existing
%      singleton*.
%
%      H = DIS_FRQCOMP returns the handle to a new DIS_FRQCOMP or the handle to
%      the existing singleton*.
%
%      DIS_FRQCOMP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DIS_FRQCOMP.M with the given input arguments.
%
%      DIS_FRQCOMP('Property','Value',...) creates a new DIS_FRQCOMP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dis_frqcomp_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dis_frqcomp_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help dis_frqcomp

% Last Modified by GUIDE v2.5 19-Sep-2008 11:42:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dis_frqcomp_OpeningFcn, ...
                   'gui_OutputFcn',  @dis_frqcomp_OutputFcn, ...
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


% --- Executes just before dis_frqcomp is made visible.
function dis_frqcomp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dis_frqcomp (see VARARGIN)

% Choose default command line output for dis_frqcomp
handles.output = hObject;
set(0,'CurrentFigure',handles.figure1);


try
    handles.file = varargin{5}.file;
    D = crc_eeg_load(handles.file);
catch
    D = crc_eeg_load;
    handles.file = fullfile(D.path,D.fname);
end

handles.D = D;

if isfield(D, 'CRC')
    if isfield(D.CRC,'score')       
        set(handles.scorerpopmenu,'String',D.CRC.score(2,:))
        set(handles.scorerpopmenu,'Value',1)
    else
    	set(handles.scorerpopmenu,'Visible','off')
        set(handles.text6,'Visible','off')
    end
else
   	set(handles.scorerpopmenu,'Visible','off')
    set(handles.text6,'Visible','off')
end

if ismember('REF2',upper(chanlabels(D)));
    Channames=[{'REF'} upper(chanlabels(D)) {'MEAN OF REF'}];
else
    Channames=[{'REF'}  upper(chanlabels(D))];   
end

if and(ismember('M1',upper(chanlabels(D))),ismember('M2',upper(chanlabels(D))));
    Channames=[Channames {'M1-M2'}];
end

set(handles.popref,'String',Channames);


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dis_frqcomp wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = dis_frqcomp_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function Fmax_Callback(hObject, eventdata, handles)
% hObject    handle to Fmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Fmax as text
%        str2double(get(hObject,'String')) returns contents of Fmax as a double


% --- Executes during object creation, after setting all properties.
function Fmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Fmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Fmin_Callback(hObject, eventdata, handles)
% hObject    handle to Fmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Fmin as text
%        str2double(get(hObject,'String')) returns contents of Fmin as a double


% --- Executes during object creation, after setting all properties.
function Fmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Fmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Duration_Callback(hObject, eventdata, handles)
% hObject    handle to Duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Duration as text
%        str2double(get(hObject,'String')) returns contents of Duration as a double


% --- Executes during object creation, after setting all properties.
function Duration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Step_Callback(hObject, eventdata, handles)
% hObject    handle to Step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Step as text
%        str2double(get(hObject,'String')) returns contents of Step as a double


% --- Executes during object creation, after setting all properties.
function Step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Cmp_frq.
function varargout=Cmp_frq_Callback(hObject, eventdata, handles)
% hObject    handle to Cmp_frq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'file')
    args.file=handles.file; 
end

if isfield(handles,'D')
    args.D=handles.D; 
end

if not(isnan(str2double(get(handles.Duration,'String'))))
    args.dur=str2double(get(handles.Duration,'String'));
end

if  not(isnan(str2double(get(handles.Step,'String'))))
    args.step=str2double(get(handles.Step,'String'));
end

if  not(isnan(str2double(get(handles.Fmax,'String'))))
    args.fmax=str2double(get(handles.Fmax,'String'));
end

if  not(isnan(str2double(get(handles.Fmin,'String'))))
    args.fmin=str2double(get(handles.Fmin,'String'));
end

if strcmp(get(handles.scorerpopmenu,'Visible'),'on')
    args.scorer = get(handles.scorerpopmenu,'Value');
end

args.ref = get(handles.popref,'Value');
D = crc_spectcompute(args);
varargout{2} = D;
delete(handles.figure1)
h = timerfind;
if ~isempty(h)
    stop(h);
end



% --- Executes on selection change in scorerpopmenu.
function scorerpopmenu_Callback(hObject, eventdata, handles)
% hObject    handle to scorerpopmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns scorerpopmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from scorerpopmenu


% --- Executes during object creation, after setting all properties.
function scorerpopmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scorerpopmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in popref.
function popref_Callback(hObject, eventdata, handles)
% hObject    handle to popref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popref contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popref


% --- Executes during object creation, after setting all properties.
function popref_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


