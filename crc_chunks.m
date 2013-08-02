function varargout = crc_chunks(varargin)
% CRC_CHUNKS M-file for crc_chunks.fig
%      CRC_CHUNKS, by itself, creates a new CRC_CHUNKS or raises the existing
%      singleton*.
%
%      H = CRC_CHUNKS returns the handle to a new CRC_CHUNKS or the handle to
%      the existing singleton*.
%
%      CRC_CHUNKS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CRC_CHUNKS.M with the given input arguments.
%
%      CRC_CHUNKS('Property','Value',...) creates a new CRC_CHUNKS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before crc_chunks_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to crc_chunks_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%__________________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by Y. Leclercq & C. Phillips, 2008.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

% Edit the above text to modify the response to help crc_chunks

% Last Modified by GUIDE v2.5 22-Apr-2009 15:56:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @crc_chunks_OpeningFcn, ...
    'gui_OutputFcn',  @crc_chunks_OutputFcn, ...
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


% --- Executes just before crc_chunks is made visible.
function crc_chunks_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to crc_chunks (see VARARGIN)

% Choose default command line output for crc_chunks
handles.output = hObject;

set(0,'CurrentFigure',handles.figure1);

% Handling loading of data file.
Dmeg = crc_eeg_load;
D = struct(Dmeg);

handles.file = fullfile(Dmeg.path,Dmeg.fname);
handles.Dmeg = Dmeg;

if ~isempty(D.trials.events) 
    if strcmp(D.trials.events(1).value,'start')
        handles.markerlist=unique([D.trials.events(2:end).value]);
    else
        handles.markerlist=unique([D.trials.events(:).value]);
    end
elseif ~isempty(D.trials.events)
    handles.markerlist=unique([D.trials.events(:).value]);
else
    handles.markerlist=[];
end

if isfield(Dmeg,'info')
    if isfield(Dmeg.info,'hour')
        Str = [' ' num2str(Dmeg.info.hour(1)) 'h ' num2str(Dmeg.info.hour(2)) 'm ' num2str(Dmeg.info.hour(3)) 's' ];
    end
end

handles.Begrtime = [0 0 0];
handles.Endrtime = crc_time_converts(nsamples(Dmeg)/fsample(Dmeg));

% If program can find a start time, then working in clock time is possible
% and the corresponding variable are setted up.

if exist('Str','var')

    set(handles.clocktime,'String',Str)
    handles.clocktime=1;
    handles.Begctime = Dmeg.info.hour;

    handles.Endctime = crc_time_converts(Dmeg.info.hour(1)*60^2 + ...
        Dmeg.info.hour(2)*60 + ...
        Dmeg.info.hour(3) + ...
        nsamples(Dmeg)/fsample(Dmeg));

    handles.Endctime(1) = mod(handles.Endctime(1),24);
else
    handles.clocktime=0;
end

update(hObject,handles)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes crc_chunks wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = crc_chunks_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in Begstyle.
function Begstyle_Callback(hObject, eventdata, handles)
% hObject    handle to Begstyle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Begstyle contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Begstyle
update(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Begstyle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Begstyle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Endstyle.
function Endstyle_Callback(hObject, eventdata, handles)
% hObject    handle to Endstyle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Endstyle contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Endstyle
update(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Endstyle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Endstyle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Begh_Callback(hObject, eventdata, handles)
% hObject    handle to Begh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Begh as text
%        str2double(get(hObject,'String')) returns contents of Begh as a double


% --- Executes during object creation, after setting all properties.
function Begh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Begh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Begm_Callback(hObject, eventdata, handles)
% hObject    handle to Begm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Begm as text
%        str2double(get(hObject,'String')) returns contents of Begm as a double


% --- Executes during object creation, after setting all properties.
function Begm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Begm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Begs_Callback(hObject, eventdata, handles)
% hObject    handle to Begs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Begs as text
%        str2double(get(hObject,'String')) returns contents of Begs as a double


% --- Executes during object creation, after setting all properties.
function Begs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Begs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Endh_Callback(hObject, eventdata, handles)
% hObject    handle to Endh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Endh as text
%        str2double(get(hObject,'String')) returns contents of Endh as a double


% --- Executes during object creation, after setting all properties.
function Endh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Endh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Endm_Callback(hObject, eventdata, handles)
% hObject    handle to Endm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Endm as text
%        str2double(get(hObject,'String')) returns contents of Endm as a double


% --- Executes during object creation, after setting all properties.
function Endm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Endm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ends_Callback(hObject, eventdata, handles)
% hObject    handle to Ends (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ends as text
%        str2double(get(hObject,'String')) returns contents of Ends as a double


% --- Executes during object creation, after setting all properties.
function Ends_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ends (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Begtimepop.
function Begtimepop_Callback(hObject, eventdata, handles)
% hObject    handle to Begtimepop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update(hObject,handles)
% Hints: contents = get(hObject,'String') returns Begtimepop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Begtimepop


% --- Executes during object creation, after setting all properties.
function Begtimepop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Begtimepop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Endtimepop.
function Endtimepop_Callback(hObject, eventdata, handles)
% hObject    handle to Endtimepop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update(hObject,handles)
% Hints: contents = get(hObject,'String') returns Endtimepop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Endtimepop


% --- Executes during object creation, after setting all properties.
function Endtimepop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Endtimepop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Begtypemarker.
function Begtypemarker_Callback(hObject, eventdata, handles)
% hObject    handle to Begtypemarker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Begtypemarker contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Begtypemarker
update(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Begtypemarker_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Begtypemarker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Begnumark.
function Begnumark_Callback(hObject, eventdata, handles)
% hObject    handle to Begnumark (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Begnumark contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Begnumark


% --- Executes during object creation, after setting all properties.
function Begnumark_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Begnumark (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Endtypemarker.
function Endtypemarker_Callback(hObject, eventdata, handles)
% hObject    handle to Endtypemarker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Endtypemarker contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Endtypemarker
update(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Endtypemarker_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Endtypemarker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Endnumark.
function Endnumark_Callback(hObject, eventdata, handles)
% hObject    handle to Endnumark (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Endnumark contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Endnumark


% --- Executes during object creation, after setting all properties.
function Endnumark_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Endnumark (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chunk.
function chunk_Callback(hObject, eventdata, handles)
% hObject    handle to chunk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

evt = events(handles.Dmeg);

if get(handles.Begstyle,'Value')==1
    %Time Mode

    Begh = str2double(get(handles.Begh,'String'));
    Begm = str2double(get(handles.Begm,'String'));
    Begs = str2double(get(handles.Begs,'String'));

    if isnan(Begh) || isnan(Begm) || isnan(Begs) || ...
            (Begm>60) || (Begs>60) || (Begh<0) || (Begm<0) || (Begs<0)
        beep
        disp('ERROR: The parameters to determine the beginning of the chunk are not correct ')
        return
    end

    if get(handles.Begtimepop,'Value')==1
        %Relative Time

        Begpts=round(min(max(1,(Begh*60^2+Begm*60+Begs)*fsample(handles.Dmeg)),nsamples(handles.Dmeg)));

    else
        %Clock Time

        lgt=(nsamples(handles.Dmeg)/fsample(handles.Dmeg))/60^2; % duration in hours
        ch_sec_rounding = round(handles.Dmeg.info.hour(3)) - handles.Dmeg.info.hour(3);
        if abs(ch_sec_rounding)<1e-4
        % Rounding off if "seconds" saved differs from rounded "seconds" by
        % less than .1ms.
            handles.Dmeg.info.hour(3) = round(handles.Dmeg.info.hour(3));
        end
        datenumstart=datenum([0 0 0 handles.Dmeg.info.hour]);
        datenumBeg=datenum([0 0 0 Begh Begm Begs]);

        datediff=datevec(datenumBeg-datenumstart);
        if datediff(4)>lgt % difference <1ms doesn't count
            datenumBeg=datenum([0 0 1 Begh Begm Begs]);
            datediff=datevec(datenumBeg-datenumstart);
        end

        Begpts=min(nsamples(handles.Dmeg),max(1,round((datediff(4)*60^2+datediff(5)*60+datediff(6))*fsample(handles.Dmeg))));

    end
else
    %Marker Mode

    Begmarknum=handles.Begcurrentlist(get(handles.Begnumark,'Value')); % Marker number x
    Begpts=round(max(min(evt(Begmarknum).time*fsample(handles.Dmeg),nsamples(handles.Dmeg)),1));

end

if get(handles.Endstyle,'Value')==1
    %Time Mode
    Endh=str2double(get(handles.Endh,'String'));
    Endm=str2double(get(handles.Endm,'String'));
    Ends=str2double(get(handles.Ends,'String'));

    if isnan(Endh) || isnan(Endm) || isnan(Ends) || ...
            (Endm>60) || (Ends>60) || (Endh<0) || (Endm<0) || (Ends<0)
        beep
        disp('ERROR: The parameters to determine the end of the chunk are not correct ')
        return
    end

    if get(handles.Endtimepop,'Value')==1
        %Relative Time
        Endpts=round(max(min((Endh*60^2+Endm*60+Ends)*fsample(handles.Dmeg),nsamples(handles.Dmeg)),1));
    else
        %Clock Time

        lgt=(nsamples(handles.Dmeg)/fsample(handles.Dmeg))/60^2;
        datenumstart=datenum([0 0 0 handles.Dmeg.info.hour]);
        datenumEnd=datenum([0 0 0 Endh Endm Ends]);

        datediff=datevec(datenumEnd-datenumstart);
        if datediff(4)>lgt
            datenumEnd=datenum([0 0 1 Endh Endm Ends]);
            datediff=datevec(datenumEnd-datenumstart);
        end

        Endpts=min(nsamples(handles.Dmeg),max(1,round((datediff(4)*60^2+datediff(5)*60+datediff(6))*fsample(handles.Dmeg))));

    end
else
    % Marker Mode

    Endmarknum=handles.Endcurrentlist(get(handles.Endnumark,'Value')); % Marker number x
    Endpts=round(max(min(evt(Endmarknum).time*fsample(handles.Dmeg),nsamples(handles.Dmeg)),1));

end

if Begpts >= Endpts
    beep
    disp('ERROR: The chosen starting point is invalid: it is placed after the chosen ending point');
    return
end

flags = struct('clockt',handles.clocktime); % use defaults for the rest
crc_process_chunk(handles.Dmeg, Begpts, Endpts, flags);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update(hObject,handles)

D = struct(handles.Dmeg);

%% Beg:Time Mode
if get(handles.Begstyle,'Value')==1

    set(handles.Begh,'Style','edit');
    set(handles.Begm,'Style','edit');
    set(handles.Begs,'Style','edit');

    set(handles.Beghtxt,'String','h');
    set(handles.Begmtxt,'String','m');
    set(handles.Begstxt,'String','s');

    if get(handles.Begtimepop,'Value')==1 % Meaning relative time
        set(handles.Begh,'String',num2str(handles.Begrtime(1)));
        set(handles.Begm,'String',num2str(handles.Begrtime(2)));
        set(handles.Begs,'String',num2str(handles.Begrtime(3)));
    else
        set(handles.Begh,'String',num2str(handles.Begctime(1)));
        set(handles.Begm,'String',num2str(handles.Begctime(2)));
        set(handles.Begs,'String',num2str(handles.Begctime(3)));
    end

    set(handles.Begtimepop,'Style','popupmenu');
    if handles.clocktime
        set(handles.Begtimepop,'String',{['Relative time'];['Clock time']})
    else
        set(handles.Begtimepop,'String','Relative time')
    end
    set(handles.Begtypemarker,'Style','text');
    set(handles.Begtypemarker,'String','Not Available');

    %Disabling

    set(handles.Begnumark,'Style','text')
    set(handles.Begnumark,'String','Not Available')

    %% Beg:Marker Mode
else
    % Do not forget to handle when there are no marker
    set(handles.Begtypemarker,'Style','popupmenu');
    set(handles.Begtypemarker,'String',num2str(handles.markerlist'));

    handles.Begidx=get(handles.Begtypemarker,'Value');

    handles.Begcurrentlist=find([D.trials.events(:).value]==handles.markerlist(handles.Begidx));

    set(handles.Begnumark,'Style','popupmenu')
    set(handles.Begnumark,'String',num2str((1:length(handles.Begcurrentlist))'));

    %Disabling

    set(handles.Begtimepop,'Style','text')
    set(handles.Begtimepop,'String','Not Available');

    set(handles.Begh,'Style','text');
    set(handles.Begm,'Style','text');
    set(handles.Begs,'Style','text');

    set(handles.Begh,'String',' ');
    set(handles.Beghtxt,'String',' ');
    set(handles.Begm,'String',' ');
    set(handles.Begmtxt,'String',' ');
    set(handles.Begs,'String',' ');
    set(handles.Begstxt,'String',' ');

end

%% End:Time Mode

if get(handles.Endstyle,'Value')==1

    set(handles.Endh,'Style','edit');
    set(handles.Endm,'Style','edit');
    set(handles.Ends,'Style','edit');

    set(handles.Endhtxt,'String','h');
    set(handles.Endmtxt,'String','m');
    set(handles.Endstxt,'String','s');


    if get(handles.Endtimepop,'Value')==1 % Meaning relative time
        set(handles.Endh,'String',num2str(handles.Endrtime(1)));
        set(handles.Endm,'String',num2str(handles.Endrtime(2)));
        set(handles.Ends,'String',num2str(handles.Endrtime(3)));
    else
        set(handles.Endh,'String',num2str(handles.Endctime(1)));
        set(handles.Endm,'String',num2str(handles.Endctime(2)));
        set(handles.Ends,'String',num2str(handles.Endctime(3)));
    end

    set(handles.Endtimepop,'Style','popupmenu');
    if handles.clocktime
        set(handles.Endtimepop,'String',{['Relative time'];['Clock time']})
    else
        set(handles.Endtimepop,'String','Relative time')
    end
    set(handles.Endtypemarker,'Style','text');
    set(handles.Endtypemarker,'String','Not Available');

    %Disabling

    set(handles.Endnumark,'Style','text')
    set(handles.Endnumark,'String','Not Available')

    %% End:Marker Mode
else
    set(handles.Endtypemarker,'Style','popupmenu');
    set(handles.Endtypemarker,'String',num2str(handles.markerlist'));

    handles.Endidx=get(handles.Endtypemarker,'Value');

    handles.Endcurrentlist=find([D.trials.events(:).value]==handles.markerlist(handles.Endidx));

    set(handles.Endnumark,'Style','popupmenu')
    set(handles.Endnumark,'String',num2str((1:length(handles.Endcurrentlist))'));

    %Disabling

    set(handles.Endtimepop,'Style','text')
    set(handles.Endtimepop,'String','Not Available');

    set(handles.Endh,'Style','text');
    set(handles.Endm,'Style','text');
    set(handles.Ends,'Style','text');

    set(handles.Endh,'String',' ');
    set(handles.Endhtxt,'String',' ');
    set(handles.Endm,'String',' ');
    set(handles.Endmtxt,'String',' ');
    set(handles.Ends,'String',' ');
    set(handles.Endstxt,'String',' ');


end

%% Update handles structure
guidata(hObject, handles);

return




