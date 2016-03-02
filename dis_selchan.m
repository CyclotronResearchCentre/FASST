function varargout = dis_selchan(varargin)
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

% Last Modified by GUIDE v2.5 08-Feb-2012 19:07:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @dis_selchan_OpeningFcn, ...
    'gui_OutputFcn',  @dis_selchan_OutputFcn, ...
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


% --- Executes just before dis_selchan is made visible.
function dis_selchan_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dis_selchan (see VARARGIN)
% Choose default command line output for dis_selchan
set(0,'CurrentFigure',handles.figure1);
handles.output = hObject;

load CRC_electrodes.mat;
handles.names     = names;
handles.pos       = pos';
handles.crc_types = crc_types;

if isempty(varargin) || ~isfield(varargin{1},'file')
    % Filter for vhdr, mat and edf files
    try and(isfield(varargin{1},'multcomp'), varargin{1}.multcomp)
        %if the multiple comparison option is checked
        prefile = spm_select(Inf, 'any', 'Select imported EEG file','' ...
            ,pwd,'\.[mMvVeErR][dDhHaA][fFDdTtwW]');
        set(handles.Selectall,'enable','off','visible','off');
        set(handles.selectMEG,'enable','off','visible','off');
        set(handles.selectEEG,'enable','off','visible','off');
        set(handles.fctother,'enable','off','visible','off');
        set(handles.selectMEGPLAN,'enable','off','visible','off');
        set(handles.desall,'enable','off','visible','off');
        set(handles.checkcICA,'enable','off','visible','off');
        set(handles.Score,'enable','off','visible','off');
%         set(handles.text4,'visible','off');
%         set(handles.text5,'visible','off');
        set(handles.Save,'visible','off');
        handles.multcomp=1;
    catch
        prefile = spm_select(1, 'any', 'Select imported EEG file','' ...
            ,pwd,'\.[mMvVeErR][dDhHaA][fFDdTtwW]');
        handles.multcomp=0;
    end
    for i=1:size(prefile,1)
        D{i} = crc_eeg_load(deblank(prefile(i,:)));
        file = fullfile(D{i}.path,D{i}.fname);
        handles.file{i} = file;
        handles.chan{i} = upper(chanlabels(D{i}));
        if isfield(D{i}, 'info')
            try
                D{i}.info.date;
            catch
                D{i}.info.date = [1 1 1];
            end
            try
                D{i}.info.hour;
            catch
                D{i}.info.hour = [0 0 0];
            end
        else
            D{i}.info = struct('date',[1 1 1],'hour',[0 0 0]);
        end
        handles.Dmeg{i} = D{i};
%         handles.Struct(i) = struct(D{i});
        handles.date{i} = zeros(1,2);
        handles.date{i}(1) = datenum([D{i}.info.date D{i}.info.hour]);
        handles.date{i}(2) = handles.date{i}(1) + ...
                        datenum([ 0 0 0 crc_time_converts(nsamples(D{i})/ ...
                                                            fsample(D{i}))] );
        handles.dates(i,:) = handles.date{i}(:);
    end
    cleargraph(handles)
else
    handles.file = varargin{1}.file;
    prefile = deblank(handles.file);
    index = varargin{1}.index;
    for i=1:size(varargin{1}.Dmeg,2)
%         handles.Dmeg{i} = crc_eeg_load([path(varargin{1}.Dmeg{i}),filesep,fname(varargin{1}.Dmeg{i})]);
        handles.Dmeg{i} = varargin{1}.Dmeg{i};
    end
    if isempty(index)
        index = 1:nchannels(handles.Dmeg{1});
    end
    set(handles.Select,'String',upper(chanlabels(handles.Dmeg{1},varargin{1}.index)));
    diff = setdiff(upper(chanlabels(handles.Dmeg{1})),upper(chanlabels(handles.Dmeg{1},varargin{1}.index)));
    set(handles.Deselect,'String',diff);
    [dumb1,dumb2,index2]=intersect(upper(chanlabels(handles.Dmeg{1},varargin{1}.index)),upper(handles.names));

    idxred  = index2(handles.crc_types(index2)<-1);
    idxblue = index2(handles.crc_types(index2)>-2);

    xred = handles.pos(1,idxred);
    yred = handles.pos(2,idxred);

    xblu = handles.pos(1,idxblue);
    yblu = handles.pos(2,idxblue);

    cleargraph(handles)
    hold on
    plot(xred,yred,'r+'), plot(xblu,yblu,'b+')
    hold off
    if and(length(xblu)==0,length(xred)==0)
        cleargraph(handles)
    end

    xlim([0 1])
    ylim([0 1])
    handles.chan{1}=get(handles.Select,'String');
end

if ~isempty(varargin) && isfield(varargin{1},'delmap')
    handles.delmap=varargin{1}.delmap;
end

if (~isempty(varargin) && isfield(varargin{1},'multcomp') && varargin{1}.multcomp) || ...
        isempty(varargin)
    chanset=handles.chan{1};
    for i=1:size(prefile,1)
        chanset=intersect(chanset,handles.chan{i});
    end
    set(handles.Deselect,'String',chanset);
    handles.chan=chanset;
end

try
    handles.Dmeg{1}.CRC.cICA;
catch
    set(handles.checkcICA,'Visible','off')
end

handles.indmeeg = meegchannels(handles.Dmeg{1});
handles.indeeg = [meegchannels(handles.Dmeg{1},'EEG') meegchannels(handles.Dmeg{1},'LFP')];
handles.indmeg = meegchannels(handles.Dmeg{1},'MEG');
handles.indmegplan = meegchannels(handles.Dmeg{1},'MEGPLANAR');
% handles.indmeeg = indchantype(handles.Dmeg{1},{'MEEG','REF','LFP'});
% handles.indeeg = [indchantype(handles.Dmeg{1},'EEG') indchantype(handles.Dmeg{1},'LFP')];
% handles.indmeg = indchantype(handles.Dmeg{1},'MEG');
% handles.indmegplan = indchantype(handles.Dmeg{1},'MEGPLANAR');
handles.namother = setdiff(chanlabels(handles.Dmeg{1}), chanlabels(handles.Dmeg{1},handles.indmeeg));

if isempty(handles.indeeg)
    set(handles.selectEEG,'enable','off');
end
if isempty(handles.indmeg)
    set(handles.selectMEG,'enable','off');
end
if isempty(handles.indmegplan)
    set(handles.selectMEGPLAN,'enable','off');
end
if isempty(handles.namother)
    set(handles.fctother,'enable','off');
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dis_selchan wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = dis_selchan_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in Select.
function Select_Callback(hObject, eventdata, handles)
% hObject    handle to Select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Select contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Select
contents = get(hObject,'String');
if isempty(contents)
else
    % Remove the "activated" item from the list "Available Channels"
    [dumb1,dumb2,index]=intersect(contents{get(hObject,'Value')},contents);
    temp=[contents(1:index-1) ; contents(index+1:length(contents))];
    set(handles.Select,'String',temp);

    [dumb1,dumb2,index2]=intersect(upper(temp),upper(handles.names));

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

    xlim([0 1])
    ylim([0 1])

    set(handles.Localizer,'XTick',[]);
    set(handles.Localizer,'YTick',[]);

    %Add the "activated" in the list "Selected Channels"
    if length(get(handles.Deselect,'String'))==0
        temp={contents{get(hObject,'Value')}};
    else
        temp=[contents{get(hObject,'Value')} ; get(handles.Deselect,'String')];
    end
    set(handles.Deselect,'String',temp);

    % Prevent crashing if the first/last item of the list is selected.
    set(handles.Select,'Value',max(index-1,1));
    set(handles.Deselect,'Value',1);

    % if multiple comparison, no more than one channel
    contents=get(handles.Select,'String');
    if isfield(handles,'multcomp') && (handles.multcomp && length(contents)>1 || isempty(contents))
        set(handles.PLOT,'enable','off');
        set(handles.PLOT,'ForegroundColor',[1 0 0]);
        beep
        disp('Select only one channel for multiple files comparison')
    elseif isfield(handles,'multcomp') && (handles.multcomp && length(contents)==1)
        set(handles.PLOT,'enable','on');
        set(handles.PLOT,'ForegroundColor',[0 0 0]);
    end
end
% Indicate this selection doesn't come from an selection uploaded
set(handles.Save,'Value',0);
set(handles.load,'Value',0);
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
    % Remove the "activated" item from the list "Available Channels"
    [dumb1,dumb2,index]=intersect(contents{get(hObject,'Value')},contents);
    temp=[contents(1:index-1) ; contents(index+1:length(contents))];
    set(handles.Deselect,'String',temp);

    % Add the "activated" in the list "Selected Channels"
    if length(get(handles.Select,'String'))==0
        temp={contents{get(hObject,'Value')}};
    else
        temp=[contents{get(hObject,'Value')} ; get(handles.Select,'String')];
    end
    set(handles.Select,'String',temp);

    % Prevent crashing if the first/last item of the list is selected.
    set(handles.Deselect,'Value',max(index-1,1));
    set(handles.Select,'Value',1);

    [dumb1,dumb2,index]=intersect(upper(temp),upper(handles.names));

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

    %if multiple comparison, no more than one channel
    contents=get(handles.Select,'String');
    if isfield(handles,'multcomp') && (handles.multcomp && length(contents)>1 || isempty(contents))
        set(handles.PLOT,'enable','off');
        set(handles.PLOT,'ForegroundColor',[1 0 0]);
        beep
        disp('Select only one channel for multiple files comparison')
    elseif isfield(handles,'multcomp') && (handles.multcomp && length(contents)==1)
        set(handles.PLOT,'enable','on');
        set(handles.PLOT,'ForegroundColor',[0 0 0]);
    end
end
% Indicate this selection doesn't come from an selection uploaded
set(handles.Save,'Value',0);
set(handles.load,'Value',0);
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

% --- Executes on button press in Selectall.
function Selectall_Callback(hObject, eventdata, handles)
% hObject    handle to Selectall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.Select,'String',upper(chanlabels(handles.Dmeg{1})));
set(handles.Deselect,'String',cell(0));
if ~isempty(handles.indeeg)
    set(handles.selectEEG,'Value',1);
end
if ~isempty(handles.indmeg)
    set(handles.selectMEG,'Value',1);
end
if ~isempty(handles.indmegplan)
    set(handles.selectMEGPLAN,'Value',1);
end
if ~isempty(handles.namother)
    set(handles.fctother,'Value',1);
end

[dumb1,dumb2,index]=intersect(upper(chanlabels(handles.Dmeg{1})),upper(handles.names));

idxred=index(find(handles.crc_types(index)<-1));
idxblue=index(handles.crc_types(index)>-2);

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

% Indicate this selection doesn't come from an selection uploaded
set(handles.Save,'Value',0);
set(handles.load,'Value',0);
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in selectEEG.
function selectEEG_Callback(hObject, eventdata, handles)
% hObject    handle to selectEEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
H = gcf;
ALL = upper(chanlabels(handles.Dmeg{1}));
valEEG = get(handles.selectEEG,'Value');   %Check if there are already EEG channels selected   
EEG = upper(chanlabels(handles.Dmeg{1}, handles.indeeg));
len = length(EEG);
YET = get(handles.Select,'String')'; %Take the channels already selected
if valEEG==1
    prompt={['Indicate in how many parts you want divide this set if there are ',num2str(len),' EEG channels:'],...
    'Indicate which part you want to see : '};
    name='Select the Signal';
    numlines=1;
    defaultanswer={'1','1'};
    answer=inputdlg(prompt,name,numlines,defaultanswer);
    if ~isempty(answer) %If it's not "cancel" selected
        answernbr = str2double(answer);
        total = length(EEG);
        n = answernbr (1);
        i = answernbr (2);
        if n > total

            h = msgbox ('There is not enough EEG channels to divide by this number , by default all the channels will be selected'); 
        elseif i>n
            h = msgbox ('You can''t do this, by default all the channels will be selected');
        else
            EEG = EEG(ceil(total/n)*(i-1)+1 : ceil(total/n)*i); 
        end
        TOT = [YET EEG];
    else 
        TOT = YET; 
        set(handles.selectEEG,'Value',0)
    end  
else       
    EEG = intersect(upper(YET),upper(EEG));
    TOT = setdiff(YET,EEG);
end
set(0,'CurrentFigure',H);
OTHER = setdiff(ALL,TOT);
set(handles.Select,'String',TOT);
set(handles.Deselect,'String',OTHER);
[dumb1,dumb2,index]=intersect(upper(TOT),upper(handles.names));

idxred=index(find(handles.crc_types(index)<-1));
idxblue=index(find(handles.crc_types(index)>-2));

xred=handles.pos(1,idxred);
yred=handles.pos(2,idxred);

xblu=handles.pos(1,idxblue);
yblu=handles.pos(2,idxblue);

cleargraph(handles);

if  or(length(xblu)~=0,length(xred) ~=0)
    hold on
    plot(xred,yred,'r+'), plot(xblu,yblu,'b+')
    hold off  
end
xlim([0 1])
ylim([0 1])

set(handles.Localizer,'XTick',[]);
set(handles.Localizer,'YTick',[]);
% To avoid warning if there no more available or Selected channels
set(handles.Select,'Value',1);
set(handles.Deselect,'Value',1);
% Indicate this selection doesn't come from an selection uploaded
set(handles.Save,'Value',0);
set(handles.load,'Value',0);
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in selectMEG.
function selectMEG_Callback(hObject, eventdata, handles)
% hObject    handle to selectMEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
H = get(gcf);
megval = get(handles.selectMEG,'Value');
YET = get(handles.Select,'String');
ALL = upper(chanlabels(handles.Dmeg{1})); 
MEG = upper(chanlabels(handles.Dmeg{1}, handles.indmeg));
len = length(MEG);
if megval
    prompt={['Indicate in how many parts you want divide this set if there are ',num2str(len),' MEG channels:'],...
        'Indicate which part you want to see : '};
        name='Select the Signal';
        numlines=1;
        defaultanswer={'1','1'};
    answer=inputdlg(prompt,name,numlines,defaultanswer); 
    if isempty(answer)
        TOT = YET; 
        set(handles.selectMEG,'Value',0);
    else
    answernbr = str2double(answer);
    total = length(MEG);
    n = answernbr (1);
    i = answernbr (2);
    if n > total
        h = msgbox ('There is not enough MEG channels to divide by this number , by default all the channels will be selected'); 
    elseif i>n
        h = msgbox ('You can''t do this, by default all the channels will be selected');
    else
        MEG = MEG(ceil(total/n)*(i-1)+1 : ceil(total/n)*i); 
    end
    TOT = [YET' MEG];
    end
else 
    MEG = intersect(upper(YET),upper(MEG));
    TOT = setdiff(YET,MEG);
end
%set(0,'CurrentFigure',H);
OTHER = setdiff(ALL,TOT);
set(handles.Select,'String',TOT);
set(handles.Deselect,'String',OTHER);

[dumb1,dumb2,index]=intersect(upper(TOT),upper(handles.names));

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
% To avoid warning if there no more available or Selected channels
set(handles.Select,'Value',1);
set(handles.Deselect,'Value',1);
% Indicate this selection doesn't come from an selection uploaded
set(handles.Save,'Value',0);
set(handles.load,'Value',0);
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in selectMEGPLAN.
function selectMEGPLAN_Callback(hObject, eventdata, handles)
% hObject    handle to selectMEGPLAN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
H = get(gcf);
megplanval = get(handles.selectMEGPLAN,'Value');
YET = get(handles.Select,'String');
ALL = upper(chanlabels(handles.Dmeg{1}));
MEGPLAN = upper(chanlabels(handles.Dmeg{1}, handles.indmegplan));
len = length(MEGPLAN);
if megplanval
    prompt={['Indicate in how many parts you want divide this set if there are ',num2str(len),' MEGPLANAR channels: '],...
        'Indicate which part you want to see : '};
        name='Select the Signal';
        numlines=1;
        defaultanswer={'1','1'};
    answer=inputdlg(prompt,name,numlines,defaultanswer);
    if isempty(answer)
        TOT = YET; 
        set(handles.selectMEGPLAN,'Value',0)
    else  
    answernbr = str2double(answer);
    total = length(MEGPLAN);
    n = answernbr (1);
    i = answernbr (2);
    if n > total
    h = msgbox ('There is not enough MEGPLANAR channels to divide by this number , by default all the channels will be selected'); 
    elseif i>n
    h = msgbox ('You can''t do this, by default all the channels will be selected');
    else
    MEGPLAN = MEGPLAN(ceil(total/n)*(i-1)+1 : ceil(total/n)*i); 
    end
    TOT = [YET' MEGPLAN];
    end
else 
    MEGPLAN = intersect(upper(YET),upper(MEGPLAN));
    TOT = setdiff(YET,MEGPLAN);
end
%set(0,'CurrentFigure',H);
OTHER = setdiff(ALL,TOT);
set(handles.Select,'String',TOT);
set(handles.Deselect,'String',OTHER);

[dumb1,dumb2,index]=intersect(upper(TOT),upper(handles.names));

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
% To avoid warning if there no more available or Selected channels
set(handles.Select,'Value',1);
set(handles.Deselect,'Value',1);
% Indicate this selection doesn't come from an selection uploaded
set(handles.Save,'Value',0);
set(handles.load,'Value',0);
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in fctother.
function other_Callback(hObject, eventdata, handles)
% hObject    handle to fctother (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fcto = get(handles.fctother,'Value');
YET = get(handles.Select,'String');
DYET = get(handles.Deselect,'String');
if fcto == 1
    if ~isempty(YET)
        ADD = setdiff(YET,handles.namother');
        SEL = [YET ; ADD];
    else 
        SEL = handles.namother';
    end 
    if ~isempty(DYET)
        DES = setdiff(DYET, handles.namother)';
    else 
        beep
        fprintf('There is no more channel')
    end
else 
    if ~isempty(YET)
        SEL = setdiff(YET, handles.namother);
    else 
        SEL = cell(0);
    end 
    if ~isempty(DYET)
        ADD = setdiff(DYET,handles.namother);
        DES = [ADD'; handles.namother'];
    else 
         DES = [handles.namother];
    end
end 
set(handles.Select,'String',SEL);
set(handles.Deselect,'String',DES);

[dumb1,dumb2,index]=intersect(upper(SEL),upper(handles.names));

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

% To avoid warning if there no more available or Selected channels
set(handles.Select,'Value',1);
set(handles.Deselect,'Value',1);
% Indicate this selection doesn't come from an selection uploaded
set(handles.Save,'Value',0);
set(handles.load,'Value',0);
% Update handles structure
guidata(hObject, handles);


% Hint: get(hObject,'Value') returns toggle state of fctother
% --- Executes on button press in desall.
function desall_Callback(hObject, eventdata, handles)
% hObject    handle to desall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.selectEEG,'Value',0);
set(handles.selectMEG,'Value',0);
set(handles.selectMEGPLAN,'Value',0);
set(handles.fctother,'Value',0);
set(handles.Deselect,'String',upper(chanlabels(handles.Dmeg{1})));
set(handles.Select,'String',cell(0));
set(handles.fctother,'Value',0);
cleargraph(handles)
xlim([0 1])
ylim([0 1])
set(handles.Localizer,'XTick',[]);
set(handles.Localizer,'YTick',[]);
% Indicate this selection doesn't come from an selection uploaded
set(handles.Save,'Value',0);
set(handles.load,'Value',0);
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in PLOT
function PLOT_Callback(hObject, eventdata, handles)
% hObject    handle to PLOT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = get(handles.Select,'String');
   
if isfield(handles,'multcomp') && (handles.multcomp && length(contents)>1 || isempty(contents))%no more than one channel if multiple file comparison
    beep
    disp('Select only one channel for multiple files comparison')
    return
elseif ~isfield(handles,'multcomp') || ~handles.multcomp
    [dumb1,index] = intersect( upper(chanlabels(handles.Dmeg{1})), ...
        upper(get(handles.Select,'String')));
    try
        flags.index = sortch(handles.Dmeg{1},index);
    catch
        flags.index = fliplr(sort(index));
    end
else
    flags.dates=handles.dates;
    flags.chanset=handles.chan;
    [dumb1,index]=intersect(upper(handles.chan),upper(get(handles.Select,'String')));
    flags.index=index;
    flags.multcomp=1;
end
flags.Dmeg=handles.Dmeg;
flags.file=handles.file;
if isfield(handles,'delmap')
    flags.delmap=handles.delmap;
end
crc_dis_main(flags);
delete(handles.figure1)

% --- Executes on button press in Score.
function Score_Callback(hObject, eventdata, handles)
% hObject    handle to Score (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[dumb1,index,dumb2]=intersect(upper(chanlabels(handles.Dmeg{1})),get(handles.Select,'String'));
try
    flags.index=sortch(handles.Dmeg{1},index);
catch
    flags.index=fliplr(sort(index));
end

flags.Dmeg=handles.Dmeg;
flags.file=handles.file;
flags.scoresleep=1;
if isfield(handles,'delmap')
    flags.delmap=handles.delmap;
end
crc_dis_main(flags);
delete(handles.figure1)

% --- Executes on button press in checkcICA.
function checkcICA_Callback(hObject, eventdata, handles)
% hObject    handle to checkcICA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[dumb1,index,dumb2] = intersect(upper(chanlabels(handles.Dmeg{1})), ...
                                get(handles.Select,'String'));
try
    flags.index=sortch(handles.Dmeg{1},index);
catch
    flags.index=fliplr(sort(index))  ;
end
flags.Dmeg = handles.Dmeg;
flags.file = handles.file;
dis_cICAcheck(flags);
delete(handles.figure1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SORT CHANNEL FX
function idx = sortch(d,index)

% Find the EOG channels
EOGchan = eogchannels(d);
EOGchan = intersect(EOGchan,index);

% Find the ECG channels
ECGchan = ecgchannels(d);
ECGchan = intersect(ECGchan,index);

% Find the EMG channels
EMGchan = emgchannels(d);
EMGchan = intersect(EMGchan,index);

% Find the M/EEG channels
EEGchan = meegchannels(d);
EEGchan = intersect(EEGchan,index);

allbad = [EOGchan ECGchan EMGchan EEGchan];
other  = setdiff(index,allbad);

otherknown    = [];
othernotknown = [];
chanstr = chantype(d);
for ff = other
    if ~strcmpi(chanstr,'Other')
        othernotknown = [othernotknown ff];
    else
        otherknown = [otherknown ff];
    end
end
other = [otherknown othernotknown];

allbad = [EOGchan ECGchan EMGchan other];
eeg = setdiff(index,allbad);

AFrontal  = intersect(find(strncmp(chanlabels(d),'AF',2) ==1),eeg);
Frontal   = intersect(find(strncmp(chanlabels(d),'F',1) ==1),eeg);
Coronal   = intersect(find(strncmp(chanlabels(d),'C',1) ==1),eeg);
Temporal  = intersect(find(strncmp(chanlabels(d),'T',1) ==1),eeg);
Parietal  = intersect(find(strncmp(chanlabels(d),'P',1) ==1),eeg);
Occipital = intersect(find(strncmp(chanlabels(d),'O',1) ==1),eeg);

neweeg = [Occipital Parietal Temporal Coronal Frontal AFrontal];
eeg2 = setdiff(eeg,neweeg);

eeg = [eeg2 neweeg];

idx = [otherknown othernotknown ECGchan EMGchan EOGchan eeg];

%%%%%%%%%%%%%%%
function cleargraph(handles)

A=get(handles.figure1,'Children');
idx=find(strcmp(get(A,'Type'),'axes')==1);
try
    delete(get(A(idx),'Children'))
end

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


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



% --- Executes on button press in Saveselection.
function Saveselection_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Saveselection
S = get(handles.Save,'Value');
[filename, pathname] = uiputfile('.txt','Save the Selection');
file = fopen(fullfile(pathname, filename), 'w');              % 'w' write
data = get(handles.Select,'String');
%write the data in file .txt
%data = {'bonjour' 'dorothe'}
for i = 1 : length(data)
    count = fwrite(file, [data{i} ' '],'char');
    if count ~= length(data{i})+1
        beep
        fprintf('There is an error, we can''t save this selection... Try again')
        i = length(data);
    end 
end
fprintf(['This selection is saved under the name:', char(filename)]);
fclose(file);
% --- Executes on button press in loadselection.
function loadselection_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of loadselection

[filename, pathname] = uigetfile('.txt','Open the Selection');
namefile = fullfile(pathname,filename);
fid = fopen (namefile,'r');
file = fread(fid,'uint8=>char')';
fclose(fid);
i = 1;
Selection = cell(0);
while i<length(file)
    init = i;
    while ~strcmpi(file(i),' ') && i<length(file)      
        i=i+1;
    end 
    Selection = [Selection file(init : i-1)];
    i=i+1;   
end
YET = get(handles.Select,'String');
if ~isempty(YET) 
    Sel = setdiff(YET,Selection);
    if ~isempty(Sel)
        Selection = [Sel'; Selection'];
    end
end
set(handles.Select,'String',Selection);
% Select all the channels there are in the file
contents = upper(chanlabels(handles.Dmeg{1}));
[dumb1,index,dumb2] = intersect(contents,Selection);

% Remove channels concerned by the Selection
DES = get(handles.Deselect,'String');
temp = setdiff(DES,Selection);
set(handles.Deselect,'String',temp);
% To avoid warning
set(handles.Select,'Value',1);
set(handles.Deselect,'Value',1);
[dumb1,dumb2,index2]=intersect(upper(Selection),upper(handles.names));
% Indicate the placement of the electrodes selected
idxred  = index2(find(handles.crc_types(index2)<-1));
idxblue = index2(find(handles.crc_types(index2)>-2));

xred = handles.pos(1,idxred);
yred = handles.pos(2,idxred);

xblu = handles.pos(1,idxblue);
yblu = handles.pos(2,idxblue);

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

set(handles.load,'Value',0)






