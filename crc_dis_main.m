function varargout = crc_dis_main(varargin)
% CRC_DIS_MAIN M-file for crc_dis_main.fig
% CRC_DIS_MAIN, by itself, creates a new CRC_DIS_MAIN or raises the
% existing singleton*.
%
%      H = CRC_DIS_MAIN returns the handle to a new CRC_DIS_MAIN or the handle to
%      the existing singleton*.
%
%      CRC_DIS_MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CRC_DIS_MAIN.M with the given input arguments.
%
%      CRC_DIS_MAIN('Property','Value',...) creates a new CRC_DIS_MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dis_main_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to crc_dis_main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%__________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by Y. Leclercq & C. Phillips, 2008.
% Modified by J. Schrouff, 2010-...
% And further updated by D. Coppieters, 2012-...
% Cyclotron Research Centre, University of Liege, Belgium
% $Id:$

% Edit the frqabv text to modify the response to help crc_dis_main

% Last Modified by GUIDE v2.5 14-Jan-2011 17:09:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @crc_dis_main_OpeningFcn, ...
    'gui_OutputFcn',  @crc_dis_main_OutputFcn, ...
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

% --- Executes just before crc_dis_main is made visible.
function crc_dis_main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to crc_dis_main (see VARARGIN)

crcdef = crc_get_defaults('one');
set(0,'CurrentFigure',handles.figure1);
warning('off')

% Choose default command line output for crc_dis_main
handles.output  =   hObject;
handles.file    =   varargin{1}.file;
handles.figz    =   0;
handles.Dmeg    =   varargin{1}.Dmeg;


if isfield(varargin{1},'showOracle')
    handles.showOracle     =   1;
else
    handles.showOracle     =   0;
end

if isfield(varargin{1},'scoresleep')
    handles.scoring     =   1;
else
    handles.scoring     =   0;
end

handles.gui_active = 1;

for i=1:numel(handles.Dmeg)
    handles.Dmeg{i} =   meeg(handles.Dmeg{i});
end
D = handles.Dmeg{1};
% distinct index to make easier to select data by type
handles.index   =   varargin{1}.index;
handles.indexMEEG   =   fliplr(intersect(meegchannels(handles.Dmeg{1}),handles.index))';
handles.indnomeeg   =   setdiff(handles.index, handles.indexMEEG)';
if size(handles.indnomeeg,1)>1
    handles.indnomeeg = handles.indnomeeg';
end
if size(handles.indexMEEG,1)>1
    handles.indexMEEG = handles.indexMEEG';
end
handles.inddis      =   [handles.indnomeeg handles.indexMEEG];
% handles.inddis      =   [handles.index];
handles.indeeg = [meegchannels(D,'EEG') meegchannels(D,'LFP')];
% Default values
set(handles.normalize,'Value',0);

% Take the window size used to score
if handles.scoring % scoring sleep data
    if isfield (handles.Dmeg{1},'CRC')
        if isfield (handles.Dmeg{1}.CRC,'score')
            handles.winsize   =   handles.Dmeg{1}.CRC.score{3,1};
        else
            handles.winsize   =   crc_get_defaults('score.winsize');
        end
    else
        handles.winsize     =   crc_get_defaults('score.winsize');
    end
else % no scoring
    handles.winsize     =   crcdef.winsize;
end

handles.scale       =   crcdef.scale;
handles.eegscale    =   [];
handles.eogscale    =   [];
handles.emgscale    =   [];
handles.ecgscale    =   [];
handles.otherscale  =   [];
handles.lfpscale    =   [];
handles.megmagscale =   [];
handles.megplanarscale  =   [];
handles.chan        =   [];

%si ca revient de la détection automatique
handles.artefacteeg = [];
if isfield(handles.Dmeg{1},'CRC')
    if isfield (handles.Dmeg{1}.CRC,'artefacteeg')
        handles.artefacteeg = handles.Dmeg{1}.CRC.artefacteeg;
    end
end
handles.artefact = [];
if isfield(D,'CRC')
    if isfield(D.CRC,'DC')
        if isfield(D.CRC.DC,'shortartf')
            if isfield(D.CRC.DC.shortartf,'total')
                handles.artefact = D.CRC.DC.shortartf.total;
            end
        end
    end
end

% badchannels
handles.badchannels = [];
if isfield(D,'CRC')
    if isfield(D.CRC,'DC')
        if isfield(D.CRC.DC,'badchannels')
            if isfield(D.CRC.DC.badchannels,'chan_defaillant')
                badchannels_def = D.CRC.DC.badchannels.chan_defaillant;
                badchannels_incoh = D.CRC.DC.badchannels.chan_incoherent;
                handles.badchannels = union(badchannels_def',badchannels_incoh','rows')';
            end
        end
    end
end
handles.chanlab = meegchannels(D);
if isempty(handles.chanlab)
    handles.chanlab = find(or(or(strncmp(chanlabels(D),'F',1),strncmp(chanlabels(D),'C',1)),or(strncmp(chanlabels(D),'O',1),strncmp(chanlabels(D),'P',1))));
end
%Spike detection
%---------------
handles.move = 0;
if  isfield(varargin{1},'type') % We come back from the "Event_menu"
    handles.type = varargin{1}.type;
else       % It's the begining
    D      = handles.Dmeg{1};
    handles.type = {};
    if isfield(D,'CRC')
        if isfield(D.CRC,'Event')&& ~isempty(events(D)) %changes made
            type_sauv = D.CRC.Event;
            tp = events(D);
            tpt = {tp(:).type};
            [goodev int_tpt int_sauv] = intersect(tpt,type_sauv(:,1));
            handles.type = type_sauv(int_sauv,:);
        else
            evs    = events(D);
            if ~isempty(evs)
                valev  = {evs(:).value};
                type = cell(0);
                evm = 1;
                if isfield(D.CRC,'score')
                    nsc    = size(D.CRC.score,2);
                    type = cell(0);
                    for insc = 1 : nsc
                        for vv = 1 : length(valev)
                            if any(strcmpi( {(D.CRC.score{2,insc})},valev{vv}))||strcmpi('Newuser',valev{vv})
                                type(evm) = {evs(vv).type};
                                evm = evm + 1;
                            end
                        end
                    end
                else
                    for vv = 1 : length(valev)
                        if any(strcmpi('Newuser',valev{vv}))
                            type(evm) = {evs(vv).type};
                            evm = evm + 1;
                        end
                    end
                end
                tevm = cellstr(unique(char(type(:)),'rows'));
                ntevm = size(tevm,1);
                if ~isempty(type(:))
                    handles.type(1:ntevm,1)  = tevm;
                    handles.type(1:ntevm,2)  = cellstr('red');
                end
            end
        end
    end
end

[path, ~,~] = fileparts([mfilename('fullpath') '.m']);
settings_path = fullfile(path,'/settings.mat');

if exist(settings_path, 'file') == 2,
    load(settings_path);
    set(handles.validity,'string',['License Validity: ' settings.validity ' UTC.'])
else
    set(handles.validity,'string','License Validity: NA')
end

delete(get(handles.manevent,'Children'));
for itp = 1:size(handles.type,1)
    uimenu(handles.manevent,'Label', char(handles.type(itp,1)),'Callback',@Define_event)
end

uimenu(handles.manevent,'Label', ...
    'New Type','Callback', @newtype_Callback,...
    'Separator','on') ;
%--------------------------------------------
%coming back from "detection" or "Event_menu"

%-------------------------------------------

if isfield(varargin{1},'base')
    handles.base     =   varargin{1}.base;
    
else
    handles.base     =   {};
end

%---------------------------------------------

if isfield(varargin{1},'delmap')
    set(handles.delaymap,'visible','on');
    set(handles.delaymap,'Value',1);
    if isfield(varargin{1}.delmap,'elpos')
        handles.delmap.elpos    =   varargin{1}.delmap.elpos;
    end
    if isfield(varargin{1}.delmap,'zebris_name')
        handles.delmap.zebris_name	=   varargin{1}.delmap.zebris_name;
    end
    if isfield(varargin{1}.delmap,'cas')
        handles.delmap.cas  =   varargin{1}.delmap.cas;
    end
else
    set(handles.delaymap,'visible','off');
    set(handles.text25,'visible','off');
    set(handles.delaymap,'Value',0);
    handles.delmap  =   [];
end

if isfield(varargin{1},'multcomp') % if comparing multiple files
    handles.multcomp    =   1;
    handles.displevt    =   0;
    handles.export      =   0;
    handles.hor_grid    =   0;
    handles.vert_grid   =   0;
    set(handles.figure1, 'windowbuttonmotionfcn', '')
    %Handles date and hour and slider
    handles.date        =   varargin{1}.dates;
    handles.chanset     =   varargin{1}.chanset;
    maxdate     = max(max(handles.date));
    handles.maxdate     =   maxdate;
    mindate     =   min(min(handles.date));
    handles.mindate     =   mindate;
    mindatevec	=   datevec(min(min(handles.date)));
    diffe        =   maxdate-mindate;
    diffvec     =   datevec(diffe);
    if diffvec(3)+(diffvec(2)-1)*30+(diffvec(1)-1)*365 > 2
        handles.offset  =   0;
    else
        handles.offset  =   mindatevec(4) * 60^2 + mindatevec(5)*60 + mindatevec(6);
    end
    handles.maxx    =   diffvec(4)*60^2+diffvec(5)*60+diffvec(6);
    set(handles.slider1,'Value',1/handles.winsize);
    set(handles.slider1,'Max',handles.maxx-handles.winsize);
    set(handles.slider1,'Min',1/handles.winsize);
    set(handles.slider1,'SliderStep',[handles.winsize/(handles.maxx) 0.1]);
    set(handles.totaltime,'String',['/ ' num2str(handles.maxx)]);
    set(handles.currenttime,'String',num2str(round(handles.winsize/2)));
    set(handles.totalpage,'String',['/ ' num2str(ceil(handles.maxx/handles.winsize))]);
    set(handles.currentpage,'String',num2str(1));
    
    %disables channel slider
    set(handles.Chanslider,'enable','off')
    
    %handles popupmenustring using common channels
    popupmenustring = handles.chanset;
    if sum(strcmp(handles.chanset,'REF2'))
        popupmenustring     =   [popupmenustring 'MEAN OF REF'];
    end
    
    %handles visibility of different options
    set(handles.Score,'enable','off','visible','off')
    set(handles.pushbutton1,'enable','off','visible','off')
    set(handles.pushbutton2,'enable','off','visible','off')
    set(handles.popupmenu10,'enable','off','visible','off')
    set(handles.popupmenu11,'enable','off','visible','off')
    set(handles.radiobutton1,'enable','off','visible','off')
    set(handles.pushbutton5,'enable','off','visible','off')
    set(handles.pushbutton6,'enable','off','visible','off')
    set(handles.pushbutton7,'enable','off','visible','off')
    set(handles.NbreChan,'enable','off','visible','off')
    set(handles.uipanel6,'visible','off')
    set(handles.text22,'visible','off')
    set(handles.text23,'visible','off')
    set(handles.NbrechanTxt,'visible','off')
    set(handles.NbreChan,'Value',1)
    set(handles.Cmp_Pwr_Sp_All,'enable','on','visible','on')
    set(handles.multfil,'enable','on')
    set(handles.multnames,'enable','on')
    set(handles.multchan,'enable','on')
    set(handles.multother,'enable','on')
    set(handles.multclose,'enable','on')
    set(handles.axes4,'visible','off')
    set(handles.addundefart,'visible','off')
    set(handles.addspecart,'visible','off')
    set(handles.addaro,'visible','off')
    set(handles.addeoi,'visible','off')
    set(handles.FPL,'visible','off')
    set(handles.OPL,'visible','off')
    set(handles.Delart,'visible','off')
    set(handles.Delaro,'visible','off')
    set(handles.Del_eoi,'visible','off')
    set(handles.Delartonlyone,'visible','off')
    %set(handles.Detection,'visible','off')
    set(handles.addonlyone,'visible','off')
    %     set(handles.others,'enable','off','visible','off')
    %     set(handles.com,'enable','off','visible','off')
    set(handles.manevent,'enable','off','visible','off')
    set(handles.event_menu,'enable','off')
    set(handles.counterspect,'enable','off','visible','off')
    for ii  =   1   :   size(handles.chanset,2)
        handles.multchanlab{ii}     =   uimenu(handles.multchan,'Label',...
            char(handles.chanset{ii}),'Callback',...
            {@Chantodisp});
        if ii   ==  handles.index
            set(handles.multchanlab{ii},'Checked','on')
            handles.Chantodis   =   ii;
        end
    end
    set(handles.figure1,'name','Multiple files comparison')
    
else %if only one file to display
    %set(handles.Detection,'visible','off') %Not ready yet
    handles.multcomp=0;
    set(handles.axes5,'visible','on');
    set(handles.axes4,'visible','on');
    
    %handles date and hour and slider
    if isfield(handles.Dmeg{1},'info')
        if isfield(handles.Dmeg{1}.info,'hour')
            handles.offset = handles.Dmeg{1}.info.hour(1) * 60^2 + ...
                handles.Dmeg{1}.info.hour(2)*60 + ...
                handles.Dmeg{1}.info.hour(3);
        end
    end
    set(handles.slider1,...
        'Max', nsamples(handles.Dmeg{1})/fsample(handles.Dmeg{1})-handles.winsize/2,...
        'Value',0,...
        'Min',0)
    try
        set(handles.slider1,'SliderStep',[(handles.winsize)/(nsamples(handles.Dmeg{1})/fsample(handles.Dmeg{1})-handles.winsize/2) 0.1])
    end
    set(handles.totaltime,'String',['/ ' ...
        num2str(round(nsamples(handles.Dmeg{1})/fsample(handles.Dmeg{1})))]);
    set(handles.currenttime,'String',num2str(round(handles.winsize/2)));
    %new index for page
    set(handles.totalpage,'String',['/ ' ...
        num2str(ceil((nsamples(handles.Dmeg{1})/fsample(handles.Dmeg{1}))/handles.winsize))]);
    set(handles.currentpage,'String',num2str(1));
    % Number of channels displayed
    Nchdisp     =   str2double(get(handles.NbreChan,'String'));
    set(handles.NbreChan,'String', num2str(min(Nchdisp, length(handles.indexMEEG) + length(handles.indnomeeg))))
    Nchdisp     =   str2double(get(handles.NbreChan,'String'));
    totind      =   numel(handles.index); %length(handles.indexMEEG) + length(handles.indnomeeg);
    set(handles.Chanslider,...
        'Min',1,...
        'Max',totind - Nchdisp+1,...
        'Value',1,...
        'SliderStep', [Nchdisp/totind Nchdisp/totind])
    
    popupmenustring     =   chanlabels(handles.Dmeg{1}) ;
    if sum(strcmp(chanlabels(handles.Dmeg{1}),'REF2'))
        popupmenustring =   [popupmenustring 'MEAN OF REF'];
    end
    if and(sum(strcmp(chanlabels(handles.Dmeg{1}),'M1')), ...
            sum(strcmp(chanlabels(handles.Dmeg{1}),'M2')))
        popupmenustring =   [popupmenustring 'M1-M2'];
    end
    set(handles.figure1,'name',handles.file{1})
end

% Define filter according to the smallest sampling frequency
handles.filter.other    =   crcdef.filtEEG;
filt=zeros(size(handles.Dmeg,2),1);
for i=1:size(handles.Dmeg,2)
    filt(i)     =   fsample(handles.Dmeg{i})/2;
end
[filt,posf] = max(filt);
% Adjust filter order to sampling rate, low at 5kHz (EEG-fMRI data)
if filt>2000
    forder = 1;
else
    forder = 3;
end
handles.minsamp     =   posf;
handles.filter.EMG  =   [crcdef.filtEMG(1) ...
    min(crcdef.filtEMG(2),filt/2)];
handles.filter.EOG = crcdef.filtEOG;
[B,A] = butter(forder,[handles.filter.EOG(1)/(fsample(handles.Dmeg{posf})/2),...
    handles.filter.EOG(2)/(fsample(handles.Dmeg{posf})/2)],'pass');
handles.filter.coeffEOG=[B;A];
[B,A] = butter(forder,[handles.filter.EMG(1)/(fsample(handles.Dmeg{posf})/2),...
    handles.filter.EMG(2)/(fsample(handles.Dmeg{posf})/2)],'pass');
handles.filter.coeffEMG=[B;A];
[B,A] = butter(forder,[handles.filter.other(1)/(fsample(handles.Dmeg{posf})/2),...
    handles.filter.other(2)/(fsample(handles.Dmeg{posf})/2)],'pass');
handles.filter.coeffother=[B;A];
set(handles.frqabv,'String',num2str(handles.filter.other(2)));
set(handles.upemg,'String',num2str(handles.filter.EMG(2)));

% popupmenu setting
popupmenustring	=   [popupmenustring 'REF1'];
set(handles.EEGpopmenu,...
    'String',popupmenustring,...
    'Value',length(popupmenustring))
set(handles.otherpopmenu,...
    'String',popupmenustring,...
    'Value',length(popupmenustring))
popupmenustring = [popupmenustring 'BIPOLAR'];
set(handles.EOGpopmenu,...
    'String',popupmenustring,...
    'Value',length(popupmenustring))
set(handles.EMGpopmenu,...
    'String',popupmenustring,...
    'Value',length(popupmenustring))

% Menu to display FFT for each window according to the channel selected
popmenustring   =   chanlabels(handles.Dmeg{1}, handles.inddis);
set(handles.fftchan, 'String', popmenustring,'Value', 1,'Visible','off')

load('CRC_electrodes.mat');
handles.names       =   names;
handles.crc_types   =   crc_types;
handles.pos         =   pos';

%events display, only if one file
if ~handles.multcomp
    handles.displevt    =   0;
    pmstring    =   [{'All'}];
    evt         =   events(handles.Dmeg{1});
    evt = evt(:);
    if iscell(evt)
        evt     =   cell2mat(evt);
        disp(['Warning: data not continuous (trials of 1s), only first epoch showed'])
    end
    if ~isempty(evt)
        for i = 1:max(size(evt,2),size(evt,1))
            if ~isempty(evt(i)) && ~any(strcmpi(evt(i).type, pmstring)) &&...
                    ~isempty(evt(i).time)
                pmstring    =   [pmstring, {evt(i).type}];
            end
        end
    end
    set(handles.popupmenu10,...
        'String',pmstring,...
        'Value',length(pmstring))
    pmstring=[{'All'}, {'Number of event'}];
    if ~isempty(evt)
        for i = 1:size(evt,2)
            if ~isempty(evt(i)) && isnumeric(evt(i).value)
                evt(i).value    =   num2str(evt(i).value);
            end
            if ~isempty(evt(i)) && ~any(strcmpi(evt(i).value, pmstring)) &&...
                    ~isempty(evt(i).time)
                pmstring    =   [pmstring, {evt(i).value}];
            end
        end
    end
    set(handles.popupmenu11,...
        'String',pmstring,...
        'Value',length(pmstring))
    handles.evt     =   evt;
    handles.chostype    =   1:max(size(evt,2),size(evt,1));
    handles.chosevt     =   1:max(size(evt,2),size(evt,1));
    if isfield(handles.Dmeg{1},'CRC') && ...
            isfield(handles.Dmeg{1}.CRC,'lastdisp')
        slidval=handles.Dmeg{1}.CRC.lastdisp;
        set(handles.slider1,'Value',slidval)
    end
    set(handles.figure1,'MenuBar','figure')
    
    %Scoring information if only one file
    handles.vert_grid   =   0;
    handles.hor_grid    =   0;
    handles.export      =   0;
    set(handles.figure1, 'windowbuttonmotionfcn', @update_powerspect)
    if handles.scoring
        Score_Callback(hObject,eventdata,handles)
        handles         =   guidata(hObject);
    else
        handles.score   =   cell(8,1);
        handles.score{4,1}  =   [1/fsample(handles.Dmeg{1}) ...
            nsamples(handles.Dmeg{1})/fsample(handles.Dmeg{1})];
        handles.score{2,1}  =   'Newuser';
        handles.score{1,1}  =   [0/0, 0/0];
        set(handles.figure1,'CurrentAxes',handles.axes4)
        %Select the scorer according to the last used
        if isfield(varargin{1},'user')
            handles.currentscore = varargin{1}.user;
        else
            handles.currentscore    =   1;
        end
        
        crc_hypnoplot(handles.axes4, ...
            handles,handles.winsize)
        
        set(handles.figure1,'CurrentAxes',handles.axes1)
        set(handles.addundefart,'visible','off')
        set(handles.addonlyone,'visible','off')
        set(handles.addspecart,'visible','off')
        set(handles.addaro,'visible','off')
        set(handles.addeoi,'visible','off')
        set(handles.FPL,'visible','off')
        set(handles.OPL,'visible','off')
        set(handles.Delart,'visible','off')
        set(handles.Delartonlyone,'visible','off')
        set(handles.Delaro,'visible','off')
        set(handles.Del_eoi,'visible','off')
        set(handles.manevent,'visible','off') %To make visible to put events manually
        %set(handles.Detection,'visible','on')
    end
end
D   =   handles.Dmeg{1};

if ~isfield(D,'info')
    D.info  =   [];     %in case field other empty set info to allow
end                        %further dot name structure assignments

if ~isfield (D,'CRC')
    D.CRC   =   [];
end

% initalize the field 'goodevents' in D.CRC and, if spindle or slow wave
% detection was performed, retrieve which events were marked as bad.
aevt = events(D);
D.CRC.goodevents    =   ones(1,numel(aevt));
if iscell(aevt)
    aevt   =   cell2mat(aevt);
end

if ~isempty(aevt)
    evtype  =   [{aevt(:).type}];
    if isfield(D.CRC,'SW')
        indsw   =   [];
        for j   =   1   :   size(evtype,2)
            if ~isempty(strfind(evtype{j},'SW'))|| ~isempty(strfind(evtype{j},'delta'))
                indsw	=   [indsw,j];
            end
        end
        for i   =  1    :   size(D.CRC.SW.SW,2)
            a   =   [aevt(indsw).time];
            [d1,b]  =   sort(abs(round(a)-round([D.CRC.SW.SW(i).negmax]/1000)));
            if ~isfield(D.CRC.SW.SW(i),'good')
                D.CRC.SW.SW(i).good=1;
            end
            D.CRC.goodevents(indsw(b(1)))   =   D.CRC.SW.SW(i).good;
        end
    end
    if isfield(D.CRC,'spindles')
        indsp   =   [];
        for j   =   1   :   size(evtype,2)
            if ~isempty(strfind(evtype{j},'SP'))
                indsp   =  [indsp,j];
            end
        end
        for i   =   1   :   size(D.CRC.spindles.bounds,1)
            a   =   [aevt(indsp).time];
            [d1,b]  =   sort(abs(round(a*10)-round(D.CRC.spindles.bounds(i,1)/fsample(handles.Dmeg{1})*10)));
            if ~isfield(D.CRC.spindles,'good')
                D.CRC.spindles.good(i)=1;
            end
            D.CRC.goodevents(indsp(b(1)))   =   D.CRC.spindles.good(i);
        end
    end
end

handles.Dmeg{1}     =   D;
save(D);

% to display the main plot
mainplot(handles)
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes crc_dis_main wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = crc_dis_main_OutputFcn(hObject, eventdata, handles)
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

handles.displevt=0;
%removed and then readd%%%%%%%%
if isfield(handles,'scoring') & handles.scoring
    set(handles.figure1,'CurrentAxes',handles.axes4)
    crc_hypnoplot(handles.axes4, ...
        handles,handles.score{3,handles.currentscore})
    set(handles.figure1,'CurrentAxes',handles.axes1)
end
%%%%%%%%%%%%%%%%%
mainplot(handles)
% Update handles structure
guidata(hObject, handles);

function confslider_Callback(hObject, eventdata, handles)
value = get(hObject,'Value');
set(handles.conftext,'String',sprintf('Confidence Threshold %.1f',value));
%removed and then readd%%%%%%%%
if isfield(handles,'scoring') & handles.scoring
    set(handles.figure1,'CurrentAxes',handles.axes4)
    crc_hypnoplot(handles.axes4, ...
        handles,handles.score{3,handles.currentscore})
    set(handles.figure1,'CurrentAxes',handles.axes1)
end
mainplot(handles)
guidata(hObject, handles);


function confslider_CreateFcn(hObject, eventdata, handles)
set(hObject,'Value',1.5);
%set(handles.conftext,'String','Confidence Threshold: 1.5');

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

%Editing the window size
function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

if handles.multcomp
    handles.winsize = min(str2double(get(hObject,'String')),handles.maxx);
    maxiwin         = handles.maxx-handles.winsize;
    miniwin         = 1/handles.winsize;
else
    handles.winsize = min(str2double(get(hObject,'String')), ...
        round(nsamples(handles.Dmeg{1})/(2*fsample(handles.Dmeg{1}))));
    maxiwin         = nsamples(handles.Dmeg{1})/ ...
        fsample(handles.Dmeg{1})-handles.winsize;
    miniwin         = 1/fsample(handles.Dmeg{1});
end
set(handles.edit1,'String',num2str(handles.winsize));
set(handles.slider1,'Max',maxiwin)
set(handles.slider1,'Min',miniwin)
set(handles.slider1,'SliderStep',[handles.winsize/maxiwin 0.1])
set(handles.totalpage,'String',['/ ' ...
    num2str(ceil((nsamples(handles.Dmeg{1})/fsample(handles.Dmeg{1}))/handles.winsize))]);
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

%editing the scale
function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double

if not(str2double(get(hObject,'String'))>0)&& not(str2double(get(hObject,'String'))<0)
    set(handles.edit2,'String',num2str(handles.scale));
else
    handles.scale=str2double(get(hObject,'String'));
end

% Update handles structure
guidata(hObject, handles);

mainplot(handles)

if handles.multcomp
    i   =   length(handles.Dmeg);
else
    Chanslidval =   get(handles.Chanslider,'Value');
    slidpos     =   Chanslidval-rem(Chanslidval,1);
    NbreChandisp    =   str2double(get(handles.NbreChan,'String'));
    i   =   min(length(handles.indexMEEG) + length(handles.indnomeeg),NbreChandisp);
end
ylim([0 handles.scale*(i+1)]);
set(handles.axes1,'YTick',[handles.scale/2:handles.scale/2:i*handles.scale+handles.scale/2]);
ylabels=[num2str(round(handles.scale/2))];
for j=1:i
    if handles.multcomp
        ylabels=[ylabels {num2str(j)}];
    else
        ylabels=[ylabels chanlabels(handles.Dmeg{1},handles.index(j))];
    end
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

%Push the 'back to channel selection' button
% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

flags.index	=   handles.index;
flags.Dmeg  =   handles.Dmeg;
flags.file  =   handles.file;
dis_selchan(flags);

delete(handles.figure1);
try
    close(handles.figz)
end

%Editing the 'current time' field
function currenttime_Callback(hObject, eventdata, handles)
% hObject    handle to currenttime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of currenttime as text
%        str2double(get(hObject,'String')) returns contents of currenttime as a double

i=handles.minsamp;
slidval = str2double(get(hObject,'String')); % str2double is faster than str2num
slidval = max(slidval - handles.winsize/2,1/fsample(handles.Dmeg{i}));

if handles.multcomp
    maxiwin=handles.maxx-handles.winsize;
else
    maxiwin=nsamples(handles.Dmeg{1})/fsample(handles.Dmeg{1})-handles.winsize;
end
slidval = min(slidval,maxiwin);
set(handles.slider1,'Value',slidval)

%Remove and then readd...to be checked
if handles.scoring
    set(handles.figure1,'CurrentAxes',handles.axes4)
    crc_hypnoplot(handles.axes4, ...
        handles,handles.score{3,handles.currentscore})
    set(handles.figure1,'CurrentAxes',handles.axes1)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%Editing the 'current time' field
function currentpage_Callback(hObject, eventdata, handles)
% hObject    handle to currenttime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of currenttime as text
%        str2double(get(hObject,'String')) returns contents of currenttime as a double

i=handles.minsamp;
pageval = floor(str2double(get(hObject,'String'))); % str2double is faster than str2num
pageval = max(pageval,1);

if handles.multcomp
    maxiwin=floor((handles.maxx-handles.winsize)/handles.winsize);
else
    maxiwin=ceil((nsamples(handles.Dmeg{1})/fsample(handles.Dmeg{1}))/handles.winsize);
end
pageval = min(pageval,maxiwin);
slidval = max((pageval-1)*handles.winsize,handles.winsize/2);
set(handles.slider1,'Value',slidval)
if handles.scoring
    set(handles.figure1,'CurrentAxes',handles.axes4)
    crc_hypnoplot(handles.axes4, ...
        handles,handles.score{3,handles.currentscore})
    set(handles.figure1,'CurrentAxes',handles.axes1)
end
mainplot(handles)


% To view another file
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

%--------------------------------------------------------------------------
%Filters for the EEG, EOG or EMG signals-----------------------------------
%--------------------------------------------------------------------------

% Edit the lowpass filter cutoff of EEG
function frqabv_Callback(hObject, eventdata, handles)
% hObject    handle to frqabv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frqabv as text
%        str2double(get(hObject,'String')) returns contents of frqabv as a double

frbe = str2double(get(hObject,'String'));
i = handles.minsamp;
if not(isnan(frbe)) && length(frbe)==1
    
    if frbe>fsample(handles.Dmeg{i})/2*.99
        handles.filter.other(2) = fsample(handles.Dmeg{i})/2*.99;
        set(hObject,'String',num2str(fsample(handles.Dmeg{i})/2));
    else
        handles.filter.other(2) = frbe;
    end
    
else
    beep
    set(hObject,'String',num2str(handles.filter.other(2)))
    return
end
if fsample(handles.Dmeg{i})>1500
    forder = 1;
else
    forder = 3;
end
[B,A] = butter(forder,[handles.filter.other(1)/(fsample(handles.Dmeg{i})/2),...
    handles.filter.other(2)/(fsample(handles.Dmeg{i})/2)],'pass');
handles.filter.coeffother=[B;A];

% Plot
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


function settings_Callback(hObject, eventdata, handles)
% Update handles structure
handles.gui_active = 0;
guidata(hObject, handles);
h = crc_z3score_settings();
uiwait(h);
handles.gui_active = 1;
[path, ~,~] = fileparts([mfilename('fullpath') '.m']);
settings_path = fullfile(path,'/settings.mat');

if exist(settings_path, 'file') == 2,
    load(settings_path);
    set(handles.validity,'string',['License Validity: ' settings.validity ' UTC.'])
else
    set(handles.validity,'string','License Validity: NA')
end

guidata(hObject, handles);


function gen_report_Callback(hObject, eventdata, handles)
handles.gui_active = 0;
guidata(hObject, handles);
flags.Dmeg=handles.Dmeg;
flags.file=handles.file;
flags.scoresleep=1;
flags.index = handles.index;
flags.currentscore = handles.currentscore;
crc_z3score_report(handles.Dmeg{1}, flags);
close(handles.figure1);
return


function autoscore_Callback(hObject, eventdata, handles)
answer = questdlg('Scoring using Z3Score-FASST is not recommended. FASST can now import scores from NEO commandline directly (press Ctrl+I). Neo-commandline supports batch scoring, scoring using reduced channels, dynamic handling of noisy epochs. Find out more at z3score.com or contact support@neurobit.io.', ...
	'Score using Neo-commandline', ...
	'Continue','Cancel','Cancel');

if strcmp(answer, 'Cancel')
    return
end

handles.gui_active = 0;
guidata(hObject, handles);
flags.Dmeg=handles.Dmeg;
flags.file=handles.file;
flags.scoresleep=1;
flags.index = handles.index;
crc_z3score_score(handles.Dmeg{1}, flags);
close(handles.figure1);
return


%Edit the highpass filter cutoff
function frqblw_Callback(hObject, eventdata, handles)
% hObject    handle to frqblw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frqblw as text
%        str2double(get(hObject,'String')) returns contents of frqblw as a double

frbe    =   str2double(get(hObject,'String'));
i       =   handles.minsamp;
if not(isnan(frbe)) && length(frbe)==1
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
if fsample(handles.Dmeg{i})>1500
    forder = 1;
else
    forder = 3;
end
[B,A] = butter(forder,[handles.filter.other(1)/(fsample(handles.Dmeg{i})/2),...
    handles.filter.other(2)/(fsample(handles.Dmeg{i})/2)],'pass');
handles.filter.coeffother=[B;A];
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


function upemg_Callback(hObject, eventdata, handles)
% hObject    handle to upemg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of upemg as text
%        str2double(get(hObject,'String')) returns contents of upemg as a double

frbe	=   str2double(get(hObject,'String'));
i   =   handles.minsamp;
if not(isnan(frbe)) && length(frbe)==1
    if frbe>fsample(handles.Dmeg{i})/2*.99
        handles.filter.EMG(2)=fsample(handles.Dmeg{i})/2*.99;
        set(hObject,'String',num2str(fsample(handles.Dmeg{i})/2));
    else
        handles.filter.EMG(2) = frbe;
    end
else
    beep
    set(hObject,'String',num2str(handles.filter.EMG(2)))
    return
end
if fsample(handles.Dmeg{i})>1500
    forder = 1;
else
    forder = 3;
end
[B,A] = butter(forder,[handles.filter.EMG(1)/(fsample(handles.Dmeg{i})/2),...
    handles.filter.EMG(2)/(fsample(handles.Dmeg{i})/2)],'pass');
handles.filter.coeffEMG=[B;A];

guidata(hObject, handles);

% Main plot
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
i=handles.minsamp;
if not(isnan(frbe)) && length(frbe)==1
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
if fsample(handles.Dmeg{i})>1500
    forder = 1;
else
    forder = 3;
end
[B,A] = butter(forder,[handles.filter.EMG(1)/(fsample(handles.Dmeg{i})/2),...
    handles.filter.EMG(2)/(fsample(handles.Dmeg{i})/2)],'pass');
handles.filter.coeffEMG=[B;A];

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

%Edit the highpass cutoff frequency for EOG signal
function downeog_Callback(hObject, eventdata, handles)
% hObject    handle to downeog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of downeog as text
%        str2double(get(hObject,'String')) returns contents of downeog as a double

frbe=str2double(get(hObject,'String'));
i=handles.minsamp;
if not(isnan(frbe)) && length(frbe)==1
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
if fsample(handles.Dmeg{i})>1500
    forder = 1;
else
    forder = 3;
end
[B,A] = butter(forder,[handles.filter.EOG(1)/(fsample(handles.Dmeg{i})/2),...
    handles.filter.EOG(2)/(fsample(handles.Dmeg{i})/2)],'pass');
handles.filter.coeffEOG=[B;A];
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

%Edit the lowpass cutoff frequency for EOG signal
function upeog_Callback(hObject, eventdata, handles)
% hObject    handle to upeog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of upeog as text
%        str2double(get(hObject,'String')) returns contents of upeog as a double
frbe = str2double(get(hObject,'String'));
i = handles.minsamp;
if not(isnan(frbe)) && length(frbe)==1
    if frbe>fsample(handles.Dmeg{i})/2*.99
        handles.filter.EOG(2) = fsample(handles.Dmeg{i})/2*.99;
        set(hObject,'String',num2str(fsample(handles.Dmeg{i})/2));
    else
        handles.filter.EOG(2) = frbe;
    end
else
    beep
    set(hObject,'String',num2str(handles.filter.EOG(2)))
    return
end
if fsample(handles.Dmeg{i})>1500
    forder = 1;
else
    forder = 3;
end
[B,A] = butter(forder,[handles.filter.EOG(1)/(fsample(handles.Dmeg{i})/2),...
    handles.filter.EOG(2)/(fsample(handles.Dmeg{i})/2)],'pass');
handles.filter.coeffEOG=[B;A];
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



%--------------------------------------------------------------------------
%Menus to change the reference in either EEG, EMG or EEG
%--------------------------------------------------------------------------
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


%--------------------------------------------------------------------------
%Compute the power spectrum on one or all channels when right click--------
%--------------------------------------------------------------------------

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

%Get which file is selected
hand    =   get(handles.axes1,'Children');
Mouse   =   get(handles.axes1,'CurrentPoint');
Chan    =   ceil((Mouse(1,2)-handles.scale/2)/handles.scale);
slidval =   get(handles.slider1,'Value');

if handles.figz~=0
    z   =   handles.figz;
else
    z   =   figure;
end
figure(z)
axs     =   get(handles.figz,'Children');
cleargraph(handles.figz) %change made

if handles.multcomp
    fil     =   min(max(1,Chan),length(handles.Dmeg));
    Ctodis  =	handles.Chantodis;
    [dumb1,dumb2,Chan]  =   intersect(handles.chanset{Ctodis}, ...
        upper(chanlabels(handles.Dmeg{fil})));
    start   =   datevec(handles.date(fil,1)-handles.mindate);
    start   =   start(4)*60^2+start(5)*60+start(6);
    beg     =   slidval - start;
    tdeb    =   round(beg*fsample(handles.Dmeg{fil}));
    temps   =   tdeb:1:min(tdeb+(fsample(handles.Dmeg{fil})*handles.winsize), ...
        nsamples(handles.Dmeg{fil}));
    toshow  =   temps;
    cmap 	=   hsv(length(handles.Dmeg));
    Col     =   fil;
else
    fil     =   1;
    Chan    =   min(max(1,Chan),length(handles.indexMEEG) + length(handles.indnomeeg));
    NbreChandisp    =   str2double(get(handles.NbreChan,'String'));
    Chanslidval     =   get(handles.Chanslider,'Value');
    slidpos         =   Chanslidval-rem(Chanslidval,1);
    index           =   [handles.indnomeeg handles.indexMEEG];
    handles.inddis  =   index(slidpos : 1 : slidpos + NbreChandisp -1);
    Chan            =   index(Chan);
    chandis         =   chanlabels(handles.Dmeg{fil},Chan);
    fs              =   fsample(handles.Dmeg{fil});
    tdeb            =   round(slidval*fs);
    temps           =   tdeb:1:min(tdeb+(fs*handles.winsize), ...
        nsamples(handles.Dmeg{fil}));
    toshow          =   temps;
    cmap            =   [0.2 0.9 0.5; 1 0 0; 0 0 1 ];
end

tdeb_w = round(slidval*fs);
tend_w = min(tdeb+(fs*handles.winsize), ...
    nsamples(handles.Dmeg{fil}));
tohid_all = [];
if ~isempty(handles.score{5,handles.currentscore})
    tdebs   =   str2double(get(handles.currenttime,'String')) - handles.winsize/2;
    tfins   =   str2double(get(handles.currenttime,'String')) + handles.winsize/2;
    art     =   find(and((or(and(handles.score{5,handles.currentscore}(:,2)>tdebs,handles.score{5,handles.currentscore}(:,2)<tfins),...
        and(handles.score{5,handles.currentscore}(:,1)<tfins,handles.score{5,handles.currentscore}(:,1)>tdebs))),...
        or(handles.score{5,handles.currentscore}(:,3) == 0,handles.score{5,handles.currentscore}(:,3) == Chan)));
    art_concerned   =   handles.score{5,handles.currentscore}(art,1:2);
    if ~isempty(art_concerned)
        a=1;
        while a <= size(art_concerned,1)
            art_concerned(a);
            begart      =   max(tdeb_w,round(art_concerned(a,1)*fs));
            endart      =   min(tend_w,round(art_concerned(a,2)*fs));
            tohid2   	=   begart : endart;
            tohid_all   =   union(tohid_all,tohid2);
            tdeb_w      =   endart;
            a   =  a+1;
        end
        toshow 	=   setdiff(toshow,tohid_all);
    end
end
if length(toshow)<(handles.winsize/2)*fs        %More than 50% of artefacts
    h = msgbox('There is too much artefact on this channel');
    close(figure(z))
else
    fs      =   fsample(handles.Dmeg{fil});
    leg     =   cell(0);
    hold on
    [dumb1,dumb2,index2] = ...
        intersect(upper(chanlabels(handles.Dmeg{fil},Chan)),handles.names);
    if abs(handles.crc_types(index2))>1
        if handles.crc_types(index2)>0
            [dumb1,index1,dumb2] = ...
                intersect(upper(chanlabels(handles.Dmeg{fil})), ...
                upper(handles.names(handles.crc_types(index2))));
            try
                X   =   handles.Dmeg{fil}(Chan,toshow) - ...
                    handles.Dmeg{fil}(index1,toshow);
                Col	= 1;
            catch
                X   = 0;
            end
        else
            range   =   max(handles.Dmeg{fil}(Chan,toshow)) - ...
                min(handles.Dmeg{fil}(Chan,toshow));
            try
                X   = 	(handles.scale)*handles.Dmeg{fil}(Chan,toshow)/range;
                Col =   2;
            catch
                X   =   0;
            end
        end
    else
        try
            X   =   handles.Dmeg{fil}(Chan,toshow);
            Col	=   3;
        catch
            X   =   0;
        end
    end
    if length(X) == 1
        text(0.75,1, 'No Signal here')
        xlim([0 2])
        ylim([0 2])
        grid off
    else
        X       =   filterforspect(handles,X,[0.001 fs/3],fil);
        [P,F]   =   pwelch(X,[],[],[],fs);
        P       =   log(P);
        plot(F,P,'Color',cmap(Col,:))
        grid on
        titre   =   chandis;
        title(titre)
        ylabel('Log of power')
        xlabel('Frequency in Hz')
        [dumb name]     =   fileparts(handles.Dmeg{fil}.fname);
        under           =   find(name=='_');
        name(under)     =   ' ';
        leg{length(leg)+1}  =	name;
        legend(leg);
        axis auto
        xlim([0 20])
    end
end

handles.figz=z;

% Update handles structure
guidata(hObject, handles);


%Compute power spectrum on all files if multiple files comparison
% --------------------------------------------------------------------
function Cmp_Pwr_Sp_All_Callback(hObject, eventdata, handles)
% hObject    handle to Cmp_Pwr_Sp_All (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cmap = hsv(length(handles.Dmeg));

if handles.figz~=0
    z=handles.figz;
else
    z=figure;
end
figure(z)
axs=get(handles.figz,'Children');
cleargraph(z)

leg=cell(0);
for ii=1:size(cmap,1)
    fs=fsample(handles.Dmeg{ii});
    slidval = get(handles.slider1,'Value');
    start = datevec(handles.date(ii,1)-handles.mindate);
    start = start(4)*60^2+start(5)*60+start(6);
    beg = slidval - start;
    tdeb=round(beg*fsample(handles.Dmeg{ii}));
    temps=tdeb:1:min(tdeb+(fs*handles.winsize), ...
        nsamples(handles.Dmeg{ii}));
    toshow=temps;
    temps=(temps)/fsample(handles.Dmeg{ii})+start;
    Ctodis = handles.Chantodis;
    [dumb1,dumb2,index] = intersect(handles.chanset{Ctodis}, ...
        upper(chanlabels(handles.Dmeg{ii})));
    hold on
    [dumb1,dumb2,index2] = ...
        intersect(upper(chanlabels(handles.Dmeg{ii},index)),handles.names);
    if abs(handles.crc_types(index2))>1
        if handles.crc_types(index2)>0
            [dumb1,index1,dumb2] = ...
                intersect(upper(chanlabels(handles.Dmeg{ii})), ...
                upper(handles.names(handles.crc_types(index2))));
            try
                X = handles.Dmeg{ii}(index,toshow)- ...
                    handles.Dmeg{ii}(index1,toshow);
            catch
                X = 0;
            end
        else
            range = max(handles.Dmeg{ii}(index,toshow))- ...
                min(handles.Dmeg{ii}(index,toshow));
            try
                X =(handles.scale)*handles.Dmeg{ii}(index,toshow)/range;
            catch
                X = 0;
            end
        end
    else
        try
            X = handles.Dmeg{ii}(index,toshow);
        catch
            X = 0;
        end
    end
    
    if length(X) == 1
    else
        [P,F] = pwelch(X,[],[],[],fs);
        P = log(P);
        hold on
        plot(F,P,'Color',cmap(ii,:))
        [dumb name] = fileparts(handles.Dmeg{ii}.fname);
        under=find(name=='_');
        name(under)=' ';
        leg{length(leg)+1}=name;
    end
end
grid on
titre = char(handles.chanset{Ctodis});
title(titre)
legend(leg)
ylabel('Log of power')
xlabel('Frequency in Hz')
xlim([0 20])

handles.figz=z;

% Update handles structure
guidata(hObject, handles);


%--------------------------------------------------------------------------
% Export in a matlab figure------------------------------------------------
%--------------------------------------------------------------------------
function Export_Callback(hObject, eventdata, handles)

% hObject    handle to Export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.export=1;
z=figure;
axs=axes;
set(z,'CurrentAxes',axs)
handles.currentfigure=z;
mainplot(handles)
ha=get(handles.currentfigure,'CurrentAxes');
%display y-labels
if handles.multcomp
    li = length(handles.Dmeg);
else
    Chanslidval     =   get(handles.Chanslider,'Value');
    slidpos         =   Chanslidval-rem(Chanslidval,1);
    NbreChandisp    =   str2double(get(handles.NbreChan,'String'));
    index           =   handles.index(slidpos:1:slidpos+NbreChandisp-1);
    li   =   length(index);
end
ylim([0 handles.scale*(li+1)])
%set(ha,'YTick',[handles.scale/2:handles.scale/2:li*handles.scale+handles.scale/2]);
%ylabels     =   [num2str(round(handles.scale/2))];
for j=1:li
    if handles.multcomp
        ylabels=[ylabels num2str(j)];
    else
        ylabels=[ylabels chanlabels(handles.Dmeg{1},index(j))];
    end
    ylabels=[ylabels num2str(round(handles.scale/2))];
end
%set(ha,'YTickLabel',ylabels);

%display x-labels
slidval=get(handles.slider1,'Value');
xlim([slidval slidval+handles.winsize])
xtick = get(handles.axes1,'XTick');
if isfield(handles,'offset')
    xtick = mod(xtick + handles.offset,24*60^2);
end
[time string] = crc_time_converts(xtick);
set(ha,'XTickLabel',string)

handles.figz=z;
handles.export=0;

% Update handles structure
guidata(hObject, handles);

return

%--------------------------------------------------------------------------
%Management of number of channels and slider-------------------------------
%--------------------------------------------------------------------------

%Slider for the channels when only one file
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


%Edit the number of channel to display
function NbreChan_Callback(hObject, eventdata, handles)

% hObject    handle to NbreChan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NbreChan as text
%        str2double(get(hObject,'String')) returns contents of NbreChan as a double

Nchdispmeeg = str2double(get(handles.NbreChan,'String'));
if ~isnan(Nchdispmeeg) && length(Nchdispmeeg)==1
    set(handles.NbreChan,'String', num2str(min(Nchdispmeeg, length(handles.indexMEEG) + length(handles.indnomeeg))));
    Nchdisp = str2double(get(handles.NbreChan,'String'));
    totind  = length(handles.indnomeeg) + length(handles.indexMEEG);
    set(handles.Chanslider,...
        'Min',1,...
        'Max', totind - Nchdisp+1,...
        'Value',1,...
        'SliderStep', [Nchdisp/totind Nchdisp/totind]);
else
    beep
    set(handles.NbreChan,'String', num2str(10))
end
mainplot(handles);
guidata(hObject,handles);

%---the MEGPLANAR channels normalized or not (to be checked)%%%%%%

function normalize_Callback(hObject, eventdata, handles)

% hObject    handle to normalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns the value of normalize

norm = get(handles.normalize,'Value');
indexMEGPLAN = intersect(meegchannels(handles.Dmeg{1}, 'MEGPLANAR'), handles.index);
if norm
    n = 0;
    S = struct(handles.Dmeg{1});
    handles.indexnorm = [];
    handles.nind = [];
    while n <=length(indexMEGPLAN)-1
        n = n+1;
        while (n<=length(indexMEGPLAN)-1)&&(S.channels(indexMEGPLAN(n)).X_plot2D == S.channels(indexMEGPLAN(n+1)).X_plot2D)&(S.channels(indexMEGPLAN(n)).Y_plot2D == S.channels(indexMEGPLAN(n+1)).Y_plot2D)
            handles.indexnorm = [handles.indexnorm indexMEGPLAN(n)];
            handles.nind = [handles.nind indexMEGPLAN(n+1)];
            n = n+2;
        end
        if n<=length(indexMEGPLAN)
            handles.nind = [handles.nind indexMEGPLAN(n)];
        end
        
    end
    [name ind] = intersect(handles.indexMEEG, handles.nind);
    while ~isempty(ind)
        handles.indexMEEG = [handles.indexMEEG(1:ind(1)-1) handles.indexMEEG(ind(1) + 1 : length(handles.indexMEEG))];
        ind = ind(2:end);
    end
else
    handles.indexMEEG = handles.index(length(handles.indnomeeg)+1:numel(handles.index));
end
%Update the number of channels available
NbreChan    =   length(handles.indexMEEG) + length(handles.indnomeeg);
Nchdisp     =   min(str2double(get(handles.NbreChan,'String')),NbreChan);
handles.inddis = [handles.indnomeeg handles.indexMEEG];

set(handles.NbreChan,'String',Nchdisp);
totind  =   length(handles.indexMEEG)+length(handles.indnomeeg);
set(handles.Chanslider,...
    'Min',1,...
    'Max',totind - Nchdisp+1,...
    'Value',1,...
    'SliderStep', [Nchdisp/totind Nchdisp/totind]);
guidata(hObject,handles);
mainplot(handles);

%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%end to be checked%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.

function NbreChan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NbreChanEeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%--------------------------------------------------------------------------
%Management of the scales -------------------------------------------------
%--------------------------------------------------------------------------

%Menu of scales
% --- Executes on selection change in menu scales.
function scales_Callback(hObject, eventdata, handles)
% hObject    handle to scales (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%EEG scale
% --- Executes on selection change in menu scales.
function scale_eeg_Callback(hObject, eventdata, handles)
% hObject    handle to scale_eeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
eegsc=spm_input('EEG scale?',1,'s','1');
close gcf
try
    handles.eegscale=eval(eegsc);
catch
    if strcmpi(eegsc,'V')
        handles.eegscale=10^-6;
    elseif strcmpi(eegsc,'µV')
        handles.eegscale=1;
    else
        beep
        disp('Enter scale in 10^-x format or V or µV')
        handles.eegscale=[];
        return
    end
end
mainplot(handles)
guidata(hObject, handles);

%EOG scale
% --- Executes on selection change in menu scales.
function scale_eog_Callback(hObject, eventdata, handles)
% hObject    handle to scale_eeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
eogsc=spm_input('EOG scale?',1,'s','1');
close gcf
try
    handles.eogscale=eval(eogsc);
catch
    if strcmpi(eogsc,'V')
        handles.eogscale=10^-6;
    elseif strcmpi(eogsc,'µV')
        handles.eogscale=1;
    else
        beep
        disp('Enter scale in 10^-x format or V or µV')
        handles.eogscale=[];
        return
    end
end
mainplot(handles)
guidata(hObject, handles);

%EMG scale
% --- Executes on selection change in menu scales.
function scale_emg_Callback(hObject, eventdata, handles)
% hObject    handle to scale_eeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
emgsc=spm_input('EMG scale?',1,'s','1');
close gcf
try
    handles.emgscale=eval(emgsc);
catch
    if strcmpi(emgsc,'V')
        handles.emgscale=10^-6;
    elseif strcmpi(emgsc,'µV')
        handles.emgscale=1;
    else
        beep
        disp('Enter scale in 10^-x format or V or µV')
        handles.emgscale=[];
        return
    end
end
mainplot(handles)
guidata(hObject, handles);

%ECG scale
% --- Executes on selection change in menu scales.
function scale_ecg_Callback(hObject, eventdata, handles)
% hObject    handle to scale_eeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ecgsc=spm_input('ECG scale?',1,'s','1');
close gcf
try
    handles.ecgscale=eval(ecgsc);
catch
    if strcmpi(ecgsc,'V')
        handles.ecgscale=10^-6;
    elseif strcmpi(ecgsc,'µV')
        handles.ecgscale=1;
    else
        beep
        disp('Enter scale in 10^-x format or V or µV')
        handles.ecgscale=[];
        return
    end
end
mainplot(handles)
guidata(hObject, handles);

%LFP scale
% --- Executes on selection change in menu scales.
function scale_lfp_Callback(hObject, eventdata, handles)
% hObject    handle to scale_eeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ecgsc=spm_input('LFP scale?',1,'s','1');
close gcf
try
    handles.lfpscale=eval(ecgsc);
catch
    if strcmpi(ecgsc,'V')
        handles.lfpscale=10^-6;
    elseif strcmpi(ecgsc,'µV')
        handles.lfpscale=1;
    else
        beep
        disp('Enter scale in 10^-x format or V or µV')
        handles.lfpscale=[];
        return
    end
end
mainplot(handles)
guidata(hObject, handles);

%Other scale
% --- Executes on selection change in menu scales.
function scale_other_Callback(hObject, eventdata, handles)
% hObject    handle to scale_eeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ecgsc=spm_input('Other scale?',1,'s','1');
close gcf
try
    handles.otherscale=eval(ecgsc);
catch
    beep
    disp('Enter scale in 10^-x format')
    handles.otherscale=[];
    return
end
mainplot(handles)
guidata(hObject, handles);


%Menu of MEG scales
% --- Executes on selection change in menu scales.
function scale_meg_Callback(hObject, eventdata, handles)
% hObject    handle to scales (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%ECG scale
% --- Executes on selection change in menu scales.
function scale_megmag_Callback(hObject, eventdata, handles)
% hObject    handle to scale_eeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ecgsc=spm_input('Magnometers scale?',1,'s','1');
close gcf
try
    handles.megmagscale=eval(ecgsc);
catch
    if strcmpi(ecgsc,'T')
        handles.megmagscale=10^-12;
    elseif strcmpi(ecgsc,'fT')
        handles.megmagscale=1;
    else
        beep
        disp('Enter scale in 10^-x format or T or fT')
        handles.megmagscale=[];
        return
    end
end
mainplot(handles)
guidata(hObject, handles);

%ECG scale
% --- Executes on selection change in menu scales.
function scale_megplanar_Callback(hObject, eventdata, handles)
% hObject    handle to scale_eeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ecgsc=spm_input('Gradiometers scale?',1,'s','1');
close gcf
try
    handles.megplanarscale=eval(ecgsc);
catch
    if strcmpi(ecgsc,'T/m')
        handles.megplanarscale=10^-11;
    elseif strcmpi(ecgsc,'fT/m')
        handles.megplanarscale=1;
    else
        beep
        disp('Enter scale in 10^-x format or T/m or fT/m')
        handles.megplanarscale=[];
        return
    end
end
mainplot(handles)
guidata(hObject, handles);



%--------------------------------------------------------------------------
%Display of events and travel in the data using types of events------------
%--------------------------------------------------------------------------

%Select the type of event
% --- Executes on selection change in popupmenu10.
function popupmenu10_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu10 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu10

evtype = get(handles.popupmenu10,'String');
evnum  = get(handles.popupmenu10,'Value');

if evnum==1
    evnum = 2 : size(evtype,1);
end

chos=[];

for i=1:max(size(handles.evt,2),size(handles.evt,1))
    if any(strcmpi(handles.evt(i).type,evtype(evnum)))
        chos=[chos, i];
    end
end

handles.chostype=chos;
pmstring=[{'All'},{'Scan number'}];
nsc = handles.Dmeg{1}.CRC.score{2,:};
for i=1:size(chos,2)
    if and(~strcmpi(handles.evt(chos(i)).value, pmstring),any(strcmpi(handles.evt(chos(i)).value,nsc)))
        pmstring = [pmstring, {handles.evt(chos(i)).value}];
    end
end

handles.chosevt     = handles.chostype;
handles.displevt    = 0;

set(handles.popupmenu11,...
    'String',pmstring,...
    'Value',1) %pmstring in old version... to be checked

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%Goes to the previous event
% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.displevt
    ind=find(handles.displevt==handles.chosevt);
    if ind-1>0
        handles.displevt=handles.chosevt(ind-1);
    else
        handles.displevt=handles.chosevt(ind);
    end
else
    difdur=zeros(1,size(handles.chosevt,2));
    slidval = get(handles.slider1,'Value');
    for i=1:size(handles.chosevt,2)
        difdur(i)=handles.evt(handles.chosevt(i)).time-slidval;
    end
    difdur(difdur>=0)=0;
    difdur(difdur<0)=1;
    diffdur=diff(difdur);
    todisp=find(diffdur==-1);
    if ~isempty(todisp)
        handles.displevt=handles.chosevt(todisp-1);
    elseif sum(difdur)==0
        handles.displevt=handles.chosevt(1);
    else
        handles.displevt=handles.chosevt(end);
    end
end
stat=handles.Dmeg{1}.CRC.goodevents(handles.displevt);
if stat==0
    set(handles.radiobutton1,'Value',0);
else
    set(handles.radiobutton1,'Value',1);
end
slidval=handles.evt(handles.displevt).time-handles.winsize/2;
if slidval<=1
    slidval=handles.evt(handles.displevt).time-2; %default value if first event is at the boundary of the file
end
set(handles.slider1,'Value',slidval)
if handles.scoring
    set(handles.figure1,'CurrentAxes',handles.axes4)
    crc_hypnoplot(handles.axes4, ...
        handles,handles.score{3,handles.currentscore})
    set(handles.figure1,'CurrentAxes',handles.axes1)
end
a=events(handles.Dmeg{1});
if or(strcmpi(a(handles.displevt).type,'SW'), ...
        strcmpi(a(handles.displevt).type,'delta'))
    set(handles.delaymap,'visible','on');
    set(handles.text25,'visible','on');
end
mainplot(handles);
delaymap_Callback(hObject, eventdata, handles);
guidata(hObject, handles);



%Goes to next event
% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.displevt
    ind=find(handles.displevt==handles.chosevt);
    if ind+1<=size(handles.chosevt,2)
        handles.displevt=handles.chosevt(ind+1);
        if handles.evt(handles.chosevt(ind+1)).time>= nsamples(handles.Dmeg{1})...
                /fsample(handles.Dmeg{1})
            handles.displevt=handles.chosevt(ind);
        end
    else
        handles.displevt=handles.chosevt(ind);
    end
else
    difdur=zeros(1,size(handles.chosevt,2));
    slidval = get(handles.slider1,'Value');
    for i=1:size(handles.chosevt,2)
        difdur(i)=handles.evt(handles.chosevt(i)).time-slidval;
    end
    difdur(difdur>=0)=0;
    difdur(difdur<0)=1;
    diffdur=diff(difdur);
    todisp=find(diffdur==-1);
    if ~isempty(todisp)&& todisp+1<size(handles.evt,2)
        handles.displevt=handles.chosevt(todisp+1);
    elseif sum(difdur)==0
        handles.displevt=handles.chosevt(1);
    else
        handles.displevt=handles.chosevt(end);
    end
end
stat=handles.Dmeg{1}.CRC.goodevents(handles.displevt);
if stat==0
    set(handles.radiobutton1,'Value',0);
else
    set(handles.radiobutton1,'Value',1);
end
slidval=handles.evt(handles.displevt).time-handles.winsize/2;
if slidval<=0
    slidval=handles.evt(handles.displevt).time-2; %default value if first event is at the boundary of the file
end
set(handles.slider1,'Value',slidval)
if handles.scoring
    set(handles.figure1,'CurrentAxes',handles.axes4)
    crc_hypnoplot(handles.axes4, ...
        handles,handles.score{3,handles.currentscore})
    set(handles.figure1,'CurrentAxes',handles.axes1)
end
a=events(handles.Dmeg{1});
if or(strcmpi(a(handles.displevt).type,'SW'), ...
        strcmpi(a(handles.displevt).type,'delta'))
    set(handles.delaymap,'visible','on');
    set(handles.text25,'visible','on');
end
mainplot(handles);
delaymap_Callback(hObject, eventdata, handles);
guidata(hObject, handles);

%Check to say if event is good (1) or bad (0)
% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1
goodevt=get(handles.radiobutton1,'Value');
D=handles.Dmeg{1};
aevt=events(D);
gdevt=D.CRC.goodevents;
if ~goodevt
    gdevt(handles.displevt)=0;
else
    gdevt(handles.displevt)=1;
end
if or(strcmpi(aevt(handles.displevt).type,'SW'), ...
        strcmpi(aevt(handles.displevt).type,'delta'))
    a=aevt(handles.displevt).time;
    [d1,b]=sort(abs(round([D.CRC.SW.SW(:).negmax]/1000)-round(a)));
    D.CRC.SW.SW(b(1)).good=goodevt;
elseif or(strcmpi(aevt(handles.displevt).type,'Post-SP'), ...
        strcmpi(aevt(handles.displevt).type,'Ant-SP')) || ...
        strcmpi(aevt(handles.displevt).type,'SP')
    a=aevt(handles.displevt).time;
    [d1,b]=sort(abs((round(D.CRC.spindles.bounds(:,1)/fsample(handles.Dmeg{1})*10)-round(a*10))));
    D.CRC.spindles.good(b(1))=goodevt;
end
D.CRC.goodevents=gdevt;
handles.Dmeg{1}=D;
save(D);
guidata(hObject,handles)

%Create the delay map of a SW if asked and when scrolling through events
% --- Executes on button press in delaymap.
function delaymap_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dispmap=get(handles.delaymap,'Value');
hand=handles;
if dispmap
    D=handles.Dmeg{1};
    aevt=events(D);
    if or(strcmpi(aevt(handles.displevt).type,'SW'), ...
            strcmpi(aevt(handles.displevt).type,'delta'))
        a=aevt(handles.displevt).time;
        asw=D.CRC.SW.SW;
        [d1,b]=sort(abs(round([asw(:).negmax])-round(a*1000)));
        toshow=b(1);
        if ~isfield(handles.delmap,'elpos')
            handles.delmap.elpos=spm_input('Electrodes positions:',+1,'b','file|auto',[0,1],0);
            close gcf
            if handles.delmap.elpos==0 && ...
                    (~isfield(handles.delmap,'zebris_name') || ...
                    isempty(handles.delmap,'zebris_name'))
                handles.delmap.zebris_name=spm_select(1, 'any', 'Select electrodes positioning file','' ,pwd,'.*');
            else
                handles.delmap.zebris_name='CRC_electrodes.mat';
            end
        end
        if ~isfield(handles.delmap,'cas')
            handles.delmap.cas=1;
        end
        hand.delmap.elpos=handles.delmap.elpos;
        hand.delmap.zebris_name=handles.delmap.zebris_name;
        hand.delmap.cas=handles.delmap.cas;
        handles=hand;
        set(handles.figure1,'CurrentAxes',handles.axes5);
        delete(get(handles.axes5,'Children'));
        fig=handles.figure1;
        ax=handles.axes5;
        [handles.delmap.hb,handles.delmap.titl]=crc_SWS_mapping(handles.Dmeg{1},toshow,handles.delmap.cas,handles.delmap.zebris_name,handles.delmap.elpos,fig,ax);
    end
else
    delete(get(handles.axes5,'Children'));
    if isfield(handles.delmap,'hb')
        try
            delete(handles.delmap.hb);
        end
    end
    if isfield(handles.delmap,'titl')
        try
            delete(handles.delmap.titl);
        end
    end
end
%update handles
guidata(hObject,handles)


%Menu to chose the value of the event within a selected type
% --- Executes on selection change in popupmenu11.
function popupmenu11_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu11 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu11
evval=get(handles.popupmenu11,'String');
evnu=get(handles.popupmenu11,'Value');
if evnu==1 || evnu==2
    evnu=3:size(evval,1);
end
chos=[];
for i=1:size(handles.chostype,2)
    if any(strcmpi(handles.evt(handles.chostype(i)).value,evval(evnu)))
        chos=[chos, handles.chostype(i)];
    end
end
handles.chosevt=chos;
handles.displevt=0;

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
%Save the position of the slider in the data-------------------------------
%--------------------------------------------------------------------------
% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

slidval=get(handles.slider1,'Value');
D=handles.Dmeg{1};
if isfield(handles.Dmeg{1},'CRC')
    D.CRC.lastdisp=slidval;
else
    D.CRC=struct('lastdisp',[]);
    D.CRC.lastdisp=slidval;
end
handles.Dmeg{1}=D;
save(D);



%--------------------------------------------------------------------------
%----------------------- Scoring options ----------------------------------
%--------------------------------------------------------------------------

% --- Executes on button press in Score.
function Score_Callback(hObject, eventdata, handles)
% hObject    handle to Score (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Handling score, arousal, artefact, opl/fpl information, window size & events of interest info
handles.scoring=1;
try
    if size(handles.Dmeg{1}.CRC.score,1)<7
        % Meaning it is an old cell array
        
        handles.score = cell(8,size(handles.Dmeg{1}.CRC.score,2));
        
        %Taking the score & username from the old cell array
        handles.score(1:2,:)=handles.Dmeg{1}.CRC.score(1:2,:);
        
        for ii=1:size(handles.Dmeg{1}.CRC.score,2)
            % Old cell array had a 20 seconds windows size
            handles.score{3,ii} = 20;
            % Updating OPL & FPL
            try
                handles.score{4,ii} = handles.Dmeg{1}.CRC.pl;
            catch
                handles.score{4,ii} = [1/fsample(handles.Dmeg{1}) ...
                    nsamples(handles.Dmeg{1})/fsample(handles.Dmeg{1})];
            end
            % Updating Artefact info
            try
                handles.score{5,ii}= handles.Dmeg{1}.CRC.art;
            catch
                handles.score{5,ii}=[];
            end
            % Updating Arousal info
            try
                handles.score{6,ii}=handles.Dmeg{1}.CRC.art;
            catch
                handles.score{6,ii}=[];
            end
            
            % Creating event of interest vector
            handles.score{7,ii}=[];
            %vector for names of artefacts
            handles.score{8,ii}=[];
        end
    else
        handles.score = handles.Dmeg{1}.CRC.score;
        if size(handles.score,1)==7
            for isc=1:size(handles.score,2)
                handles.score{8,isc}=cell(size(handles.score{5,isc},1),1);
            end
        end
        
    end
    handles.Dmeg{1}.CRC.score = handles.score;
    handles.currentscore = 1;
    if(handles.showOracle)
        for isc=1:size(handles.score,2)
            if strcmpi('Oracle',handles.score{2,isc})
                handles.currentscore=isc;
                break;
            end
        end
    end
    
    if ~isempty (handles.score{5,handles.currentscore})&&size(handles.score{5,handles.currentscore},2)<3
        handles.score{5,1}(:,3)=0;
    end
catch
    %Creating from scratch the right structure
    handles.score=cell(8,1);
    %%%%%%%%%%%%%%%%
    % From 1 to 7: %
    %%%%%%%%%%%%%%%%
    % 1: Score
    % 2: Name of scorer
    % 3: window size chosen to score
    % 4: OPL&FPL
    % 5: Artefacts
    % 6: Arousals
    % 7: Events of Interests
    % 8: Labels of the artefacts (empty if unspecified)
    
    %Defining user Name
    prompt = {'Please enter your name'};
    def= {'Newuser'};
    num_lines = 1;
    dlg_title = 'Name of the new scorer';
    handles.score(2,1) = inputdlg(prompt,dlg_title,num_lines,def);
    
    %Choosing size of window to score
    prompt = {'Please choose the size of the scoring windows'};
    def= {'30'};
    num_lines = 1;
    dlg_title = 'Size of the scoring windows (in sec)';
    handles.score(3,1) = inputdlg(prompt,dlg_title,num_lines,def);
    handles.score{3,1} = str2double(handles.score{3});
    
    % Creating FPL & OPL
    handles.score{4,1} = [1/fsample(handles.Dmeg{1}) ...
        nsamples(handles.Dmeg{1})/fsample(handles.Dmeg{1})];
    
    % Creating artefacts
    handles.score{5,1} = [];
    
    % Creating arousals
    handles.score{6,1} = [];
    
    % Creating event of interest
    handles.score{7,1} = [];
    
    % Creating labels associated with artefacts
    handles.score{8,1} = [];
    
    % Putting NaN in the score
    handles.score{1,1} = ...
        0/0*ones(1,ceil(nsamples(handles.Dmeg{1}) / ...
        (handles.score{3,1}*fsample(handles.Dmeg{1}))));
    if isfield(handles.Dmeg{1}, 'CRC')
        handles.Dmeg{1}.CRC.score=handles.score;
    else
        handles.Dmeg{1}.CRC=struct('score',[]);
        handles.Dmeg{1}.CRC.score=handles.score;
    end
    handles.currentscore=1;
end
if ~isempty (handles.score{5,1})&&size(handles.score{5,1},2)<3
    handles.score{5,1}(:,3)=0;
end
%check that last page has a score, if not put nan
sc=handles.Dmeg{1}.CRC.score;
for i=1:size(sc,2)
    theosz=ceil(nsamples(handles.Dmeg{1})/fsample(handles.Dmeg{1})/sc{3,i}); % Theosz contents the number of windows
    if theosz>size(sc{1,i},2)
        handles.Dmeg{1}.CRC.score{1,i}=[handles.Dmeg{1}.CRC.score{1,i}, 0/0];
    end
    save(handles.Dmeg{1});
end

% Setting up the add artefact/arousal/event of interest.
handles.adddeb = ones(1,size(handles.score,2));
handles.addardeb = ones(1,size(handles.score,2));
handles.add_eoi = ones(1,size(handles.score,2));
% Defining window size according the current scorer
handles.winsize=handles.score{3,handles.currentscore};

%Changing the visibility of axes and scoring menu
set(handles.axes4,'Visible','on')
set(handles.fftchan,'Visible','off') %to plot fft
set(handles.score_menu,'enable','on')
%set(handles.score_stats,'enable','on')
set(handles.score_import,'enable','on')
set(handles.score_compare,'enable','on')
%set(handles.score_check,'enable','on')
set(handles.score_user,'enable','on')
set(handles.del_user,'enable','on')
set(handles.vertgrid,'enable','on')
set(handles.horgrid,'enable','on')
set(handles.grid,'enable','on')
set(handles.addundefart,'enable','on','visible','on')
set(handles.addonlyone,'enable','on','visible','on') %artifact on single channel
set(handles.addspecart,'enable','on','visible','on')
set(handles.addaro,'enable','on','visible','on')
set(handles.addeoi,'enable','on','visible','on')
set(handles.FPL,'enable','on','visible','on')
set(handles.OPL,'enable','on','visible','on')
set(handles.Delart,'enable','on','visible','on')
set(handles.Delartonlyone,'enable','on','visible','on') %delete artifact on single channel
set(handles.Delaro,'enable','on','visible','on')
set(handles.Del_eoi,'enable','on','visible','on')
set(handles.manevent,'enable','on','visible','on'); %To put manual event

handles.namesc = [];
delete(get(handles.score_user,'Children'));
for isc=1:size(handles.score,2)
    handles.scorers{isc} = uimenu(handles.score_user,'Label', ...
        char(handles.score(2,isc)),'Callback',@defined_scorer) ;
    handles.namesc{isc}=char(handles.score(2,isc));
end

delete(get(handles.del_user,'Children'));
for isc=1:size(handles.score,2)
    handles.scorers{isc} = uimenu(handles.del_user,'Label', ...
        char(handles.score(2,isc)),'Callback',@delete_scorer) ;
    handles.namesc{isc}=char(handles.score(2,isc));
end

handles.num_scorers=isc;
handles.scorers{isc+1}=uimenu(handles.score_user,'Label', ...
    'New scorer','Callback',@new_scorer,...
    'Separator', 'on') ;
set(handles.scorers{handles.currentscore},'Checked','on');
delete(get(handles.addspecart,'Children'));
crcdef = crc_get_defaults('score');
for iart=1:size(crcdef.lab_art,2)
    uimenu(handles.addspecart,'Label', ...
        char(crcdef.lab_art{iart}),'Callback',{@addtypeart,handles}) ;
end

delete(get(handles.axes4,'Children'));
set(handles.figure1,'CurrentAxes',handles.axes4);
crc_hypnoplot(handles.axes4, ...
    handles,handles.score{3,handles.currentscore});
set(handles.figure1,'CurrentAxes',handles.axes1);
delete(get(handles.addonlyone,'Children'));
uimenu(handles.addonlyone,'Label', ...
    ('Start an undefined artefact on this channel'),'Callback',{@addstart,handles}) ;
set(handles.edit1,'String',num2str(handles.winsize)); % To update the size of the window according to the scorer used
% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------
%Menu for the events (to be checked)
%--------------------------------------------------------------------------

function event_menu_Callback (hObject,eventdata,handles)
% hObject    handle to event_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function Properties_Callback (hObject,eventdata,handles)
% hObject    handle to Properties (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.figure1, 'windowbuttonmotionfcn', '')

slidval=get(handles.slider1,'Value');
D = handles.Dmeg{1};
if isfield(handles.Dmeg{1},'CRC')
    D.CRC.lastdisp=slidval;
else
    D.CRC = struct('lastdisp',[]);
    D.CRC.lastdisp = slidval;
end
handles.Dmeg{1} = D;
save(D);

flags.Dmeg  =   handles.Dmeg;
flags.type  =   handles.type;
flags.file  =   handles.file;
flags.index =   handles.index;
flags.user  =   handles.currentscore;

crc_modif_events(flags);


% --------------------------------------------------------------------
function score_menu_Callback(hObject, eventdata, handles)
% hObject    handle to score_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function score_stats_Callback(hObject, eventdata, handles)
% hObject    handle to score_stats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

crcdef = crc_get_defaults('score');
legendm=cell(0);

% Define colormap
cmapref = [0.2 0.75 0.6; 0 0.8 1; 0.1 0.5 0.9; 0.1 0.2 0.8; 0.1 0.15 0.5; 0.5 0.5 0.9; 0.9 0.4 0.4; 0.9 0.6 0.3]; %the 8 colors for all sleep stages
cmap=[];

FPL = handles.score{4,handles.currentscore}(1);
OPL = handles.score{4,handles.currentscore}(2);

%Time allowed to sleep
TRS = OPL - FPL;
[TRStime TRSstring] = crc_time_converts(TRS);

%Invalidation
adapted=1:length(handles.Dmeg{1}.CRC.score{1,handles.currentscore});
nottobescored = find(adapted <...
    handles.score{4,handles.currentscore}(1)/handles.score{3,handles.currentscore} | ...
    adapted > ...
    handles.score{4,handles.currentscore}(2)/handles.score{3,handles.currentscore});

stat = cell(0);
iisc = handles.Dmeg{1}.CRC.score{1,handles.currentscore};

iisc(nottobescored)=-1;

Zero = find(iisc==0);   % awake
One = find(iisc==1);    % Stage 1
Two = find(iisc==2);    % Stage 2
Three = find(iisc==3);  % Stage 3
Four = find(iisc==4);   % Stage 4
Five = find(iisc==5);   % REM
Six = find(iisc==6);    % MT
Seven=find(iisc==7);    % Awake

%Time of the sleeping period
TPS = (max([Two Three Four Five])-1)*handles.score{3,handles.currentscore}...
    - (min([Two Three Four Five])-1)*handles.score{3,handles.currentscore};

[TPStime TPSstring] = crc_time_converts(TPS);

%Total Sleep Time
TST = length([Two Three Four Five])*handles.score{3,handles.currentscore};

[TSTtime TSTstring] = crc_time_converts(TST);

%Sleep Efficiency
SEff = TST/TRS;

%Latency St1
if isempty(One)
    LatS1string = ['No ',crcdef.stnames_L{2},' Scored'];
else
    LatS1 = (One(1)-1)*handles.score{3,handles.currentscore} - FPL;
    [LatS1time LatS1string] = crc_time_converts(LatS1);
end
%Lantency St2
if isempty(Two)
    LatS2string = ['No ',crcdef.stnames_L{3},' Scored'];
else
    LatS2 = (Two(1)-1)*handles.score{3,handles.currentscore} - FPL;
    [LatS2time LatS2string] = crc_time_converts(LatS2);
end

%Lantency REM
if isempty(Five)
    LatREMstring = ['No ',crcdef.stnames_L{6},' Scored'];
else
    LatREM = (Five(1)-1)*handles.score{3,handles.currentscore} - FPL;
    [LatREMtime LatREMstring] = crc_time_converts(LatREM);
end

%Min & Percentage of W
W = length([Zero])*handles.score{3,handles.currentscore};
[Wtime Wstring] = crc_time_converts(W);

PW = W/TST;
% if (PW ~= 0)
legendm = [legendm crcdef.stnames_S{1}]
cmap=[cmap;cmapref(1,:)];
% end

%Min & Percentage of S1
S1 = length([One])*handles.score{3,handles.currentscore};
[S1time S1string] = crc_time_converts(S1);
PS1 = S1/TST;
% if (PS1 ~= 0)
legendm = [legendm crcdef.stnames_S{2}];
cmap=[cmap;cmapref(2,:)];
% end

%Min & Percentage of S2
S2 = length([Two])*handles.score{3,handles.currentscore};
[S2time S2string] = crc_time_converts(S2);
PS2 = S2/TST;
% if (PS2 ~= 0)
legendm = [legendm crcdef.stnames_S{3}];
cmap=[cmap;cmapref(3,:)];
% end

%Min & Percentage of S3
S3 = length([Three])*handles.score{3,handles.currentscore};
[S3time S3string] = crc_time_converts(S3);
PS3 = S3/TST;
% if (PS3 ~= 0)
legendm = [legendm crcdef.stnames_S{4}];
cmap=[cmap;cmapref(4,:)];
% end

%Min & Percentage of S4
S4 = length([Four])*handles.score{3,handles.currentscore};
[S4time S4string] = crc_time_converts(S4);
PS4 = S4/TST;
% if (PS4 ~= 0)
legendm = [legendm crcdef.stnames_S{5}];
cmap=[cmap;cmapref(5,:)];
% end

%Min & Percentage of SWS
SWS = length([Three Four])*handles.score{3,handles.currentscore};
[SWStime SWSstring] = crc_time_converts(SWS);

PSWS = SWS/TST;

%Min & Percentage of REM
REM = length([Five])*handles.score{3,handles.currentscore};
[REMtime REMstring] = crc_time_converts(REM);
PREM = REM/TST;
% if (PREM ~= 0)
legendm = [legendm crcdef.stnames_S{6}];
cmap=[cmap;cmapref(6,:)];
% end

%Min & Percentage of MT
MT = length([Six])*handles.score{3,handles.currentscore};
[MTtime MTstring] = crc_time_converts(MT);
PMT = MT/TST;
%if (PMT ~= 0)
legendm = [legendm crcdef.stnames_S{7}];
cmap=[cmap;cmapref(7,:)];
%end

%percentage of 'unscorable'
Unsc = length([Seven])*handles.score{3,handles.currentscore};
[Unsctime Unscstring] = crc_time_converts(Unsc);
legendm = [legendm crcdef.stnames_S{8}];
cmap=[cmap;cmapref(8,:)];
PUnsc = Unsc/TST;

%Begining display

if handles.figz~=0
    fig_z = handles.figz;
else
    fig_z = 0;
end
% call stat & display routine
[out_V,fig_z] = crc_statgen(handles.Dmeg{1},handles.currentscore,fig_z);
handles.figz = fig_z;
% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function score_export_Callback(hObject, eventdata, handles)
% hObject    handle to score_import (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[~,name,~] = fileparts(fullfile(handles.Dmeg{1}));
[file, path, ~] = uiputfile([name '_' handles.Dmeg{1}.CRC.score{2,handles.currentscore} '_scores.csv']);
scores = handles.Dmeg{1}.CRC.score{1,handles.currentscore};
csvwrite([path file], scores');

% --------------------------------------------------------------------
function score_importneo_Callback(hObject, eventdata, handles)
% hObject    handle to score_import (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file,path,indx] = uigetfile('*.neo', '*.json');
if ~indx
    return;
end

try
    fileID = fopen([path file],'r');
    response = fscanf(fileID,'%s');
    fclose(fileID);
    response = loadjson(response);
    scores = response.message;
    
    if response.status == 0,
        close(h);
        errordlg(['Error message:' response.message],'Error importing data.'); 
        return;
    end

    
    %check if structure exits
    if(~isfieldRecursive(handles.Dmeg{1},'CRC','score')),
        n = 1;
        handles.Dmeg{1}.CRC.score = {};
    else
        n = find(ismember(handles.Dmeg{1}.CRC.score(2,:),'Oracle'),1,'first');
        if(isempty(n))
            n = size(handles.Dmeg{1}.CRC.score,2)+1;
        end
    end
    
    %save relative confidence
    handles.Dmeg{1}.relConfidence = [];
    handles.Dmeg{1}.relConfidence = scores(:,2);
    
    handles.Dmeg{1}.CRC.score{1,n} = scores(:,1)';
    handles.Dmeg{1}.CRC.score{2,n} = 'Oracle';
    handles.Dmeg{1}.CRC.score{3,n} = 30;
    handles.Dmeg{1}.CRC.score{4,n} = [0.005,length(scores(:,1))*30];
    handles.Dmeg{1}.CRC.score{5,n} = [];
    handles.Dmeg{1}.CRC.score{6,n} = [];
    handles.Dmeg{1}.CRC.score{7,n} = [];
    
    %Save artifacts if available
    if isfield(response,'artifact'),
        artifacts = response.artifact;
        
        n_art = find(ismember(handles.Dmeg{1}.CRC.score(2,:),'Oracle_Artifact'),1,'first');
        if(isempty(n_art))
            n_art = size(handles.Dmeg{1}.CRC.score,2)+1;
        end
        
        handles.Dmeg{1}.CRC.score{1,n_art} = artifacts';
        handles.Dmeg{1}.CRC.score{2,n_art} = 'Oracle_Artifact';
        handles.Dmeg{1}.CRC.score{3,n_art} = 5;
        handles.Dmeg{1}.CRC.score{4,n_art} = [0.005,length(artifacts)*5];
        handles.Dmeg{1}.CRC.score{5,n_art} = [];
        handles.Dmeg{1}.CRC.score{6,n_art} = [];
        handles.Dmeg{1}.CRC.score{7,n_art} = [];
        
    end
    
    handles.flags.Dmeg{1} = save(handles.Dmeg{1});
    % Update handles structure
    
    %update the uimenu to add this new scorer
    delete(get(handles.score_user,'Children'));
    handles.score = handles.Dmeg{1}.CRC.score;
    for isc=1:size(handles.score,2)
        handles.scorers{isc} = uimenu(handles.score_user,'Label', ...
            char(handles.score{2,isc}),'Callback',@defined_scorer) ;
        handles.namesc{isc}=char(handles.score{2,isc});
    end
    handles.num_scorers=isc;
    handles.currentscore=isc;
    handles.scorers{isc+1}=uimenu(handles.score_user,'Label', ...
        'New scorer','Callback', @new_scorer,...
        'Separator','on') ;
    set(handles.scorers{handles.currentscore},'Checked','on');
    set(handles.figure1,'CurrentAxes',handles.axes4);
    crc_hypnoplot(handles.axes4, ...
        handles,handles.score{3,handles.currentscore});
    set(handles.figure1,'CurrentAxes',handles.axes1);
    % Setting up the add artefact/arousal/event of interest.
    handles.adddeb = ones(1,size(handles.score,2));
    handles.addardeb = ones(1,size(handles.score,2));
    handles.add_eoi = ones(1,size(handles.score,2));
    mainplot(handles)
    % Update del user menu
    delete(get(handles.del_user,'Children'));
    handles.score = handles.Dmeg{1}.CRC.score;
    for isc=1:size(handles.score,2)
        handles.scorers{isc} = uimenu(handles.del_user,'Label', ...
            char(handles.score{2,isc}),'Callback',@delete_scorer) ;
        handles.namesc{isc}=char(handles.score{2,isc});
    end
    % Update handles structure
    guidata(hObject, handles);
    msgbox('Sleep scores imported successfully','Success');
catch e
    disp(e.message);
    msgbox('Error importing sleep scores.','Error');
    return
end



% --------------------------------------------------------------------
function score_importcsv_Callback(hObject, eventdata, handles)
% hObject    handle to score_import (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path,indx] = uigetfile('*.csv');
if ~indx
    return;
end
try
    
    scores = csvread([path file]);
    prompt = {'Enter scorer name:','Epoch size (secs):'};
    defaultans = {file,'30'};
    input = inputdlg(prompt,'Import Sleep Scores',[1 40], defaultans);
    if isempty(input{1}) || isempty(input{2}),
        return;
    end
    name = input{1};
    winsize = str2double(input{2});
    n_imp = find(ismember(handles.Dmeg{1}.CRC.score(2,:),name),1,'first');
    if(isempty(n_imp))
        n_imp = size(handles.Dmeg{1}.CRC.score,2)+1;
    end
    
    handles.Dmeg{1}.CRC.score{1,n_imp} = scores';
    handles.Dmeg{1}.CRC.score{2,n_imp} = name;
    handles.Dmeg{1}.CRC.score{3,n_imp} = winsize;
    handles.Dmeg{1}.CRC.score{4,n_imp} = [0.005,length(scores)*winsize];
    handles.Dmeg{1}.CRC.score{5,n_imp} = [];
    handles.Dmeg{1}.CRC.score{6,n_imp} = [];
    handles.Dmeg{1}.CRC.score{7,n_imp} = [];
    
    handles.Dmeg{1} = save(handles.Dmeg{1});
    
    %update the uimenu to add this new scorer
    delete(get(handles.score_user,'Children'));
    handles.score = handles.Dmeg{1}.CRC.score;
    for isc=1:size(handles.score,2)
        handles.scorers{isc} = uimenu(handles.score_user,'Label', ...
            char(handles.score{2,isc}),'Callback',@defined_scorer) ;
        handles.namesc{isc}=char(handles.score{2,isc});
    end
    handles.num_scorers=isc;
    handles.currentscore=isc;
    handles.scorers{isc+1}=uimenu(handles.score_user,'Label', ...
        'New scorer','Callback', @new_scorer,...
        'Separator','on') ;
    set(handles.scorers{handles.currentscore},'Checked','on');
    set(handles.figure1,'CurrentAxes',handles.axes4);
    crc_hypnoplot(handles.axes4, ...
        handles,handles.score{3,handles.currentscore});
    set(handles.figure1,'CurrentAxes',handles.axes1);
    % Setting up the add artefact/arousal/event of interest.
    handles.adddeb = ones(1,size(handles.score,2));
    handles.addardeb = ones(1,size(handles.score,2));
    handles.add_eoi = ones(1,size(handles.score,2));
    mainplot(handles)
    % Update del user menu
    delete(get(handles.del_user,'Children'));
    handles.score = handles.Dmeg{1}.CRC.score;
    for isc=1:size(handles.score,2)
        handles.scorers{isc} = uimenu(handles.del_user,'Label', ...
            char(handles.score{2,isc}),'Callback',@delete_scorer) ;
        handles.namesc{isc}=char(handles.score{2,isc});
    end
    % Update handles structure
    guidata(hObject, handles);
    msgbox('Sleep scores imported successfully','Success');
catch
    msgbox('Error importing sleep scores.','Error');
    return
end


% --------------------------------------------------------------------
function score_import_Callback(hObject, eventdata, handles)
% hObject    handle to score_import (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

impfile = spm_select(1, 'any', 'Select mat file with other score','' ...
    ,pwd,'\.[mM][aA][Tt]');

D = crc_eeg_load(impfile);

try
    D.CRC.score;
    okscore = 1;
catch
    okscore = 0;
end

oksamples = D.nsamples==handles.Dmeg{1}.nsamples;
okchannels = length(chanlabels(D))==length(chanlabels(handles.Dmeg{1}));

dlg={};
if ~okscore
    if isempty(dlg);
        dlg = {'No score found in the selected file'};
    else
        dlg = {dlg{:} 'No score found in the selected file'};
    end
end

if ~oksamples
    if isempty(dlg)
        dlg = {'The selected file does not have the same amount of samples than the viewed one'};
    else
        dlg ={ dlg{:} 'The selected file does not have the same amount of samples than the viewed one'};
    end
end

if ~okchannels
    if isempty(dlg)
        dlg = {'The selected file does not have the same amount of channels than the viewed one'};
    else
        dlg = {dlg{:} 'The selected file does not have the same amount of channels than the viewed one'};
    end
end

if oksamples && okchannels && okscore
    %Check that the imported scores have the correct dimensions
    if size(D.CRC.score,2)==7
        for i=1:size(D.CRC.score,2)
            D.CRC.score{8,i}=[];
        end
    end
    % Create new menu without Guide
    
    if handles.figz~=0
        z = handles.figz;
    else
        z = figure;
    end
    
    figure(z)
    handles.figz = z;
    clf(z);
    
    set(z,...% The main GUI figure
        'MenuBar','none', ...
        'Toolbar','none', ...
        'HandleVisibility','callback', ...
        'Color', get(0,...
        'defaultuicontrolbackgroundcolor'));
    set(z,'FileName',impfile);
    set(z,'Position',[218 210 990 632])
    
    hPlotAxes = axes(...    % Axes for plotting the selected plot
        'Parent', z, ...
        'Units', 'normalized', ...
        'HandleVisibility','callback', ...
        'Position',[0.11 0.13 0.80 0.67]);
    
    Mainhandle = gcbo;
    
    sc=[];
    for i=1:size(D.CRC.score,2)
        sc=[sc, {D.CRC.score{2,i}}];
    end
    sc=[sc {'All'}];
    hPlotsPopupmenu = uicontrol(... % List of available types of plot
        'Parent', z, ...
        'Units','normalized',...
        'Position',[0.11 0.85 0.45 0.1],...
        'HandleVisibility','callback', ...
        'String',sc,...
        'Callback', {@plothypnoimport,handles,Mainhandle},...
        'Style','popupmenu');
    
    hUpdateButton = uicontrol(... % Button for updating selected plot
        'Parent', z, ...
        'Units','normalized',...
        'HandleVisibility','callback', ...
        'Position',[0.6 0.85 0.3 0.1],...
        'String','Import score',...
        'Callback', {@saveimportscore,handles,Mainhandle});
    
    crc_hypnoplot(hPlotAxes,handles,D.CRC.score{3,1}, ...
        D.CRC.score{1,1})
    
else
    errordlg(dlg, 'Error');
end

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function saveimportscore(hObject, eventdata, handles, gcbo)

handles = guidata(gcbo);

fig=get(hObject,'Parent');
Child = get(fig, 'Children');

for ii = 1:length(Child)
    if  strcmp(get(Child(ii),'Type'),'axes')
        axesidx = Child(ii);
    elseif strcmp(get(Child(ii),'Type'),'uicontrol')
        if strcmp(get(Child(ii),'Style'),'popupmenu')
            popidx = Child(ii);
        end
    end
end

impfile=get(fig,'FileName');
D=crc_eeg_load(impfile);

popVal=get(popidx,'Value');
if size(D.CRC.score,2)>=popVal
    
    D.CRC.score(:,popVal)
    if size(handles.Dmeg{1}.CRC.score,1) == size(D.CRC.score(:,popVal),1)
        handles.Dmeg{1}.CRC.score=[handles.Dmeg{1}.CRC.score D.CRC.score(:,popVal)];
    else
        errordlg('Scores do not have the same sizes', 'Error');
    end
    
    handles.addardeb=[handles.addardeb 1];
    handles.adddeb = [handles.adddeb 1];
    handles.add_eoi = [handles.add_eoi 1];
else
    for ii=1:size(D.CRC.score,2)
        if size(handles.Dmeg{1}.CRC.score,1) == size(D.CRC.score(:,ii),1)
            handles.Dmeg{1}.CRC.score=[handles.Dmeg{1}.CRC.score D.CRC.score(:,ii)];
        else
            errordlg('Scores do not have the same sizes', 'Error');
        end
        handles.addardeb=[handles.addardeb 1];
        handles.adddeb = [handles.adddeb 1];
        handles.add_eoi = [handles.add_eoi 1];
    end
end

handles.score=handles.Dmeg{1}.CRC.score;
D=handles.Dmeg{1};
save(D);

close(fig)

guidata(gcbo, handles);

% --------------------------------------------------------------------
function plothypnoimport(hObject, eventdata, handles, gcbo)

handles = guidata(gcbo);

fig=get(hObject,'Parent');
Child = get(fig, 'Children');

for ii = 1:length(Child)
    if  strcmp(get(Child(ii),'Type'),'axes')
        axesidx = Child(ii);
    end
end

impfile=get(fig,'FileName');
D=crc_eeg_load(impfile);
Valpop=get(hObject,'Value');

if size(D.CRC.score,2)>=Valpop
    crc_hypnoplot(axesidx, handles,...
        D.CRC.score{3,Valpop},D.CRC.score{1,Valpop},D.CRC.score{2,Valpop})
end

guidata(gcbo, handles);

% --------------------------------------------------------------------
function score_compare_Callback(hObject, eventdata, handles)
% hObject    handle to score_compare (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.figz~=0
    z = handles.figz;
else
    z = figure;
end

if size(handles.Dmeg{1}.CRC.score,2)<2
    disp('Only one score in the file, can not be compared')
    close(z)
    return
end

% figure(z)
% close(z)
figure(z)
handles.figz = z;
clf(z);

set(z,...% The main GUI figure
    'MenuBar','none', ...
    'Toolbar','none', ...
    'HandleVisibility','callback', ...
    'Color', get(0,...
    'defaultuicontrolbackgroundcolor'));
set(z,'Position',[218 210 990 632])

Mainhandle = gcbo;

handles.popcmp(1) = uicontrol(... % List of available types of plot
    'Parent', z, ...
    'Units','normalized',...
    'Position',[0.125 0.485 0.2 0.5],...
    'HandleVisibility','callback', ...
    'String',handles.score(2,:),...
    'Callback', {@update_hypnocomp,handles,Mainhandle},...
    'Style','popupmenu');

handles.popcmp(2) = uicontrol(... % List of available types of plot
    'Parent', z, ...
    'Units','normalized',...
    'Position',[0.675 0.485 0.2 0.5],...
    'HandleVisibility','callback', ...
    'String',handles.score(2,:),...
    'Callback', {@update_hypnocomp,handles,Mainhandle},...
    'Style','popupmenu');
set(handles.popcmp(2),'Value',2);

hUpdateButton = uicontrol(... % Button for updating selected plot
    'Parent', z, ...
    'Units','normalized',...
    'HandleVisibility','callback', ...
    'Position',[0.35 0.025 0.3 0.05],...
    'String','Merge Score',...
    'Callback', {@crc_hypnomerge,handles,Mainhandle});

handles.subcmp(1) = subplot(311);
handles.subcmp(2) = subplot(312);
handles.subcmp(3) = subplot(313);

% Plot subaxes
set(z,'CurrentAxes',handles.subcmp(1))
crc_hypnoplot(handles.subcmp(1),...
    handles,...
    handles.Dmeg{1}.CRC.score{3,1},...
    handles.Dmeg{1}.CRC.score{1,1},handles.Dmeg{1}.CRC.score{2,1})

set(z,'CurrentAxes',handles.subcmp(2))
crc_hypnoplot(handles.subcmp(2),...
    handles,...
    handles.Dmeg{1}.CRC.score{3,2},...
    handles.Dmeg{1}.CRC.score{1,2},handles.Dmeg{1}.CRC.score{2,2})

[mtch,perc_match] = crc_hypnocompare([handles.subcmp(3) z], [1 2], ...
    handles.Dmeg{1}.CRC.score, handles);
handles.match = mtch;
handles.perc_match = perc_match;

% Update handles structure
guidata(hObject, handles);

%% Visual comparison of hypnograms
function update_hypnocomp(hObject, kk1, kk2, gcbo)

handles = guidata(gcbo);

fig = get(hObject,'Parent');

Valpop=get(hObject,'Value');

if handles.popcmp(1)==hObject
    cla(handles.subcmp(1))
    set(fig,'CurrentAxes',handles.subcmp(1))
    crc_hypnoplot(handles.subcmp(1),...
        handles,...
        handles.Dmeg{1}.CRC.score{3,Valpop},...
        handles.Dmeg{1}.CRC.score{1,Valpop},...
        handles.Dmeg{1}.CRC.score{2,Valpop})
    
elseif handles.popcmp(2)==hObject
    cla(handles.subcmp(2))
    set(fig,'CurrentAxes',handles.subcmp(2))
    crc_hypnoplot(handles.subcmp(2),...
        handles,...
        handles.Dmeg{1}.CRC.score{3,Valpop},...
        handles.Dmeg{1}.CRC.score{1,Valpop},...
        handles.Dmeg{1}.CRC.score{2,Valpop})
end

Val1 = get(handles.popcmp(1),'Value');
Val2 = get(handles.popcmp(2),'Value');

[mtch,perc_match] = crc_hypnocompare([handles.subcmp(3) fig], ...
    [Val1 Val2], ...
    handles.Dmeg{1}.CRC.score, handles);
handles.match = mtch;
handles.perc_match = perc_match;

guidata(gcbo, handles);

%% Merging of hypnograms
function crc_hypnomerge(hObject, kk1, handles, gcbo)
% FORMAT crc_mergehypno(hand,pl,handles,windowsize,score,scorer)
% Merging 2 hypnograms into 1.

handles = guidata(gcbo);

fig = get(hObject,'Parent');

notmatchidx = find(handles.match~=0);

Val1 = get(handles.popcmp(1),'Value');
Val2 = get(handles.popcmp(2),'Value');

% create a third column for artefact on single channel if there are only
% column
if size(handles.Dmeg{1}.CRC.score{5,Val1},2)<3
    handles.Dmeg{1}.CRC.score{5,Val1}(:,3) = 0;
end
if size(handles.Dmeg{1}.CRC.score{5,Val2},2)<3
    handles.Dmeg{1}.CRC.score{5,Val2}(:,3) = 0;
end

if handles.Dmeg{1}.CRC.score{3,Val1}==handles.Dmeg{1}.CRC.score{3,Val2}
    
    tmpscore = handles.Dmeg{1}.CRC.score{1,Val1};
    tmpscore(notmatchidx) = NaN; %#ok<*FNDSB> % put a NaN where there are mismatch
    
    handles.Dmeg{1}.CRC.score=...
        [handles.Dmeg{1}.CRC.score handles.Dmeg{1}.CRC.score(:,Val1)];
    
    handles.Dmeg{1}.CRC.score{1,end} = tmpscore;
    handles.Dmeg{1}.CRC.score{2,end} = ...
        ['Merge ' handles.Dmeg{1}.CRC.score{2,Val1} '-', ...
        handles.Dmeg{1}.CRC.score{2,Val2} ];
    handles.Dmeg{1}.CRC.score{4,end}(1) = ...
        min(handles.Dmeg{1}.CRC.score{4,Val1}(1), ...
        handles.Dmeg{1}.CRC.score{4,Val2}(1));
    handles.Dmeg{1}.CRC.score{4,end}(2) = ...
        max(handles.Dmeg{1}.CRC.score{4,Val1}(2), ...
        handles.Dmeg{1}.CRC.score{4,Val2}(2));
    handles.Dmeg{1}.CRC.score{5,end} = ...
        [handles.Dmeg{1}.CRC.score{5,Val1}; ...
        handles.Dmeg{1}.CRC.score{5,Val2}];
    handles.Dmeg{1}.CRC.score{6,end} = ...
        [handles.Dmeg{1}.CRC.score{6,Val1} ; ...
        handles.Dmeg{1}.CRC.score{6,Val2}];
    handles.Dmeg{1}.CRC.score{7,end} = ...
        [handles.Dmeg{1}.CRC.score{7,Val1} ; ...
        handles.Dmeg{1}.CRC.score{7,Val2}];
    
    if isfield(handles.Dmeg{1}.CRC,'sc_merge')
        Nmerge = size(handles.Dmeg{1}.CRC.sc_merge,1);
        handles.Dmeg{1}.CRC.sc_merge{Nmerge+1,1} = handles.Dmeg{1}.CRC.score{2,end};
        handles.Dmeg{1}.CRC.sc_merge{Nmerge+1,2} = handles.perc_match;
    else
        handles.Dmeg{1}.CRC.sc_merge{1,1} = handles.Dmeg{1}.CRC.score{2,end};
        handles.Dmeg{1}.CRC.sc_merge{1,2} = handles.perc_match;
    end
    
    handles.addardeb = [handles.addardeb 1];
    handles.adddeb   = [handles.adddeb 1];
    handles.add_eoi  = [handles.add_eoi 1];
    
    handles.score = handles.Dmeg{1}.CRC.score;
    D = handles.Dmeg{1};
    save(D);
    
    delete(get(handles.score_user,'Children'));
    for isc=1:size(handles.score,2)
        handles.scorers{isc} = uimenu(handles.score_user,'Label', ...
            char(handles.score(2,isc)),'Callback',@defined_scorer) ;
        handles.namesc{isc} = char(handles.score(2,isc));
    end
    handles.num_scorers = isc;
    handles.scorers{isc+1} = uimenu(handles.score_user,'Label', ...
        'New scorer','Callback',@new_scorer,...
        'Separator', 'on') ;
    set(handles.scorers{handles.currentscore},'Checked','on');
    delete(get(handles.axes4,'Children'));
    set(handles.figure1,'CurrentAxes',handles.axes4);
    crc_hypnoplot(handles.axes4, handles, ...
        handles.score{3,handles.currentscore}, ...
        handles.score{1,handles.currentscore}, ...
        handles.score{2,handles.currentscore})
    set(handles.figure1,'CurrentAxes',handles.axes1);
    
    delete(get(handles.del_user,'Children'));
    for isc=1:size(handles.score,2)
        handles.scorers{isc} = uimenu(handles.del_user,'Label', ...
            char(handles.score(2,isc)),'Callback',@delete_scorer) ;
        handles.namesc{isc} = char(handles.score(2,isc));
    end
    close(fig);
else
    errordlg({'The window size of the score are not the same. The scores cannot be merged'}, 'Error');
end

guidata(gcbo, handles);

% --------------------------------------------------------------------
function score_check_Callback(hObject, eventdata, handles)
% hObject    handle to score_import (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get stuff out
score  = handles.score;
scorer = handles.currentscore;
i_win  = str2double(get(handles.currentpage,'string'));
n_disp = crc_get_defaults('score.check_disnum');

% Do the check
[n_empty,l_empty] = crc_check_hypno(score, scorer, i_win, n_disp);

% jumpt to next empty window or not
if crc_get_defaults('score.check_jump') && ~isempty(n_empty)
    set(handles.slider1,'val',(n_empty-1)*handles.winsize)
    mainplot(handles)
end

% % Update handles structure
% guidata(hObject, handles);

%--------------------------------------------------------------------------
% menu on right click
%--------------------------------------------------------------------------
% --------------------------------------------------------------------
function ArtefactMenu_Callback(hObject, eventdata, handles)
% hObject    handle to ArtefactMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Mouse = get(handles.axes1,'CurrentPoint');
guidata(hObject,handles)

% --------------------------------------------------------------------
function Delar_Callback(hObject, eventdata, handles)
% hObject    handle to Delar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function Deleoi_Callback(hObject, eventdata, handles)
% hObject    handle to Deleoi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function Delart_Callback(hObject, eventdata, handles)
% hObject    handle to Delart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Mouse=get(handles.axes1,'CurrentPoint');
minimum=min(min(abs(handles.score{5,handles.currentscore}-Mouse(1,1))));
[row,col]=find((abs(handles.score{5,handles.currentscore}-Mouse(1,1))-minimum)==0);

if length(row)==2
    handles.adddeb(handles.currentscore) = 1;
    if handles.unspecart
        set(handles.addundefart,'Label','Add "start undefined Artefact"');
    else
        set(handles.addspecart,'Label','Add "start specified Artefact"');
    end
end

handles.score{5,handles.currentscore} = ...
    [handles.score{5,handles.currentscore}(1:row-1,:) ; ...
    handles.score{5,handles.currentscore}(row+1:size(handles.score{5,handles.currentscore},1),:)];

handles.score{8,handles.currentscore} = ...
    [handles.score{8,handles.currentscore}(1:row-1) ; ...
    handles.score{8,handles.currentscore}(row+1:size(handles.score{8,handles.currentscore},1))];

set(handles.axes1,'Color',[1 1 1]);
% Save the changes
handles.Dmeg{1}.CRC.score  =handles.score;
D = handles.Dmeg{1};
save(D);

set(handles.figure1,'CurrentAxes',handles.axes4)
crc_hypnoplot(handles.axes4, ...
    handles,handles.score{3,handles.currentscore})
set(handles.figure1,'CurrentAxes',handles.axes1)

mainplot(handles)

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function Delaro_Callback(hObject, eventdata, handles)
% hObject    handle to Delaro (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Mouse=get(handles.axes1,'CurrentPoint');
minimum=min(min(abs(handles.score{6,handles.currentscore}-Mouse(1,1))));
[row,col]=find((abs(handles.score{6,handles.currentscore}-Mouse(1,1))-minimum)==0);

if length(row)==2
    handles.addardeb(handles.currentscore) = 1;
    set(handles.addaro,'Label','Add "start Arousal" point');
end

handles.score{6,handles.currentscore} = ...
    [handles.score{6,handles.currentscore}(1:row-1,:) ;...
    handles.score{6,handles.currentscore}(row+1:size(handles.score{6,handles.currentscore},1),:)];
set(handles.axes1,'Color',[1 1 1]);
%Save the changes
handles.Dmeg{1}.CRC.score=handles.score;
D=handles.Dmeg{1};
save(D);

set(handles.figure1,'CurrentAxes',handles.axes4)
crc_hypnoplot(handles.axes4, ...
    handles,handles.score{3,handles.currentscore})
set(handles.figure1,'CurrentAxes',handles.axes1)

mainplot(handles)

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function Del_eoi_Callback(hObject, eventdata, handles)
% hObject    handle to Del_eoi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Mouse=get(handles.axes1,'CurrentPoint');
minimum=min(min(abs(handles.score{7,handles.currentscore}-Mouse(1,1))));
[row,col]=find((abs(handles.score{7,handles.currentscore}-Mouse(1,1))-minimum)==0);

if length(row)==2
    handles.add_eoi(handles.currentscore) = 1;
    set(handles.addeoi,'Label','Add "start Event of interest" point');
end

handles.score{7,handles.currentscore} = ...
    [handles.score{7,handles.currentscore}(1:row-1,:) ;...
    handles.score{7,handles.currentscore}(row+1:size(handles.score{7,handles.currentscore},1),:)];

%Save the changes
handles.Dmeg{1}.CRC.score=handles.score;
D=handles.Dmeg{1};
save(D);

set(handles.figure1,'CurrentAxes',handles.axes4)
crc_hypnoplot(handles.axes4, ...
    handles,handles.score{3,handles.currentscore})
set(handles.figure1,'CurrentAxes',handles.axes1)

mainplot(handles)

% Update handles structure
guidata(hObject, handles);

function Delartonlyone_Callback(hObject, eventdata, handles)
% hObject    handle to Delartonlyone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Mouse   =   get(handles.axes1,'CurrentPoint');
chan    =   ceil((Mouse(1,2)-handles.scale/2)/handles.scale);
channel =   handles.inddis(chan);
art_concerned   =   find(handles.score{5,handles.currentscore}(:,3) == channel);
minimum =   min(min(abs(handles.score{5,handles.currentscore}(art_concerned,:)-Mouse(1,1))));
[row,col]   =   find((abs(handles.score{5,handles.currentscore}-Mouse(1,1))-minimum)==0);
handles.score{5,handles.currentscore} = ...
    [handles.score{5,handles.currentscore}(1:row-1,:) ; ...
    handles.score{5,handles.currentscore}(row+1:size(handles.score{5,handles.currentscore},1),:)];

handles.score{8,handles.currentscore} = ...
    [handles.score{8,handles.currentscore}(1:row-1) ; ...
    handles.score{8,handles.currentscore}(row+1:size(handles.score{8,handles.currentscore},1))];
[dum1,r,dum] = intersect(handles.chan, channel);
handles.chan = [handles.chan(1:r-1) handles.chan(r+1:length(handles.chan))];
set(handles.axes1,'Color',[1 1 1]);

%Save the changes
handles.Dmeg{1}.CRC.score   =   handles.score;
D   =   handles.Dmeg{1};
save(D);

set(handles.figure1,'CurrentAxes',handles.axes4)
crc_hypnoplot(handles.axes4, ...
    handles,handles.score{3,handles.currentscore})
set(handles.figure1,'CurrentAxes',handles.axes1)

mainplot(handles)
% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function addundefart_Callback(hObject, eventdata, handles)
% hObject    handle to addundefart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.figure1,'CurrentAxes',handles.axes1)
Mouse   =   handles.Mouse;
timedbtart  =   Mouse(1,1);
handles.unspecart   =   1;
lab     =   'unspecified';

% add the eighth line
if size(handles.score,1)<8
    for isc=1:size(handles.score,2)
        handles.score{8,isc}=cell(size(handles.score{5,isc},1),1);
    end
end

% add the zeros on the last column (to adapt with configuration for artifacts on single channel
if size(handles.score{5,handles.currentscore},2)<3
    for ii = 1 : size(handles.score{5,handles.currentscore},1)
        handles.score{5,handles.currentscore}(:,3)  =   0;
    end
end

if handles.adddeb(handles.currentscore) == 1
    if ~isempty(handles.score{5,handles.currentscore})
        if sum(and(and(timedbtart>handles.score{5,handles.currentscore}(:,1),...
                timedbtart<handles.score{5,handles.currentscore}(:,2)),handles.score{5,handles.currentscore}(:,3)==0))
            beep
            disp('Invalid "start Artefact" point')
            disp(' ')
        else
            handles.score{5,handles.currentscore}= ...
                [handles.score{5,handles.currentscore}; timedbtart timedbtart 0];
            handles.score{8,handles.currentscore}=[handles.score{8,handles.currentscore};...
                {lab}];
            set(handles.addundefart,'Label','Add "end undefined Artefact"');
            handles.adddeb(handles.currentscore) = 0;
            fact    =   min(str2double(get(handles.NbreChan,'String')),length(handles.indnomeeg) + length(handles.indexMEEG));
            plot(ones(1,2)*timedbtart,[0 handles.scale*(fact+1)], ...
                'UIContextMenu',handles.Deletemenu,'Color',[0 0 0])
            text(timedbtart,handles.scale*(fact+6/8), ...
                'S.Art','Color',[0 0 0],'FontSize',14)
            set(handles.axes1,'Color',[0.9 0.7 0.7])
            
        end
        
    else
        handles.score{5,handles.currentscore}=...
            [handles.score{5,handles.currentscore}; timedbtart timedbtart 0];
        handles.score{8,handles.currentscore}=[handles.score{8,handles.currentscore};...
            {lab}]; %for each artefact, we save the name of the artefact
        set(handles.addundefart,'Label','Add "end undefined Artefact"');
        handles.adddeb(handles.currentscore) = 0;
        fact    =  min(str2double(get(handles.NbreChan,'String')),length(handles.indnomeeg) + length(handles.indexMEEG));
        plot(ones(1,2)*timedbtart,[0 handles.scale*(fact+1)], ...
            'UIContextMenu',handles.Deletemenu,'Color',[0 0 0])
        
        text(timedbtart,handles.scale*(fact+6/8), ...
            'S.Art','Color',[0 0 0],'FontSize',14)
        set(handles.axes1,'Color',[0.9 0.7 0.7])
    end
else
    Mouse   =    get(handles.axes1,'CurrentPoint');
    timefinart  =   Mouse(1,1);
    [row]   =   find(and(and(timefinart>handles.score{5,handles.currentscore}(:,1),...
        timefinart>handles.score{5,handles.currentscore}(:,2)),handles.score{5,handles.currentscore}(:,3)==0));
    test3   =   sum(and(handles.score{5,handles.currentscore}...
        (size(handles.score{5,handles.currentscore},1),1)<...
        handles.score{5,handles.currentscore}(row,1), ...
        handles.score{5,handles.currentscore}...
        (size(handles.score{5,handles.currentscore},1),2)<...
        handles.score{5,handles.currentscore}(row,2)));
    if or(or(sum(and(and(timefinart>handles.score{5,handles.currentscore}(:,1),...
            timefinart<handles.score{5,handles.currentscore}(:,2)),handles.score{5,handles.currentscore}(:,3)==0))>0, ...
            timefinart<handles.score{5,handles.currentscore}...
            (size(handles.score{5,handles.currentscore},1),2)),test3)
        beep
        disp('Invalid "end Artefact" point')
        disp(' ')
    else
        set(handles.addundefart,'Label','Add "start undefined Artefact"');
        handles.score{5,handles.currentscore}...
            (size(handles.score{5,handles.currentscore},1),2)= timefinart;
        handles.adddeb(handles.currentscore) = 1;
        fact    =   min(str2double(get(handles.NbreChan,'String')),length(handles.indnomeeg) + length(handles.indexMEEG));
        plot(ones(1,2)*timedbtart,[0 handles.scale*(fact+1)], ...
            'UIContextMenu',handles.Deletemenu,'Color',[0 0 0])
        text(timedbtart,handles.scale*(fact+6/8), ...
            'E.Art','Color',[0 0 0],'FontSize',14)
        set(handles.axes1,'Color',[1 1 1])
        set(handles.figure1,'CurrentAxes',handles.axes4)
        plot(handles.score{5,handles.currentscore}(size(handles.score{5,handles.currentscore},1),1),ones(1,1)*8,...
            '+','MarkerSize',8,'MarkerFaceColor',[0.4 0.4 0.4],...
            'MarkerEdgeColor',[0.4 0.4 0.4],'tag','undart')
        set(handles.figure1,'CurrentAxes',handles.axes1)
    end
end
handles.Dmeg{1}.CRC.score   =   handles.score;
D   =   handles.Dmeg{1};
save(D);

% Update handles structure
guidata(hObject, handles);
fftchan_Callback(hObject, [],handles)

%--------------------------------------------------------------------------
%Edit state of one channels when right click--------
%--------------------------------------------------------------------------

% --------------------------------------------------------------------
function addonlyone_Callback(hObject,eventdata,handles)
% hObject    handle to addonlyone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function addstart(hObject,eventdata,handles)

set(handles.figure1,'CurrentAxes',handles.axes1)

a=get(hObject,'Parent');
b=get(a,'Parent');
c=get(b,'Parent');
handles=guidata(c);

% add the eighth line
if size(handles.score,1)<8
    for isc=1:size(handles.score,2)
        handles.score{8,isc}=cell(size(handles.score{5,isc},1),1);
    end
end

% add the zeros on the last column
if size(handles.score{5,handles.currentscore},2)<3
    for ii = 1 : size(handles.score{5,handles.currentscore},1)
        handles.score{5,handles.currentscore}(:,3)=0;
    end
end

%Time pointed
chan = ceil((handles.Mouse(1,2)-handles.scale/2)/handles.scale);
time = handles.Mouse(1,1);
Chanslidval = get(handles.Chanslider,'Value');
slidpos     = Chanslidval-rem(Chanslidval,1);
NbreChandisp    = str2double(get(handles.NbreChan,'String'));
index   = [handles.indnomeeg handles.indexMEEG];
handles.inddis = index(slidpos : 1 : slidpos + NbreChandisp -1);
chantodel   =   get(gcbo,'Label');
channel     =   handles.inddis(chan);
addend      =   intersect(handles.chan, channel);
artchan     =   find(handles.score{5,:}(:,3)==channel);
in_art = find(handles.score{5,:}(artchan,1)<time & handles.score{5,:}(artchan,2)>time);
%label noted in the last line to describe the artefacts
lab = 'unspecified on only one channel';

if ~isempty(strfind(chantodel,'Start'))
    if ~isempty(handles.score{5,handles.currentscore})
        if or(~isempty(addend),~isempty(in_art))
            beep
            disp('Invalid "start Artefact" point, you must first end the last artefact on this channel')
            disp(' ')
        else
            [ch, line, dum]=intersect(handles.score{5,handles.currentscore}(:,3), channel);
            teststart =0;
            for l=1:length(line)
                teststart = sum(or(and(handles.score{5,handles.currentscore}(line(l),2)>time, ...
                    handles.score{5,handles.currentscore}(line(l),1)<time),...
                    handles.score{5,handles.currentscore}(line(l),2)==handles.score{5,handles.currentscore}(line(l),1)));
            end
            if teststart~=0
                beep
                disp('Invalid "start Artefact" point, this zone is already been recovered')
                disp(' ')
            else
                tfinal = nsamples(handles.Dmeg{1})/fsample(handles.Dmeg{1});
                handles.score{5,handles.currentscore}=[handles.score{5,handles.currentscore}; time tfinal channel];
                handles.score{8,handles.currentscore}=[handles.score{8,handles.currentscore}; {lab}];
                handles.chan =[handles.chan channel];
            end
        end
    else
        tfinal = nsamples(handles.Dmeg{1})/fsample(handles.Dmeg{1});
        handles.score{5,handles.currentscore}=[handles.score{5,handles.currentscore}; time tfinal channel];
        handles.score{8,handles.currentscore}=[handles.score{8,handles.currentscore}; {lab}];
        handles.chan =[handles.chan channel];
    end
else
    if isempty(addend)
        beep
        disp('Invalid "end Artefact" point, you must first determine the "Start artefact on this channel"')
        disp(' ')
    else
        [ch,check,dum] = intersect(handles.score{5,handles.currentscore}(:,3),channel);
        if length(check)>1
            [ch, in, dum] = intersect(handles.score{5,handles.currentscore}(:,3), channel);
            for i = 1:length(check)-1
                chk = sum(or(handles.score{5,handles.currentscore}(in(end),1)>timefin),...
                    and(handles.score{5,handles.currentscore}(check(i),1)<timefin,...
                    handles.score{5,handles.currentscore}(check(i),2)<timefin));
            end
            if chk  ~= 0
                beep
                disp('This end is not available')
                disp(' ')
            else
                handles.score{5,handles.currentscore}(in(end),2) = time;
                handles.chan = setdiff(handles.chan, channel);
            end
        else
            [ch, in, dum] = intersect(handles.score{5,handles.currentscore}(:,3), channel);
            if handles.score{5,handles.currentscore}(in(end),1)>time
                beep
                disp('This end is not available')
                disp(' ')
            else
                handles.score{5,handles.currentscore}(in(end),2) = time;
                handles.chan = setdiff(handles.chan, channel);
            end
        end
    end
end
% ---- save data ----
handles.Dmeg{1}.CRC.score = handles.score;
D = handles.Dmeg{1};
save(D);

% ---- update menu ----
delete(get(handles.addonlyone,'Children'));
uimenu(handles.addonlyone,'Label', ...
    ('Start artefact on this channel'),'Callback',{@addstart,handles}) ;
if ~isempty(handles.chan)
    uimenu(handles.addonlyone,'Label', ...
        'End artefact on this channel ','Callback',{@addstart,handles});
end

set(handles.figure1,'CurrentAxes',handles.axes4)
crc_hypnoplot(handles.axes4,handles,handles.score{3,handles.currentscore})
set(handles.figure1,'CurrentAxes',handles.axes1)
% Update handles structure
guidata(c, handles);
fftchan_Callback(hObject,[],handles)
mainplot(handles)
% --------------------------------------------------------------------
function addspecart_Callback(hObject, eventdata, handles)
% hObject    handle to addspecart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.figure1,'CurrentAxes',handles.axes1)

%--------------------------------------------------------------------------
function addtypeart(hObject,eventdata,handles)
% hObject    handle to addspecart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.figure1,'CurrentAxes',handles.axes1)
a   =   get(hObject,'Parent');
b   =   get(a,'Parent');
c   = 	get(b,'Parent');
handles	=   guidata(c);
handles.unspecart   =   0;

if size(handles.score,1)<8
    for isc=1:size(handles.score,2)
        handles.score{8,isc}=cell(size(handles.score{5,isc},1),1);
    end
end

if handles.adddeb(handles.currentscore) == 1
    
    Mouse   =   handles.Mouse;
    timedbtart=Mouse(1,1);
    
    if ~isempty(handles.score{5,handles.currentscore})
        if sum(and(timedbtart>handles.score{5,handles.currentscore}(:,1),...
                timedbtart<handles.score{5,handles.currentscore}(:,2)))
            
            beep
            disp('Invalid "start Artefact" point')
            disp(' ')
            
        else
            lab=get(hObject,'Label');
            
            handles.score{5,handles.currentscore}= ...
                [handles.score{5,handles.currentscore};...
                timedbtart timedbtart 0];
            handles.score{8,handles.currentscore}=[handles.score{8,handles.currentscore};...
                {lab}];
            
            set(handles.addspecart,'Label','Add "end specified Artefact"');
            
            handles.adddeb(handles.currentscore) = 0;
            fact    =   min(str2double(get(handles.NbreChan,'String')),length(handles.indnomeeg) + length(handles.indexMEEG));
            plot(ones(1,2)*timedbtart,[0 handles.scale*(fact+1)], ...
                'UIContextMenu',handles.Deletemenu,'Color',[0 0 0])
            
            text(timedbtart,handles.scale*(fact+6/8), ...
                'S.Art','Color',[0 0 0],'FontSize',14)
            set(handles.axes1,'Color',[0.9 0.7 0.7])
        end
    else
        lab=get(hObject,'Label');
        handles.score{5,handles.currentscore}=...
            [handles.score{5,handles.currentscore};...
            timedbtart timedbtart 0];
        handles.score{8,handles.currentscore}=[handles.score{8,handles.currentscore};...
            {lab}];
        
        set(handles.addspecart,'Label','Add "end specified Artefact"');
        
        handles.adddeb(handles.currentscore) = 0;
        fact    =   min(str2double(get(handles.NbreChan,'String')),length(handles.indnomeeg) + length(handles.indexMEEG));
        plot(ones(1,2)*timedbtart,[0 handles.scale*(fact+1)], ...
            'UIContextMenu',handles.Deletemenu,'Color',[0 0 0])
        
        text(timedbtart,handles.scale*(fact+6/8), ...
            'S.Art','Color',[0 0 0],'FontSize',14)
        set(handles.axes1,'Color',[0.9 0.7 0.7])
        
    end
end
handles.Dmeg{1}.CRC.score=handles.score;
D=handles.Dmeg{1};
save(D);

crcdef = crc_get_defaults('score');
delete(get(handles.addspecart,'Children'));

if strfind(get(handles.addspecart,'Label'),'start')
    for iart=1:size(crcdef.lab_art,2)
        uimenu(handles.addspecart,'Label', ...
            char(crcdef.lab_art{iart}),'Callback',{@addtypeart,handles}) ;
    end
elseif strfind(get(handles.addspecart,'Label'),'end')
    uimenu(handles.addspecart,'Label', ...
        'End artefact','Callback',{@endtypeart,handles}) ;
end

% Update handles structure
guidata(c, handles);

%--------------------------------------------------------------------------
function endtypeart(hObject,eventdata,handles)

a  	=   get(hObject,'Parent');
b   =   get(a,'Parent');
c   =   get(b,'Parent');
handles	=   guidata(c);
Mouse   =   handles.Mouse;
timefinart  =   Mouse(1,1);
handles.unspecart   =   0;
[row]   =   find(and(timefinart>handles.score{5,handles.currentscore}(:,1),...
    timefinart>handles.score{5,handles.currentscore}(:,2)));

test3   =  sum(and(handles.score{5,handles.currentscore}...
    (size(handles.score{5,handles.currentscore},1),1)<...
    handles.score{5,handles.currentscore}(row,1), ...
    handles.score{5,handles.currentscore}...
    (size(handles.score{5,handles.currentscore},1),2)<...
    handles.score{5,handles.currentscore}(row,2)));

if or(or(sum(and(timefinart>handles.score{5,handles.currentscore}(:,1),...
        timefinart<handles.score{5,handles.currentscore}(:,2)))>0, ...
        timefinart<handles.score{5,handles.currentscore}...
        (size(handles.score{5,handles.currentscore},1),2)),test3)
    
    beep
    disp('Invalid "end Artefact" point')
    disp(' ')
else
    set(handles.addspecart,'Label','Add "start specified Artefact"');
    
    handles.score{5,handles.currentscore}...
        (size(handles.score{5,handles.currentscore},1),2)= timefinart;
    
    handles.adddeb(handles.currentscore) = 1;
    fact    =   min(str2double(get(handles.NbreChan,'String')),length(handles.indnomeeg) + length(handles.indexMEEG));
    plot(ones(1,2)*timefinart,[0 handles.scale*(fact+1)], ...
        'UIContextMenu',handles.Deletemenu,'Color',[0 0 0])
    text(timefinart,handles.scale*(fact+6/8), ...
        'E.Art','Color',[0 0 0],'FontSize',14)
    set(handles.axes1,'Color',[1 1 1])
    
end
handles.Dmeg{1}.CRC.score=handles.score;
D=handles.Dmeg{1};
save(D);

crcdef = crc_get_defaults('score');
delete(get(handles.addspecart,'Children'));

if strfind(get(handles.addspecart,'Label'),'start')
    for iart=1:size(crcdef.lab_art,2)
        uimenu(handles.addspecart,'Label', ...
            char(crcdef.lab_art{iart}),'Callback',{@addtypeart,handles}) ;
    end
elseif strfind(get(handles.addspecart,'Label'),'end')
    uimenu(handles.addspecart,'Label', ...
        'End artefact','Callback',{@endtypeart,handles}) ;
end
% Update handles structure
guidata(c, handles);

% --------------------------------------------------------------------
function addaro_Callback(hObject, eventdata, handles)
% hObject    handle to addaro (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Mouse   =   handles.Mouse;
if handles.addardeb(handles.currentscore) == 1
    timedbtaro  =   Mouse(1,1);
    if ~isempty(handles.score{6,handles.currentscore})
        
        if sum(and(timedbtaro>handles.score{6,handles.currentscore}(:,1),...
                timedbtaro<handles.score{6,handles.currentscore}(:,2)))
            
            beep
            disp('Invalid "start Arousal" point')
            disp(' ')
            
        else
            
            handles.score{6,handles.currentscore}=...
                [handles.score{6,handles.currentscore}; ...
                timedbtaro timedbtaro];
            
            set(handles.addaro,'Label','Add "end Arousal" point');
            handles.addardeb(handles.currentscore) = 0;
            fact    =   min(str2double(get(handles.NbreChan,'String')),length(handles.indnomeeg) + length(handles.indexMEEG));
            plot(ones(1,2)*timedbtaro,[0 handles.scale*(fact+1)], ...
                'UIContextMenu',handles.Delar,'Color',[0 0 0])
            text(timedbtaro,handles.scale*(fact+6/8), ...
                'S.Aro.','Color',[1 0 0],'FontSize',14)
            set(handles.axes1,'Color',[0.8 0.8 0.95]);
            
        end
        
    else
        
        handles.score{6,handles.currentscore}=...
            [handles.score{6,handles.currentscore};...
            timedbtaro timedbtaro];
        set(handles.addaro,'Label','Add "end Arousal" point');
        handles.addardeb(handles.currentscore) = 0;
        fact    =   min(str2double(get(handles.NbreChan,'String')),length(handles.indnomeeg) + length(handles.indexMEEG));
        plot(ones(1,2)*timedbtaro,[0 handles.scale*(fact+1)], ...
            'UIContextMenu',handles.Delar,'Color',[0 0 0])
        text(timedbtaro,handles.scale*(fact+6/8), ...
            'S.Aro','Color',[1 0 0],'FontSize',14)
        set(handles.axes1,'Color',[0.8 0.8 0.95]);
        
    end
else
    timefinaro=Mouse(1,1);
    [row]=find(and(timefinaro>handles.score{6,handles.currentscore}(:,1),...
        timefinaro>handles.score{6,handles.currentscore}(:,2)));
    
    test3=sum(and(handles.score{6,handles.currentscore}...
        (size(handles.score{6,handles.currentscore},1),1)<...
        handles.score{6,handles.currentscore}(row,1), ...
        handles.score{6,handles.currentscore}...
        (size(handles.score{6,handles.currentscore},1),2)<...
        handles.score{6,handles.currentscore}(row,2)));
    
    if or(or(sum(and(timefinaro>handles.score{6,handles.currentscore}(:,1),...
            timefinaro<handles.score{6,handles.currentscore}(:,2)))>0, ...
            timefinaro<handles.score{6,handles.currentscore}...
            (size(handles.score{6,handles.currentscore},1),2)),test3)
        
        beep
        disp('Invalid "end Arousal" point')
        disp(' ')
    else
        set(handles.addaro,'Label','Add "start Arousal" point');
        
        handles.score{6,handles.currentscore}...
            (size(handles.score{6,handles.currentscore},1),2)= timefinaro;
        
        handles.addardeb(handles.currentscore) = 1;
        fact    =   min(str2double(get(handles.NbreChan,'String')),length(handles.indnomeeg) + length(handles.indexMEEG));
        plot(ones(1,2)*timefinaro,[0 handles.scale*(fact+1)], ...
            'UIContextMenu',handles.Delar,'Color',[0 0 0])
        text(timefinaro,handles.scale*(fact+6/8), ...
            'E.Aro','Color',[1 0 0],'FontSize',14)
        set(handles.axes1,'Color',[1 1 1]);
        set(handles.figure1,'CurrentAxes',handles.axes4)
        plot(handles.score{6,handles.currentscore}(size(handles.score{6,handles.currentscore},1),1),ones(1,1)*8,...
            '+','MarkerSize',8,'MarkerFaceColor',[1 0 0],...
            'MarkerEdgeColor',[1 0 0],'tag','arou')
        set(handles.figure1,'CurrentAxes',handles.axes1)
    end
end

% Save the changes
handles.Dmeg{1}.CRC.score = handles.score;
save(handles.Dmeg{1});

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function addeoi_Callback(hObject, eventdata, handles)
% hObject    handle to addeoi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Mouse   =   handles.Mouse;
if handles.add_eoi(handles.currentscore) == 1
    timedbteoi=Mouse(1,1);
    if ~isempty(handles.score{7,handles.currentscore})
        if sum(and(timedbteoi>handles.score{7,handles.currentscore}(:,1),...
                timedbteoi<handles.score{7,handles.currentscore}(:,2)))
            
            beep
            disp('Invalid "start Event of interest" point')
            disp(' ')
            
        else
            handles.score{7,handles.currentscore}=...
                [handles.score{7,handles.currentscore}; ...
                timedbteoi timedbteoi];
            
            set(handles.addeoi,'Label','Add "end Event of interest" point');
            handles.add_eoi(handles.currentscore) = 0;
            fact    =   min(str2double(get(handles.NbreChan,'String')),length(handles.indnomeeg) + length(handles.indexMEEG));
            plot(ones(1,2)*timedbteoi,[0 handles.scale*(fact+1)], ...
                'UIContextMenu',handles.Deleoi,'Color',[0 0 0])
            text(timedbteoi,handles.scale*(fact+6/8), ...
                'S.EOI.','Color',[0.75 0.2 0.2],'FontSize',14)
        end
    else
        handles.score{7,handles.currentscore}=...
            [handles.score{7,handles.currentscore};...
            timedbteoi timedbteoi];
        
        set(handles.addeoi,'Label','Add "end Event of interest" point');
        handles.add_eoi(handles.currentscore) = 0;
        fact    =   min(str2double(get(handles.NbreChan,'String')),length(handles.indnomeeg) + length(handles.indexMEEG));
        plot(ones(1,2)*timedbteoi,[0 handles.scale*(fact+1)], ...
            'UIContextMenu',handles.Deleoi,'Color',[0 0 0])
        text(timedbteoi,handles.scale*(fact+6/8), ...
            'S.EOI','Color',[0.75 0.2 0.2],'FontSize',14)
    end
else
    timefineoi=Mouse(1,1);
    
    [row]=find(and(timefineoi>handles.score{7,handles.currentscore}(:,1),...
        timefineoi>handles.score{7,handles.currentscore}(:,2)));
    
    test3=sum(and(handles.score{7,handles.currentscore}...
        (size(handles.score{7,handles.currentscore},1),1)<...
        handles.score{7,handles.currentscore}(row,1), ...
        handles.score{7,handles.currentscore}...
        (size(handles.score{7,handles.currentscore},1),2)<...
        handles.score{7,handles.currentscore}(row,2)));
    
    if or(or(sum(and(timefineoi>handles.score{7,handles.currentscore}(:,1),...
            timefineoi<handles.score{7,handles.currentscore}(:,2)))>0, ...
            timefineoi<handles.score{7,handles.currentscore}...
            (size(handles.score{7,handles.currentscore},1),2)),test3)
        
        beep
        disp('Invalid "end Event of interest" point')
        disp(' ')
    else
        set(handles.addeoi,'Label','Add "start Event of interest" point');
        
        handles.score{7,handles.currentscore}...
            (size(handles.score{7,handles.currentscore},1),2)= timefineoi;
        
        handles.add_eoi(handles.currentscore) = 1;
        fact    =   min(str2double(get(handles.NbreChan,'String')),length(handles.indnomeeg) + length(handles.indexMEEG));
        plot(ones(1,2)*timefineoi,[0 handles.scale*(fact+1)], ...
            'UIContextMenu',handles.Deleoi,'Color',[0 0 0])
        text(timefineoi,handles.scale*(fact+6/8), ...
            'E.EOI','Color',[0.75 0.2 0.2],'FontSize',14)
        set(handles.figure1,'CurrentAxes',handles.axes4)
        plot(handles.score{7,handles.currentscore}(size(handles.score{7,handles.currentscore},1),1),ones(1,1)*8,...
            '+','MarkerSize',8,'MarkerFaceColor',[0.75 0.2 0.2],...
            'MarkerEdgeColor',[0.75 0.2 0.2],'tag','eoi')
        set(handles.figure1,'CurrentAxes',handles.axes1)
        
    end
end

% Save the changes
handles.Dmeg{1}.CRC.score = handles.score;
save(handles.Dmeg{1});

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function FPL_Callback(hObject, eventdata, handles)
% hObject    handle to FPL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Mouse	=  	handles.Mouse;
cpoint  =   Mouse(1,1);
handles.score{4,handles.currentscore}(1)=cpoint;

% Display everything on the axes
mainplot(handles);

set(handles.figure1,'CurrentAxes',handles.axes4)

crc_hypnoplot(handles.axes4,...
    handles,handles.score{3,handles.currentscore});

set(handles.figure1,'CurrentAxes',handles.axes1);

% Save the changes
handles.Dmeg{1}.CRC.score = handles.score;
save(handles.Dmeg{1});

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function OPL_Callback(hObject, eventdata, handles)
% hObject    handle to OPL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Mouse   =   handles.Mouse;
cpoint  =   Mouse(1,1);

handles.score{4,handles.currentscore}(2)=cpoint;

mainplot(handles)

set(handles.figure1,'CurrentAxes',handles.axes4)

crc_hypnoplot(handles.axes4,...
    handles,handles.score{3,handles.currentscore})

set(handles.figure1,'CurrentAxes',handles.axes1);

% Save the changes
handles.Dmeg{1}.CRC.score = handles.score;
save(handles.Dmeg{1});

% Update handles structure
guidata(hObject, handles);

%To be checked
% --------------------------------------------------------------------
function newtype_Callback(hObject, eventdata, handles)

% hObject    handle to spike (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

a = gcbo;
b = get(a,'Parent');
handles = guidata(gcbo);

guidata(b, handles);

set(handles.figure1, 'windowbuttonmotionfcn', '')

%create new type structure
prompt = {'Please enter a name for this type of spike'};
def = {'Newtype'};
num_lines = 1;
dlg_title = 'Name of a new kind of spike';
type = inputdlg(prompt,dlg_title,num_lines,def);

%update the uimenu to add this new TYPE
ntype = size(handles.type,1);
handles.type(ntype+1,1) = type;
handles.type(ntype+1,2) = cellstr('red');

delete(get(handles.manevent,'Children'));

for i = 1:size(handles.type,1)
    uimenu(handles.manevent,'Label', char(handles.type(i,1)),'Callback',@Define_event)
end

uimenu(handles.manevent,'Label', ...
    'New Type','Callback', @newtype_Callback,...
    'Separator','on') ;

%update the popupmenu10 (events)
current_events = get(handles.popupmenu10,'String');
events_updated = [current_events; type];
set(handles.popupmenu10,'String',events_updated,'Value',length(events_updated))

%update the popupmenu11 (Value of the event = scorer)
current_value = get(handles.popupmenu11,'String');
scorer  = char(handles.score{2,handles.currentscore});
check   = intersect(current_value,scorer);
if isempty(check)
    current_value{end+1} = scorer;
end
set(handles.popupmenu11,'String',current_value,'Value',length(current_value))

set(handles.figure1, 'windowbuttonmotionfcn', @update_powerspect)

guidata(b, handles);
% Update handles structure
mainplot(handles)

function Define_event(hObject, eventdata)

b = get(gcbo,'Parent');
c = get(b,'Parent');
d = get(c,'Parent');

handles = guidata(gcbo);

%Parameters of the new event
type    =   get(gcbo,'Label');
Mouse	=  	handles.Mouse;
cpoint  =   Mouse(1,1);

%Update the events
D   =   handles.Dmeg{1};
ev  =   events(D);
Nev =   size(ev,2);

ev(Nev+1).type  = char(type);
ev(Nev+1).value = char(handles.score{2,handles.currentscore});
ev(Nev+1).time  = cpoint;
ev(Nev+1).duration  = 0;
ev(Nev+1).offset    = 0;

%save
D = events(D,1,ev);
D.CRC.goodevents    =   ones(1,numel(ev));
save(D);
handles.Dmeg{1} = D;

%redefine handles.evt :
if ~isempty(ev)
    for i   =   1   :   size(ev,2)
        if ~isempty(ev(i)) && isnumeric(ev(i).value)
            ev(i).value    =   num2str(ev(i).value);
        end
    end
end
handles.evt = ev;

%add the scorer in value for events
current_value = get(handles.popupmenu11,'String');
scorer  = char(handles.score{2,handles.currentscore});
check   = intersect(current_value,scorer);
if isempty(check)
    current_value{end+1} = scorer;
end
set(handles.popupmenu11,'String',current_value,'Value',1);

%redefine handles.chosevt
evnum  = get(handles.popupmenu10,'Value');
evtype  = get(handles.popupmenu10,'String');
todisp = evtype(evnum);

if evnum==1
    evnum = 2 : size(evtype,1);
end

chos=[];

for i=1:size(handles.evt,2)
    if any(strcmpi(handles.evt(i).type,evtype(evnum)))
        chos=[chos, i];
    end
end

handles.chosevt = chos;

%new types
evtype = [{'All'} unique({ev(:).type})];
[valtype intevt intdis] = intersect(evtype, todisp);
if isempty(intevt)
    intevt = 1;
elseif intevt>length(evtype)
    intevt = 1;
end
set(handles.popupmenu10,'String',evtype,'Value',intevt);

%new values
pmstring = [{'All'},{'Scan number'},unique({ev(:).value})];

for i=1:size(chos,2)
    if ~strcmpi(handles.evt(chos(i)).value, pmstring)
        pmstring = [pmstring, {handles.evt(chos(i)).value}];
    end
end

handles.displevt    = 0;

set(handles.popupmenu11,...
    'String',pmstring,...
    'Value',length(pmstring))

%Add an landmark on the top axes
set(handles.figure1,'CurrentAxes',handles.axes4);
crc_hypnoplot(handles.axes4, ...
    handles,handles.score{3,handles.currentscore});
set(handles.figure1,'CurrentAxes',handles.axes1);

%Store data in the figure's application
guidata(d, handles);
%plot the new event
mainplot(handles)

function Delete_Event_Callback(hObject, eventdata, handles)

Mouse = get(handles.axes1,'CurrentPoint');
locx  = Mouse(1,1);

%search the event to delete
ev = events(handles.Dmeg{1});
locev = cell2mat({(ev(:).time)});
loc = locev - locx;

[m indm]= min(abs(loc(:)));

Stock1 = ev(1 : indm-1);
Stock2 = ev(indm+1 : size(ev,2));

ev = [Stock1 Stock2];
%save the modifications
D = handles.Dmeg{1};
D = events(D,1,ev);
D.CRC.goodevents    =   ones(1,numel(ev));
save(D);
handles.Dmeg{1} = D;

%redefine handles.evt :
if ~isempty(ev)
    for i   =   1   :   size(ev,2)
        if ~isempty(ev(i)) && isnumeric(ev(i).value)
            ev(i).value    =   num2str(ev(i).value);
        end
    end
end
handles.evt = ev;

%redefine handles.chosevt
evnum  = get(handles.popupmenu10,'Value');
evtype  = get(handles.popupmenu10,'String');
todisp = evtype(evnum);

if evnum==1
    evnum = 2 : size(evtype,1);
end

chos=[];

for i=1:size(handles.evt,2)
    if any(strcmpi(handles.evt(i).type,evtype(evnum)))
        chos=[chos, i];
    end
end

handles.chosevt = chos;

%new types
evtype = [{'All'} unique({ev(:).type})];
[valtype intevt intdis] = intersect(evtype, todisp);
if isempty(intevt)
    intevt = 1;
elseif intevt>length(evtype)
    intevt = 1;
end
set(handles.popupmenu10,'String',evtype,'Value',intevt);

%new values
pmstring = [{'All'},{'Scan number'},unique({ev(:).value})];

for i=1:size(chos,2)
    if ~strcmpi(handles.evt(chos(i)).value, pmstring)
        pmstring = [pmstring, {handles.evt(chos(i)).value}];
    end
end

handles.displevt    = 0;

set(handles.popupmenu11,...
    'String',pmstring,...
    'Value',length(pmstring))

set(handles.figure1,'CurrentAxes',handles.axes4);
crc_hypnoplot(handles.axes4, ...
    handles,handles.score{3,handles.currentscore});
set(handles.figure1,'CurrentAxes',handles.axes1);

%update and save
guidata(hObject,handles)
mainplot(handles)

function manevent_Callback(hObject,eventdata,handles)
% hObject    handle to DetectSpike (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function Deletemenu_Callback(hObject, eventdata, handles)
% hObject    handle to Deletemenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --------------------------------------------------------------------

%To be checked%%%%%%%%%%%%%%
% --------------------------------------------------------------------
function DeletEvents_Callback(hObject, eventdata, handles)
% hObject    handle to DeletEvents (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --------------------------------------------------------------------
handles.Mouse = get(handles.axes1,'CurrentPoint');
guidata(hObject,handles)


%To be checked%%%%%%%%%%%%%%

% --------------------------------------------------------------------
function Deletemenuone_Callback(hObject, eventdata, handles)
% hObject    handle to Deletemenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --------------------------------------------------------------------

function grid_Callback(hObject, eventdata, handles)
% hObject    handle to grid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function vertgrid_Callback(hObject, eventdata, handles)
% hObject    handle to vertgrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmpi(get(gcbo,'Checked'),'on')
    set(gcbo,'Checked','off')
    handles.vert_grid=0;
else
    set(gcbo,'Checked','on')
    handles.vert_grid=1;
end
mainplot(handles)

guidata(gcbo, handles);

% --------------------------------------------------------------------
function horgrid_Callback(hObject, eventdata, handles)
% hObject    handle to horgrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmpi(get(gcbo,'Checked'),'on')
    set(gcbo,'Checked','off')
    handles.hor_grid=0;
else
    set(gcbo,'Checked','on')
    handles.hor_grid=1;
end

mainplot(handles)

guidata(gcbo, handles);

% --------------------------------------------------------------------
function score_user_Callback(hObject, eventdata, handles)
% hObject    handle to score_user (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function new_scorer(hObject, eventdata)
% hObject    handle to score_user (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
%get handles of the main parent
a = get(gcbo,'Parent');
b = get(a,'Parent');
c = get(b,'Parent');
%to be checked = old version contained : get(c); What does it means?
handles=guidata(gcbo);

%create new scoring structure
prompt = {'Please enter your name'};
def = {'Newuser'};
num_lines = 1;
dlg_title = 'Name of the new scorer';
handles.score(2,handles.num_scorers +1) = inputdlg(prompt,dlg_title,num_lines,def);

%Choosing size of window to score
prompt = {'Please choose the size of the scoring windows'};
def = {num2str(crc_get_defaults('score.winsize'))};
num_lines = 1;
dlg_title = 'Size of the scoring windows (in sec)';
handles.score{3,handles.num_scorers +1} = inputdlg(prompt,dlg_title,num_lines,def);
handles.score{3,handles.num_scorers +1} = str2double(handles.score{3,handles.num_scorers +1});

% Creating FPL & OPL
handles.score{4,handles.num_scorers +1} = [1/fsample(handles.Dmeg{1}) ...
    nsamples(handles.Dmeg{1})/fsample(handles.Dmeg{1})];

% Creating artefacts
handles.score{5,handles.num_scorers +1} = [];

% Creating arousals
handles.score{6,handles.num_scorers +1} = [];

% Creating event of interest
handles.score{7,handles.num_scorers +1} = [];

% Putting NaN in the score
handles.score{1,handles.num_scorers +1} = ...
    0/0*ones(1,ceil(nsamples(handles.Dmeg{1}) / ...
    (handles.score{3,handles.num_scorers +1}*fsample(handles.Dmeg{1}))));

D = handles.Dmeg{1};
if isfield(D, 'CRC')
    D.CRC.score=handles.score;
else
    D.CRC = struct('score',[]);
    D.CRC.score=handles.score;
end
handles.num_scorers = handles.num_scorers + 1;
handles.Dmeg{1} = D;
%update the uimenu to add this new scorer
delete(get(handles.score_user,'Children'));
for isc=1:size(handles.score,2)
    handles.scorers{isc} = uimenu(handles.score_user,'Label', ...
        char(handles.score{2,isc}),'Callback',@defined_scorer) ;
    handles.namesc{isc}=char(handles.score{2,isc});
end

delete(get(handles.del_user,'Children'));
for isc=1:size(handles.score,2)
    handles.scorers{isc} = uimenu(handles.del_user,'Label', ...
        char(handles.score{2,isc}),'Callback',@delete_scorer) ;
    handles.namesc{isc}=char(handles.score{2,isc});
end

handles.num_scorers=isc;
handles.currentscore=isc;
handles.scorers{isc+1}=uimenu(handles.score_user,'Label', ...
    'New scorer','Callback', @new_scorer,...
    'Separator','on') ;
set(handles.scorers{handles.currentscore},'Checked','on');
set(handles.figure1,'CurrentAxes',handles.axes4);
crc_hypnoplot(handles.axes4, ...
    handles,handles.score{3,handles.currentscore});
set(handles.figure1,'CurrentAxes',handles.axes1);
% Setting up the add artefact/arousal/event of interest.
handles.adddeb = ones(1,size(handles.score,2));
handles.addardeb = ones(1,size(handles.score,2));
handles.add_eoi = ones(1,size(handles.score,2));
mainplot(handles)
% Update handles structure
guidata(c, handles);
return


function delete_scorer(hObject, eventdata)
% hObject    handle to score_user (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

a = get(gcbo,'Parent');
b = get(a,'Parent');
c = get(b,'Parent');
%To be checked = old version contained : get(c);

handles = guidata(gcbo);
namuser = get(gcbo,'Label');

handles.score = handles.Dmeg{1}.CRC.score;
for isc=1:size(handles.score,2)
    handles.namesc{isc}=char(handles.score{2,isc});
end

for isc=1:size(handles.namesc,2)
    if strcmpi(namuser,handles.namesc{isc})
        delete_id=isc;
    end
end

handles.score = handles.Dmeg{1}.CRC.score;
if size(handles.score,2) <= 1,
    msgbox('You cannot delete the last scorer.','Success');
    return;
end

%keyboard
answer = questdlg(sprintf('Are you sure you want to permanently delete the scorer: %s?',namuser) , ...
    'Confirm delete', ...
    'Yes','No', 'No');
% Handle response
switch answer
    case 'No'
        return
end

handles.Dmeg{1}.CRC.score(:,delete_id) = [];
handles.Dmeg{1} = save(handles.Dmeg{1});

handles.score = handles.Dmeg{1}.CRC.score;
%update the uimenu to add this new scorer
delete(get(handles.score_user,'Children'));

for isc=1:size(handles.score,2)
    handles.scorers{isc} = uimenu(handles.score_user,'Label', ...
        char(handles.score{2,isc}),'Callback',@defined_scorer) ;
    handles.namesc{isc}=char(handles.score{2,isc});
end
handles.num_scorers=isc;
handles.currentscore=isc;
handles.scorers{isc+1}=uimenu(handles.score_user,'Label', ...
    'New scorer','Callback', @new_scorer,...
    'Separator','on') ;
set(handles.scorers{handles.currentscore},'Checked','on');
set(handles.figure1,'CurrentAxes',handles.axes4);
crc_hypnoplot(handles.axes4, ...
    handles,handles.score{3,handles.currentscore});
set(handles.figure1,'CurrentAxes',handles.axes1);
% Setting up the add artefact/arousal/event of interest.
handles.adddeb = ones(1,size(handles.score,2));
handles.addardeb = ones(1,size(handles.score,2));
handles.add_eoi = ones(1,size(handles.score,2));
mainplot(handles)
% Update del user menu
delete(get(handles.del_user,'Children'));
handles.score = handles.Dmeg{1}.CRC.score;
for isc=1:size(handles.score,2)
    handles.scorers{isc} = uimenu(handles.del_user,'Label', ...
        char(handles.score{2,isc}),'Callback',@delete_scorer) ;
    handles.namesc{isc}=char(handles.score{2,isc});
end

% Update handles structure
guidata(c, handles);
return


function defined_scorer(hObject, eventdata)
% hObject    handle to score_user (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

a = get(gcbo,'Parent');
b = get(a,'Parent');
c = get(b,'Parent');
%To be checked = old version contained : get(c);

handles = guidata(gcbo);
namuser = get(gcbo,'Label');

for isc=1:size(handles.namesc,2)
    if strcmpi(namuser,handles.namesc{isc})
        handles.currentscore=isc;
        set(gcbo, 'Checked', 'on');
    else
        set(handles.scorers{isc},'Checked','off')
    end
end

set(handles.figure1,'CurrentAxes',handles.axes4);
crc_hypnoplot(handles.axes4, ...
    handles,handles.score{3,handles.currentscore});
set(handles.figure1,'CurrentAxes',handles.axes1);
mainplot(handles)
% Update handles structure
guidata(c, handles);
guidata(gcbo, handles);
return

% -------------------------------------------------------------------
%To be checked : new subfunctions in click
% --- Click
function click(hObject, eventdata, handles)

Mouse = get(handles.axes4,'CurrentPoint');

if Mouse(1,2) > 0   %Even without score sleep, possibility to access from axes4 where we want on main screen
    slidval = Mouse(1);
    slidval = floor(slidval/handles.winsize)*handles.winsize;
    slidval = max(slidval,1/fsample(handles.Dmeg{1}));
    if slidval==nsamples(handles.Dmeg{1})/fsample(handles.Dmeg{1})
        slidval=str2double(get(handles.currenttime,'String'));
    else
        slidval = min(slidval, ...
            nsamples(handles.Dmeg{1})/fsample(handles.Dmeg{1})-handles.winsize/2);
    end
    set(handles.slider1,'Value',slidval)
    set(handles.figure1,'CurrentAxes',handles.axes4)
    a=get(handles.axes4,'Children');
    for i=1:size(a,1)
        if strcmpi(get(a(i),'tag'),'cursor')
            delete(a(i))
        end
    end
    hold on
    handles.cursor = plot(slidval+handles.winsize/2,8,'v','Color',[0.2 0.2 0.2],'LineWidth',2.5,'tag','cursor');
    handles.displevt=0;
    guidata(hObject, handles);
    mainplot(handles)
end

% To plot the sensors on the localizer
Mouse2 = get(handles.axes1,'CurrentPoint');
X = get(handles.axes1,'XTick');
Y = get(handles.axes1,'YTick');

if ~get(handles.counterspect,'Value') && and(and(Mouse2(1,1)<X(end),Mouse2(1,1)>X(1)),and(Mouse2(1,2)<Y(end),Mouse2(1,2)>Y(1)))
    if isfield(handles,'namesc')
        C = Mouse2(1,1);
        E = [handles.evt(:).time];
        [ve handles.move]=min(abs(E-C));
        if and(any(abs(C-E) < 1),any(strcmpi(handles.evt(handles.move).value,handles.namesc)))
            fprintf(['You are going to move the event ' char(handles.evt(handles.move).type) ' of ' char(handles.evt(handles.move).value) '\n'])
        else
            handles.move = 0;
        end
    else
        handles.move = 0;
    end
    guidata(hObject,handles);
    
    set(handles.figure1,'CurrentAxes',handles.axes5)
    delete(findobj('tag', 'powerspctrm'))   %Effacer le denier power spectrum
    delete(findobj('tag', 'error'))
    chan = get(handles.fftchan,'Value');
    list = cellstr(get(handles.fftchan,'String'));
    listcompl = chanlabels(handles.Dmeg{1},handles.index);
    [dumb1,dumb2,index]	=   intersect(upper(list),upper(handles.names));
    chantodis  = list(chan);
    [dumb1,dumb2,indextodis] = intersect(upper(chantodis),upper(handles.names));
    
    if isempty(indextodis)
        
        id = indchannel(handles.Dmeg{1},chantodis);
        iddis = indchannel(handles.Dmeg{1},list);
        idtot = indchannel(handles.Dmeg{1},listcompl);
        
        xytot = coor2D(handles.Dmeg{1},idtot);
        xbl = xytot(1,:);
        ybl = xytot(2,:);
        
        xytodis = coor2D(handles.Dmeg{1},id);
        xblack = xytodis(1);
        yblack = xytodis(2);
        
        xy = coor2D(handles.Dmeg{1},iddis);
        xblue = xy(1,:);
        yblue = xy(2,:);
        
        cleargraph(handles.axes5)
        
        if (xblack)==0
            cleargraph(handles.figure1,'axes5')
        end
        
        hold on
        plot(xblack,yblack,'kv','tag','localizer'), plot(xblue,yblue,'b+','tag','localizer'),plot(xbl,ybl,'b.','tag','localizer')
        hold off
        
    else
        idxred	=   index(find(handles.crc_types(index)<-1));
        idxblue	=   index(find(handles.crc_types(index)>-2));
        
        xred  	=   handles.pos(1,idxred);
        yred    =   handles.pos(2,idxred);
        
        xblu    =   handles.pos(1,idxblue);
        yblu    =   handles.pos(2,idxblue);
        
        xblack  =   handles.pos(1,indextodis);
        yblack  =   handles.pos(2,indextodis);
        
        cleargraph(handles.axes5)
        
        hold on
        plot(xred,yred,'r+','tag','localizer'), plot(xblu,yblu,'b+','tag','localizer'), plot(xblack,yblack,'kv','tag','localizer')
        hold off
        if and(length(xblu)==0,length(xred)==0)
            cleargraph(handles.figure1,'axes5')
        end
    end
    xlim([0 1])
    ylim([0 1])
    
    set(handles.fftchan,'visible','off');
end

set(handles.figure1,'CurrentAxes',handles.axes5);

% Plot the power spectrum on a bigger figure

%Mouse3 = abs(get(handles.axes5,'CurrentPoint'));
%X = abs(get(handles.axes5,'XTick'));
%Y = abs(get(handles.axes5,'YTick'));
%grid off
%if ~or(isempty(X),isempty(Y))
%    if and(and(Mouse3(1,1)<X(end),Mouse3(1,1)>X(1)),or(and(Mouse3(1,2)<Y(end),Mouse3(1,2)>Y(1)),and(Mouse3(1,2)<Y(1),Mouse3(1,2)>Y(end))))
%        set(handles.figure1,'CurrentAxes',handles.axes5);
%        pl = get(findobj('tag', 'powerspctrm'));
%        if handles.figz~=0
%            z   =   handles.figz;
%        else
%            z   =   figure;
%        end
%        figure(z)
%        axs     =   get(handles.figz,'Children');
%        cleargraph(z)
%            set(handles.figure1,'YData',pl.Ydata)
%            handles.notdetect = 0;
%            handles.figz = z;
%    end
%end


% --- E
function keypress(hObject, eventdata, handles)

crcdef  =   crc_get_defaults('score');
key     =   str2double(get(handles.figure1,'CurrentCharacter'));
currentwindow   =	floor(str2double(get(handles.currenttime,'String')) ...
    /handles.winsize)+1;
%winsize_select  =   str2double(get(handles.edit1,'String'));
%winsize_score   =   handles.Dmeg{1}.CRC.score{3,handles.currentscore};

if handles.scoring && (key>-1 && key<crcdef.nrStage)
    handles.move = 0;
    % Update Score
    maxwindow   =   ceil(nsamples(handles.Dmeg{1})/fsample(handles.Dmeg{1})/...
        handles.score{3,handles.currentscore});
    if currentwindow>maxwindow
        currentwindow    =   maxwindow;
    end
    handles.score{1,handles.currentscore}(currentwindow)    =   key;
    
    % Save Score
    handles.Dmeg{1}.CRC.score = handles.score;
    save(handles.Dmeg{1});
    
    hypjustplot = 0;
    
    set(handles.figure1,'CurrentAxes',handles.axes1);
    % Goes to next window
    slidval     =   str2double(get(handles.currenttime,'String'));
    slidval     =   (floor(slidval/handles.winsize)+1)*(handles.winsize);
    slidval     =   max(slidval,1/fsample(handles.Dmeg{1}));
    
    if slidval>=nsamples(handles.Dmeg{1})/fsample(handles.Dmeg{1})-handles.winsize/2
        slidval     =   currentwindow*handles.winsize;
    else
        slidval     =	min(slidval, currentwindow*handles.winsize);
        %         nsamples(handles.Dmeg{1})/fsample(handles.Dmeg{1})-handles.winsize/2);
    end
    set(handles.slider1,'Value',slidval);
    set(handles.figure1,'CurrentAxes',handles.axes4);
    crc_hypnoplot(handles.axes4,...
        handles,handles.score{3,handles.currentscore});
    %Goes to the previous window
elseif or(get(handles.figure1,'CurrentCharacter')=='B',get(handles.figure1,'CurrentCharacter')=='b')
    handles.move = 0;
    slidval = str2double(get(handles.currenttime,'String'));
    slidval = (floor(slidval/handles.winsize)-1)*handles.winsize;
    slidval = max(slidval,1/fsample(handles.Dmeg{1}));
    slidval = min(slidval, currentwindow*handles.winsize);
    %         nsamples(handles.Dmeg{1})/fsample(handles.Dmeg{1})-handles.winsize/2);
    
    set(handles.slider1,'Value',slidval);
    set(handles.figure1,'CurrentAxes',handles.axes4);
    crc_hypnoplot(handles.axes4,...
        handles,handles.score{3,handles.currentscore});
    hypjustplot=0;
    % Goes to next window
elseif or(get(handles.figure1,'CurrentCharacter')=='F',get(handles.figure1,'CurrentCharacter')=='f')
    handles.move = 0;
    slidval     =   str2double(get(handles.currenttime,'String'));
    slidval     =   (floor(slidval/handles.winsize)+1)*handles.winsize;
    slidval     =   max(slidval,1/fsample(handles.Dmeg{1}));
    
    if slidval>=nsamples(handles.Dmeg{1})/fsample(handles.Dmeg{1})
        slidval     =   currentwindow*handles.winsize;
    else
        slidval     =   min(slidval, currentwindow*handles.winsize);
        %         nsamples(handles.Dmeg{1})/fsample(handles.Dmeg{1})-handles.winsize/2);
    end
    set(handles.slider1,'Value',slidval);
    set(handles.figure1,'CurrentAxes',handles.axes4);
    crc_hypnoplot(handles.axes4,...
        handles,handles.score{3,handles.currentscore});
    hypjustplot=0;
else
    touch = int2str(get(handles.figure1,'CurrentCharacter'));
    if and(handles.move~=0,strcmpi(touch,'28'))
        ev = events(handles.Dmeg{1});
        ev(handles.move).time = ev(handles.move).time - 1/fsample(handles.Dmeg{1}); % On avance par échantillon
        handles.evt = ev;
        handles.Dmeg{1} = events(handles.Dmeg{1},1,ev);
        save(handles.Dmeg{1});
        fprintf(['The event ' char(handles.evt(handles.move).type) ' of ' char(handles.evt(handles.move).value) ' is moving on the left : ' num2str(handles.evt(handles.move).time) 'sec \n'])
        mainplot(handles)
    elseif and(handles.move~=0,strcmpi(touch,'29'))
        ev = events(handles.Dmeg{1});
        ev(handles.move).time + 1/fsample(handles.Dmeg{1});
        ev(handles.move).time = ev(handles.move).time + 1/fsample(handles.Dmeg{1}); % On avance par échantillon
        handles.evt = ev;
        handles.Dmeg{1} = events(handles.Dmeg{1},1,ev);
        save(handles.Dmeg{1});
        fprintf(['The event ' char(handles.evt(handles.move).type) ' of ' char(handles.evt(handles.move).value) ' is moving on the right : ' num2str(handles.evt(handles.move).time) 'sec \n'])
        mainplot(handles)
    else
        beep;
        fprintf('Wrong key pressed! Please choose press a number between 0 and %d.\n',crcdef.nrStage-1)
        for ii = 1    :   crcdef.nrStage
            fprintf('%d : %s\n',ii-1,crcdef.stnames_L{ii})
        end
        fprintf('\n')
    end
    hypjustplot     =  1;
end

if hypjustplot
    slidval = get(handles.slider1,'Value');
    set(handles.figure1,'CurrentAxes',handles.axes4);
    a   =   get(handles.axes4,'Children');
    for i   =   1   :   size(a,1)
        if strcmpi(get(a(i),'Type'),'line')
            delete(a(i));
        end
    end
    hold on
    handles.cursor  =    plot(slidval+handles.winsize/2,8,'v','Color',[0.6 0.6 0.6],'LineWidth',2.5);
    hold off
else
    set(handles.figure1,'CurrentAxes',handles.axes1);
    mainplot(handles)
end
set(handles.figure1,'CurrentAxes',handles.axes1);
% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------
%Options for multiple files comparison-------------------------------------
%--------------------------------------------------------------------------
% --------------------------------------------------------------------
function multfil_Callback(hObject, eventdata, handles)
% hObject    handle to multfil (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function multnames_Callback(hObject, eventdata, handles)
% hObject    handle to multnames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(gcbo, 'Checked'),'on')
    set(gcbo, 'Checked', 'off');
else
    set(gcbo, 'Checked', 'on');
end
mainplot(handles)

% --------------------------------------------------------------------
function multchan_Callback(hObject, eventdata, handles)
% hObject    handle to multchan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function Chantodisp(hObject, eventdata)
% hObject    handle to multchan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
namchan=get(hObject,'Label');
a=get(hObject,'Parent');
b=get(a,'Parent');
c=get(b,'Parent');
handles=guidata(c);
for ii=1:size(handles.chanset,2)
    if strcmpi(namchan,handles.chanset{ii})
        set(handles.multchanlab{ii}, 'Checked', 'on');
        handles.Chantodis=ii;
    else
        set(handles.multchanlab{ii},'Checked','off')
    end
end
mainplot(handles)
% Update handles structure
guidata(c, handles);

% --------------------------------------------------------------------
function multother_Callback(hObject, ~, handles)
% hObject    handle to multother (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1)
try
    close(handles.figz)
end
flags=struct('multcomp',[]);
flags.multcomp=1;
dis_selchan(flags)

% --------------------------------------------------------------------
function multclose_Callback(hObject, eventdata, handles)
% hObject    handle to multclose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.figure1)
try
    close(handles.figz)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%------------------------------ MAIN PLOT FUNCTION ------------------------
%--------------------------------------------------------------------------

function mainplot(handles)

if handles.displevt
    handles.evt(handles.displevt).time;
    slidval=handles.evt(handles.displevt).time-handles.winsize/2;
    if slidval<=1
        slidval=max(handles.evt(handles.displevt).time-2,1/fsample(handles.Dmeg{1})); %default value if first event is at the boundary of the file
    elseif slidval>=(nsamples(handles.Dmeg{1})/fsample(handles.Dmeg{1})-handles.winsize/2)
        slidval=min(handles.evt(handles.displevt).time-2,nsamples(handles.Dmeg{1})/fsample(handles.Dmeg{1})-handles.winsize/2); %default value if last event is at the boundary of the file
    end
else
    slidval = floor((get(handles.slider1,'Value')/handles.winsize))*handles.winsize;
    if slidval<1/fsample(handles.Dmeg{1})
        slidval=1/fsample(handles.Dmeg{1});
        set(handles.slider1,'Value',slidval);
    elseif slidval>nsamples(handles.Dmeg{1})/fsample(handles.Dmeg{1})-handles.winsize/2;
        slidval=nsamples(handles.Dmeg{1})/fsample(handles.Dmeg{1})-handles.winsize/2;
        aval=get(handles.slider1,'Max');
        set(handles.slider1,'Value',aval);
    end
end
set(handles.figure1,'CurrentAxes',handles.axes4)
a=get(handles.axes4,'Children');
for i=1:size(a,1)
    if strcmpi(get(a(i),'tag'),'cursor')
        delete(a(i))
    end
end
hold on
handles.cursor = plot(slidval+handles.winsize/2,8,'v','Color',[0.2 0.2 0.2],'LineWidth',2.5,'tag','cursor');
set(handles.figure1,'CurrentAxes',handles.axes1)

cmap = hsv(length(handles.Dmeg));

if handles.export
    axs=get(handles.currentfigure,'CurrentAxes');
    h=axs;
else
    cleargraph(handles.figure1,'axes1');
    h=handles.axes1;
end
axes(h);

if handles.multcomp
    maxi=length(handles.Dmeg); %number of files
    maxj=1;                    %one channel/file
else
    maxi=1;                    %one file
end
set(handles.currenttime,'String',num2str(round(slidval+handles.winsize/2)));
set(handles.edit1,'String',num2str(handles.winsize));
set(handles.currentpage,'String',num2str(ceil((slidval+handles.winsize/2)/handles.winsize)));
for i=1:maxi
    nosign=1;
    if handles.multcomp
        %Write the names of the files
        if strcmpi(get(handles.multnames,'Checked'),'on')
            [dumb name] = fileparts(handles.Dmeg{i}.fname);
            under       =   find(name=='_');
            name(under) =   ' ';
            text(slidval+4*handles.winsize/10,handles.scale*(i+1/4),name, ...
                'Color',[0.8 0.8 0.8],'FontSize',14)
        end
        %handles timing problem
        start   =   datevec(handles.date(i,1)-handles.mindate);
        start   =   start(4)*60^2+start(5)*60+start(6);
        ending  =   datevec(handles.date(i,2)-handles.mindate);
        ending  =   ending(4)*60^2 + ending(5)*60 + ending(6);
        if or(slidval<start,slidval>ending)
            nosign=0;
        else
            beg     =   slidval - start;
            tdeb    =   round(beg*fsample(handles.Dmeg{i}));
            temps   =   tdeb:1:min(tdeb+(fsample(handles.Dmeg{i})*handles.winsize)-1, ...
                nsamples(handles.Dmeg{i}));
            toshow  =   temps;
            temps   =   temps/fsample(handles.Dmeg{i})+start;
        end
        %handles channel to display
        Ctodis  =   handles.Chantodis;
        [dumb1,dumb2,index]     =   intersect(handles.chanset{Ctodis}, ...
            upper(chanlabels(handles.Dmeg{i})));
    else
        %timing
        tdeb = min(round(slidval*fsample(handles.Dmeg{i})),nsamples(handles.Dmeg{i})-10);
        %channels
        Chanslidval = get(handles.Chanslider,'Value');
        slidpos     = Chanslidval-rem(Chanslidval,1);
        NbreChandisp    = str2double(get(handles.NbreChan,'String'));
        index   = [handles.indnomeeg handles.indexMEEG];
        handles.inddis = index(slidpos : 1 : slidpos + NbreChandisp -1);
        index   =   handles.inddis;
        set(handles.fftchan,'String',chanlabels(handles.Dmeg{1},handles.inddis))
        maxj    =   length(handles.inddis); %number of channels
        temps   =   tdeb:1:min(tdeb+(fsample(handles.Dmeg{i})*handles.winsize)-1, ...
            nsamples(handles.Dmeg{i}));
        toshow  =   temps;
        temps   =   temps/fsample(handles.Dmeg{i});
        win =  floor((str2double(get(handles.currenttime,'String')) - handles.winsize/2)/handles.winsize + 1);
    end
    
    if ~nosign  %if multiple comparison and no signal in the considered file
        text(slidval+1*handles.winsize/50,handles.scale*i, ...
            'No Signal here','Color',[0 0 0],'FontSize',14);
    else
        for j=1:maxj
            hold on
            [dumb1,dumb2,index2] = ...
                intersect(upper(chanlabels(handles.Dmeg{i},index(j))),handles.names);
            %To be checked = old version contained : testcond = chanlabels(handles.Dmeg{i},index(j));
            if handles.multcomp
                factscale = handles.scale*i; %cycle on files if mult comp
            else
                factscale = handles.scale*j; %cycle on channels if one file
            end
            if  any(index(j)==emgchannels(handles.Dmeg{i}))  %strfind(testcond{:},'EMG')
                contents = get(handles.EMGpopmenu,'String');
                selchan=upper(contents{get(handles.EMGpopmenu,'Value')});
                if isfield(handles,'emgscale') && ~isempty(handles.emgscale)
                    scal=handles.emgscale;
                else
                    unemg=units(handles.Dmeg{i},index(j));
                    try
                        scal=eval(unemg{1});
                    catch
                        if strcmpi(unemg{1},'V')
                            scal=10^-6;
                        else
                            scal=1;
                        end
                    end
                end
                filtparam=handles.filter.coeffEMG;
            elseif any(index(j)==eogchannels(handles.Dmeg{i})) %strfind(testcond{:},'EOG')
                contents = get(handles.EOGpopmenu,'String');
                selchan=upper(contents{get(handles.EOGpopmenu,'Value')});
                if isfield(handles,'eogscale') && ~isempty(handles.eogscale)
                    scal=handles.eogscale;
                else
                    unemg=units(handles.Dmeg{i},index(j));
                    try
                        scal=eval(unemg{1});
                    catch
                        if strcmpi(unemg{1},'V')
                            scal=10^-6;
                        else
                            scal=1;
                        end
                    end
                end
                filtparam=handles.filter.coeffEOG;
            elseif any(index(j)==ecgchannels(handles.Dmeg{i})) %strfind(testcond{:},'ECG'
                contents = get(handles.otherpopmenu,'String');
                selchan=upper(contents{get(handles.otherpopmenu,'Value')});
                if isfield(handles,'ecgscale') && ~isempty(handles.ecgscale)
                    scal=handles.ecgscale;
                else
                    unemg=units(handles.Dmeg{i},index(j));
                    try
                        scal=eval(unemg{1});
                    catch
                        if strcmpi(unemg{1},'V')
                            scal=3*10^-5;
                        else
                            scal=1;
                        end
                    end
                end
                filtparam=handles.filter.coeffother;
            elseif any(index(j)==meegchannels(handles.Dmeg{i})) %'EEG', 'MEGMAG', MEGPLANAR'
                contents = get(handles.EEGpopmenu,'String');
                selchan=upper(contents{get(handles.EEGpopmenu,'Value')});
                chtyp=chantype(handles.Dmeg{i},index(j));
                if strcmpi(chtyp,'EEG')
                    if isfield(handles,'eegscale') && ~isempty(handles.eegscale)
                        scal=handles.eegscale;
                    else
                        unemg=units(handles.Dmeg{i},index(j));
                        try
                            scal=eval(unemg{1});
                        catch
                            if strcmpi(unemg{1},'V')
                                scal=10^-6;
                            else
                                scal=1;
                            end
                        end
                    end
                elseif strcmpi(chtyp,'MEGMAG')
                    selchan = 'MEG';
                    if isfield(handles,'megmagscale') && ~isempty(handles.megmagscale)
                        scal=handles.megmagscale;
                    else
                        unemg=units(handles.Dmeg{i},index(j));
                        try
                            scal=eval(unemg{1});
                        catch
                            if strcmpi(unemg{1},'T')
                                scal=5*10^-14;
                            else
                                scal=1;
                            end
                        end
                    end
                elseif strcmpi(chtyp,'MEGPLANAR')
                    selchan = 'MEG';
                    if isfield(handles,'megplanarscale') && ~isempty(handles.megplanarscale)
                        scal=handles.megplanarscale;
                    else
                        unemg=units(handles.Dmeg{i},index(j));
                        try
                            scal=eval(unemg{1});
                        catch
                            if strcmpi(unemg{1},'T/m')
                                scal=5*10^-13;
                            else
                                scal=1;
                            end
                        end
                    end
                elseif strcmpi(chtyp,'LFP')
                    if isfield(handles,'lfpscale') && ~isempty(handles.lfpscale)
                        scal=handles.lfpscale;
                    else
                        unemg=units(handles.Dmeg{i},index(j));
                        try
                            scal=eval(unemg{1});
                        catch
                            if strcmpi(unemg{1},'V')
                                scal=10^-5;
                            else
                                scal=10 ;
                            end
                        end
                    end
                end
                filtparam=handles.filter.coeffother;
            else
                contents = get(handles.EEGpopmenu,'String');
                selchan=upper(contents{get(handles.EEGpopmenu,'Value')});
                chtyp=chantype(handles.Dmeg{i},index(j));
                if strcmpi(chtyp,'other') || strcmpi(chtyp,'unknown')
                    if isfield(handles,'otherscale') && ~isempty(handles.otherscale)
                        scal=handles.otherscale;
                    else
                        unemg=units(handles.Dmeg{i},index(j));
                        try
                            scal=eval(unemg{1});
                        catch
                            if strcmpi(unemg{1},'V')
                                scal=10^-6;
                            else
                                scal=1;
                            end
                            
                        end
                    end
                end
                filtparam=handles.filter.coeffother;
            end
            %Plot data
            switch  selchan
                case 'MEG'
                    % MEG data
                    normalize   =   get(handles.normalize,'Value');
                    chtyp = chantype(handles.Dmeg{i},index(j));
                    if  strcmpi(chtyp,'MEGMAG')
                        plt     =   plot(temps,factscale+(handles.Dmeg{i}(index(j),toshow))/scal,'Color',[0.25 0.25 0.25]);
                    elseif  strcmpi(chtyp,'MEGPLANAR')
                        if (~normalize)
                            plt     =   plot(temps,factscale+(handles.Dmeg{i}(index(j),toshow))/scal,'Color',[0.5 0.5 0]);
                        else
                            plt     =   plot(temps,factscale+(((handles.Dmeg{i}(index(j),toshow)).^2+...
                                (handles.Dmeg{i}(index(j)+1,toshow)).^2).^0.5)/scal,'Color',[0 0.5 0]);
                        end
                    end
                case    'REF1'
                    plt         = 	plot(temps,factscale + (handles.Dmeg{i}(index(j),toshow))/scal);
                case    'MEAN OF REF'
                    basedata    =   handles.Dmeg{i}(index(j),toshow);
                    ref2idx     =   find(strcmp(chanlabels(handles.Dmeg{i}),'REF2'));
                    scnddata    =   handles.Dmeg{i}(index(j),toshow) - handles.Dmeg{i}(ref2idx,toshow);
                    toplotdat   =   mean([basedata ; scnddata]);
                    plt         =   plot(temps,factscale+(toplotdat)/scal);
                    
                case    'M1-M2'
                    basedata    =   handles.Dmeg{i}(index(j),toshow);
                    M1idx       =   find(strcmp(chanlabels(handles.Dmeg{i}),'M1'));
                    M2idx       =   find(strcmp(chanlabels(handles.Dmeg{i}),'M2'));
                    meanM       =   mean([handles.Dmeg{i}(M1idx,toshow) ; ...
                        handles.Dmeg{i}(M2idx,toshow)]);
                    toplotdat   = 	basedata - meanM;
                    plt         =   plot(temps,factscale+(toplotdat)/scal);
                    
                case    'BIPOLAR'
                    if handles.crc_types(index2)>0
                        [dumb1,index3] = ...
                            intersect(upper(chanlabels(handles.Dmeg{i})), ...
                            upper(handles.names(handles.crc_types(index2))));
                    else
                        index3 = [];
                    end
                    if ~isempty(index3)
                        bipolar	=   handles.Dmeg{i}(index(j),toshow) - handles.Dmeg{i}(index3,toshow);
                        plt  	=   plot(temps,factscale+(bipolar)/scal,'Color','k');
                    else
                        plt     =   plot(temps,factscale+(handles.Dmeg{i}(index(j),toshow))/scal,'Color','k');
                    end
                    
                otherwise
                    [dumb1,index3]  =   intersect(upper(chanlabels(handles.Dmeg{i})),selchan);
                    basedata        =   handles.Dmeg{i}(index(j),toshow);
                    toplotdat       =   basedata - handles.Dmeg{i}(index3,toshow);
                    plt  	=   plot(temps,factscale+(toplotdat)/scal);
            end
            filterlowhigh(plt,i,handles,filtparam,factscale)
            if plt~=0 && any(index(j)==ecgchannels(handles.Dmeg{i}))
                set(plt,'Color',[1 0.1 0.1])
            end
            if plt~=0 && handles.multcomp
                set(plt,'Color',cmap(i,:))
            end
            
            %artefacts on single channel  (To be checked)
            if  ~isempty(handles.artefact)
                deb_epoch   =   find(and(handles.artefact>temps(1),handles.artefact<=temps(end)+1/fsample(handles.Dmeg{1})));%1 seconde = fenetre d'analyse pour la détection des artefacts de courtes durées
                for epo     =   1 : length(deb_epoch)
                    X   =   get(plt,'YData');
                    deb =  (handles.artefact(deb_epoch(epo)) - (win-1)*handles.winsize/1 - 1)*fsample(handles.Dmeg{1}) + 1; %1=temps d'une époque artefactée
                    fin =  min(length(X)+win*handles.winsize*fsample(handles.Dmeg{1}),(handles.artefact(deb_epoch(epo)) - (win-1)*handles.winsize/1)*fsample(handles.Dmeg{1}));
                    time_epo = [(handles.artefact(deb_epoch(epo))-1) + 1/fsample(handles.Dmeg{1}) : 1/fsample(handles.Dmeg{1}) : ((handles.artefact(deb_epoch(epo))-1) + 1/fsample(handles.Dmeg{1}))+length(deb:fin)/fsample(handles.Dmeg{1})-1/fsample(handles.Dmeg{1})];
                    X   =   X(deb:fin);%fs est le nombre d'échantillons contenus dans une seconde
                    plot(time_epo,X,'-.','color','r');
                end
            end
            %artefacts channels: automatic
            if ~isempty(handles.badchannels)&& ~isempty(find(handles.badchannels == win)) %isfield(handles.Dmeg{1}.CRC,'badchannels')
                [electr gum]=find(handles.badchannels == win);
                electr = handles.chanlab(electr);
                if any(index(j) ==  electr)
                    X   =   get(plt,'YData');
                    plot(temps,X,'-r');
                end
            end
            %artefacts on single channel  (To be checked)
            if handles.scoring && ~isempty(handles.score{5,handles.currentscore}) && any(index(j)== handles.score{5,handles.currentscore}(:,3))
                [indice gum]= find(handles.score{5,handles.currentscore}(:,3)==index(j));
                l=1;
                tdebs   =   str2double(get(handles.currenttime,'String'))-handles.winsize/2;
                tfins   =	str2double(get(handles.currenttime,'String'))+handles.winsize/2;
                while l<=length(indice)
                    if(or(or(handles.score{5,handles.currentscore}(indice(l),2)>tdebs,handles.score{5,handles.currentscore}(indice(l),1)<tdebs),...
                            or(handles.score{5,handles.currentscore}(indice(l),1)<tfins,handles.score{5,handles.currentscore}(indice(l),2)>tfins)))
                        [dum ind]=intersect(handles.inddis,index(j));
                        if handles.score{5,handles.currentscore}(indice(l),1)>tdebs
                            text(handles.score{5,handles.currentscore}(indice(l),1),...
                                (ind*handles.scale + handles.scale/2), ...
                                'Start Bad Chan','Color',[0 0 0],'FontSize',12) ;
                            plot(ones(1,2)*handles.score{5,handles.currentscore}(indice(l),1), ...
                                [(ind*handles.scale-handles.scale/2) (ind*handles.scale + handles.scale/2)], ...
                                'Color',[0 0 0])
                        end
                        if handles.score{5,handles.currentscore}(indice(l),2)<tfins
                            text(handles.score{5,handles.currentscore}(indice(l),2),...
                                (ind*handles.scale + handles.scale/2), ...
                                'End Bad Chan','Color',[0 0 0],'FontSize',12)
                            plot(ones(1,2)*handles.score{5,handles.currentscore}(indice(l),2), ...
                                [(ind*handles.scale-handles.scale/2) (ind*handles.scale + handles.scale/2)], ...
                                'Color',[0 0 0])
                        end
                        if handles.score{5,handles.currentscore}(indice(l),2) ~= handles.score{5,handles.currentscore}(indice(l),1)
                            sfin = min(tfins*fsample(handles.Dmeg{1}),handles.score{5,handles.currentscore}(indice(l),2)*fsample(handles.Dmeg{1}));
                        else
                            sfin = tfins*fsample(handles.Dmeg{1});
                        end
                        %tdeb = min(round(slidval*fsample(handles.Dmeg{1})),nsamples(handles.Dmeg{1})-10);
                        sdebut = max(handles.score{5,handles.currentscore}(indice(l),1)*fsample(handles.Dmeg{1}), tdebs*fsample(handles.Dmeg{1}));
                        time = sdebut+1 : sfin;
                        time = time/fsample(handles.Dmeg{1});
                        X = get(plt,'YData');
                        
                        %Plot
                        X = X(sdebut-tdebs*fsample(handles.Dmeg{1})+1:sfin-tdebs*fsample(handles.Dmeg{1}));
                        if handles.export
                            plot(time,X,'Color',[0.8 0.8 0.8]);
                        else
                            plot(time,X,'UIContextMenu',handles.Deletemenuone,'Color',[0.8 0.8 0.8]);
                        end
                    end
                    l=l+1;
                end
            end
            plt = 0;
        end
    end
end

if handles.displevt
    timevt=handles.evt(handles.displevt).time;
    szevt=handles.scale*size(index,2);
    hold on
    plot(timevt*ones(1,szevt+1),handles.scale/2:handles.scale/2+szevt,'-r','Linewidth',2)
end

%display the labels on the y-axis
if handles.multcomp
    li=length(handles.Dmeg);
else
    li=length(index);
end

ylim([0 handles.scale*(li+1)])
set(handles.axes1,'YTick',[handles.scale/2 :handles.scale/2:li*handles.scale+handles.scale/2]);
ylabels=[num2str(round(handles.scale/2))];

for j = 1 : li
    if handles.multcomp
        ylabels =[ylabels {num2str(j)}];
    else
        if and(get(handles.normalize,'Value'),strcmpi(chantype(handles.Dmeg{1},index(j)),'MEGPLANAR'))
            stringn = char(chanlabels(handles.Dmeg{1},index(j)+2));
            string = [stringn(1:end-1), 'N'];
            ylabels  = [ylabels {string}];
        else
            ylabels = [ylabels chanlabels(handles.Dmeg{1},index(j))];
        end
    end
    ylabels = [ylabels num2str(round(handles.scale/2))];
end
set(handles.axes1,'YTickLabel',ylabels);

xlim([temps(1) temps(1)+handles.winsize])
xtick = get(handles.axes1,'XTick');
if isfield(handles,'offset')
    xtick = mod(xtick + handles.offset,24*60^2);
end
[time string] = crc_time_converts(xtick);
set(handles.axes1,'XTickLabel',string)

% display horizontal grid
if handles.hor_grid
    for i = 1:li
        plot([temps(1) temps(end)],[(35+handles.scale*i) (35+handles.scale*i)], ...
            ':','Color',[0.6 0.6 0.6])
        plot([temps(1) temps(end)],[(handles.scale*i-35) (handles.scale*i-35)], ...
            ':','Color',[0.6 0.6 0.6])
    end
end

% display vertical grid
if handles.vert_grid
    Ax=get(handles.axes1,'XTick');
    for ii = Ax(1)-1:1:Ax(end)+1
        plot([ii ii],[0 (handles.scale*(1+li))],':','Color',[0.6 0.6 0.6])
    end
end

%Update the score of the current window
if handles.scoring && handles.winsize == handles.score{3,handles.currentscore} % New gardian to check the window size corresponds to those is used to score file
    ll  =   str2double(get(handles.NbreChan,'String'));
    fact    =   min(ll,length(handles.indexMEEG) + length(handles.indnomeeg));
    currentwindow = floor(str2double(get(handles.currenttime,'String'))/handles.winsize)+1;
    if currentwindow>size(handles.score{1,handles.currentscore},2)
        currentwindow = currentwindow-1;
    end
    curscore = handles.score{1,handles.currentscore}(currentwindow);
    text(temps(end-round(0.8*fsample(handles.Dmeg{1}))),handles.scale*(fact+5/8),num2str(curscore),'Color','k','FontSize',14)
    switch  curscore
        case    0
            plot(temps,((li+1)*handles.scale-1/10000000)*ones(1,length(temps)), ...
                'linewidth',10,'color',[0.2 0.75 0.6])
        case    1
            plot(temps,((li+1)*handles.scale-1/10000000)*ones(1,length(temps)), ...
                'linewidth',10,'color',[0 0.8 1])
        case    2
            plot(temps,((li+1)*handles.scale-1/10000000)*ones(1,length(temps)), ...
                'linewidth',10,'color',[0.1 0.5 0.9])
        case    3
            plot(temps,((li+1)*handles.scale-1/10000000)*ones(1,length(temps)), ...
                'linewidth',10,'color',[0.1 0.2 0.8])
        case    4
            plot(temps,((li+1)*handles.scale-1/10000000)*ones(1,length(temps)), ...
                'linewidth',10,'color',[0.1 0.15 0.5])
        case    5
            plot(temps,((li+1)*handles.scale-1/10000000)*ones(1,length(temps)), ...
                'linewidth',10,'color',[0.5 0.5 0.9])
        case    6
            plot(temps,((li+1)*handles.scale-1/10000000)*ones(1,length(temps)), ...
                'linewidth',10,'color',[0.9 0.4 0.4])
        case    7
            plot(temps,((li+1)*handles.scale-1/10000000)*ones(1,length(temps)), ...
                'linewidth',10,'color',[0.9 0.6 0.3])
    end
    
    % Plot opl & fpl
    tdebs = str2double(get(handles.currenttime,'String'))-handles.winsize/2;
    tfins = str2double(get(handles.currenttime,'String'))+handles.winsize/2;
    fpl = find(and(handles.score{4,handles.currentscore}(:,1)>tdebs,...
        handles.score{4,handles.currentscore}(:,1)<tfins));
    opl = find(and(handles.score{4,handles.currentscore}(:,2)>tdebs,...
        handles.score{4,handles.currentscore}(:,2)<tfins));
    for i=1:length(fpl)
        plot(ones(1,2)*handles.score{4,handles.currentscore}(fpl(i),1),[0 handles.scale*(fact+1)],'Color',[0 0 0])
        text(handles.score{4,handles.currentscore}(fpl(i),1),handles.scale*(fact+6/8),'FPL','Color',[0 0 0.9],'FontSize',14)
    end
    for i=1:length(opl)
        plot(ones(1,2)*handles.score{4,handles.currentscore}(opl(i),2),[0 handles.scale*(fact+1)],'Color',[0 0 0])
        text(handles.score{4,handles.currentscore}(opl(i),2),handles.scale*(fact+6/8),'OPL','Color',[0 0 0.9],'FontSize',14)
    end
    % Display art
    if ~isempty(handles.score{5,handles.currentscore})
        tdebs=str2double(get(handles.currenttime,'String'))-handles.winsize/2;
        tfins=str2double(get(handles.currenttime,'String'))+handles.winsize/2;
        %---artefact on all channels
        startart=find(and(and(handles.score{5,handles.currentscore}(:,1)>tdebs,...
            handles.score{5,handles.currentscore}(:,1)<tfins),handles.score{5,handles.currentscore}(:,3)==0));
        endart=find(and(and(handles.score{5,handles.currentscore}(:,2)>tdebs,...
            handles.score{5,handles.currentscore}(:,2)<tfins),handles.score{5,handles.currentscore}(:,3)==0));
        for i=1:length(startart)
            if handles.export
                plot(ones(1,2)*handles.score{5,handles.currentscore}(startart(i),1), ...
                    [0 handles.scale*(fact+1)], ...
                    'Color',[0 0 0])
                text(handles.score{5,handles.currentscore}(startart(i),1),...
                    handles.scale*(fact+6/8), ...
                    'S.Art','Color',[0 0 0],'FontSize',14)
            else
                plot(ones(1,2)*handles.score{5,handles.currentscore}(startart(i),1), ...
                    [0 handles.scale*(fact+1)], ...
                    'UIContextMenu',handles.Deletemenu,'Color',[0 0 0])
                text(handles.score{5,handles.currentscore}(startart(i),1),...
                    handles.scale*(fact+6/8), ...
                    'S.Art','Color',[0 0 0],'FontSize',14)
            end
        end
        for i=1:length(endart)
            if handles.export
                plot(ones(1,2)*handles.score{5,handles.currentscore}(endart(i),2), ...
                    [0 handles.scale*(fact+1)], ...
                    'Color',[0 0 0])
                text(handles.score{5,handles.currentscore}(endart(i),1),...
                    handles.scale*(fact+6/8), ...
                    'E.Art','Color',[0 0 0],'FontSize',14)
            else
                plot(ones(1,2)*handles.score{5,handles.currentscore}(endart(i),2), ...
                    [0 handles.scale*(fact+1)],'UIContextMenu', ...
                    handles.Deletemenu,'Color',[0 0 0])
                text(handles.score{5,handles.currentscore}(endart(i),2),...
                    handles.scale*(fact+6/8), ...
                    'E.Art','Color',[0 0 0],'FontSize',14)
            end
        end
    end
    
    % Display arousal
    if ~isempty(handles.score{6,handles.currentscore})
        tdebs=str2double(get(handles.currenttime,'String'))-handles.winsize/2;
        tfins=str2double(get(handles.currenttime,'String'))+handles.winsize/2;
        startaro=find(and(handles.score{6,handles.currentscore}(:,1)>tdebs,...
            handles.score{6,handles.currentscore}(:,1)<tfins));
        endaro=find(and(handles.score{6,handles.currentscore}(:,2)>tdebs,...
            handles.score{6,handles.currentscore}(:,2)<tfins));
        for i=1:length(startaro)
            if handles.export
                plot(ones(1,2)*handles.score{6,handles.currentscore}(startaro(i),1), ...
                    [0 handles.scale*(fact+1)],...
                    'Color',[0 0 0])
            else
                plot(ones(1,2)*handles.score{6,handles.currentscore}(startaro(i),1), ...
                    [0 handles.scale*(fact+1)],'UIContextMenu',...
                    handles.Delar,'Color',[0 0 0])
            end
            text(handles.score{6,handles.currentscore}(startaro(i),1),...
                handles.scale*(fact+6/8), ...
                'S.Aro','Color',[1 0 0],'FontSize',14)
        end
        for i=1:length(endaro)
            if handles.export
                plot(ones(1,2)*handles.score{6,handles.currentscore}(endaro(i),2), ...
                    [0 handles.scale*(fact+1)],...
                    'Color',[0 0 0])
            else
                plot(ones(1,2)*handles.score{6,handles.currentscore}(endaro(i),2), ...
                    [0 handles.scale*(fact+1)],'UIContextMenu',...
                    handles.Delar,'Color',[0 0 0])
            end
            text(handles.score{6,handles.currentscore}(endaro(i),2),...
                handles.scale*(fact+6/8), ...
                'E.Aro','Color',[1 0 0],'FontSize',14)
        end
    end
    
    % Display EOI
    if ~isempty(handles.score{7,handles.currentscore})
        tdebs=str2double(get(handles.currenttime,'String'))-handles.winsize/2;
        tfins=str2double(get(handles.currenttime,'String'))+handles.winsize/2;
        starteoi=find(and(handles.score{7,handles.currentscore}(:,1)>tdebs,...
            handles.score{7,handles.currentscore}(:,1)<tfins));
        endeoi=find(and(handles.score{7,handles.currentscore}(:,2)>tdebs,...
            handles.score{7,handles.currentscore}(:,2)<tfins));
        for i=1:length(starteoi)
            if handles.export
                plot(ones(1,2)*handles.score{7,handles.currentscore}(starteoi(i),1), ...
                    [0 handles.scale*(fact+1)],...
                    'Color',[0 0 0])
            else
                plot(ones(1,2)*handles.score{7,handles.currentscore}(starteoi(i),1), ...
                    [0 handles.scale*(fact+1)],'UIContextMenu',...
                    handles.Deleoi,'Color',[0 0 0])
            end
            text(handles.score{7,handles.currentscore}(starteoi(i),1),...
                handles.scale*(fact+6/8), ...
                'S.EOI','Color',[0.75 0.2 0.2],'FontSize',14)
        end
        for i=1:length(endeoi)
            if handles.export
                plot(ones(1,2)*handles.score{7,handles.currentscore}(endeoi(i),2), ...
                    [0 handles.scale*(fact+1)],...
                    'Color',[0 0 0])
            else
                plot(ones(1,2)*handles.score{7,handles.currentscore}(endeoi(i),2), ...
                    [0 handles.scale*(fact+1)],'UIContextMenu',...
                    handles.Deleoi,'Color',[0 0 0])
            end
            text(handles.score{7,handles.currentscore}(endeoi(i),2),...
                handles.scale*(fact+6/8), ...
                'E.EOI','Color',[0.75 0.2 0.2],'FontSize',14)
        end
    end
end
grid on

% Display trigger
if ~handles.multcomp
    ev = events(handles.Dmeg{1});
    
    if iscell(ev)
        ev=cell2mat(ev);
    end
    if isempty(ev)
        ev=struct('time', -500,'value', -2000);
    end
    try
        [int indextoshow indextrig] = intersect(toshow,round([ev(:).time]*fsample(handles.Dmeg{1})));
    catch
        [int indextoshow indextrig] = intersect(toshow,round([ev{:}.time]*fsample(handles.Dmeg{1})));
    end
    itrigger=[];
    if ~isempty(handles.base) && ~isempty(intersect(handles.base(:,1),{ev(indextrig).type}))
        for tr = 1 : max(size(indextrig))
            itrigger = [itrigger find([ev(:).time] == ev(indextrig(tr)).time)];
        end
    end
    
    %     %Check if some events happens at the same time (je ne pense pas que
    %     %cette partie de code fonctionne!)
    %     if ~isempty(indextrig) && indextrig(1)~=1
    %         if and(ev(indextrig(1)).time == ev(indextrig(1)-1).time,length(indextrig)~=length(indextrig(1):indextrig(end)))
    %             trou=diff(indextrig);
    %             trouve=find(trou==2)+1;
    %             int = sort([int(1) int int(trouve)]);
    %             indextrig = sort([(indextrig(1)-1) indextrig (indextrig(trouve)-1)]);
    %         elseif ev(indextrig(1)).time == ev(indextrig(1)-1).time
    %             int = sort([int(1) int]);
    %             indextrig = sort([(indextrig(1)-1) indextrig]);
    %         end
    %     end
    %     if ~isempty(indextrig) && length(indextrig)~=length(indextrig(1):indextrig(end))
    %         trou=diff(indextrig);
    %         trouve=find(trou==2)+1;
    %         int = sort([int int(trouve)]);
    %         indextrig = sort([indextrig (indextrig(trouve)-1)]);
    %     end
    int = int/fsample(handles.Dmeg{1});
    Nev_dis = length(int);
    if Nev_dis % do all this if there are several triggers to be displayed !
        try
            tmp_val = ev(indextrig(1)).value;
        catch
            tmp_val = ev{indextrig(1)}.value;
        end
        if ischar(tmp_val)
            use_numv = 0;
        else
            use_numv = 1;
        end
        etype=cell(Nev_dis,1);
        for i=1:Nev_dis
            if use_numv
                etype{i}=num2str(ev(indextrig(i)).value);
            else
                etype{i}=ev(indextrig(i)).value;
            end
        end
        etpv = {ev(indextrig).type};
        plot(int,0.5*ones(1,length(int))*handles.scale/50*(NbreChandisp), ...
            'k^','LineWidth',2.5)
        if ~isempty(handles.type)
            [manev inte intm] = intersect(etpv,handles.type(:,1));
        else
            manev = [];
        end
        nme = length(manev);
        for me = 1 : nme
            inte = find(strcmpi(etpv,manev(me)));
            col = char(handles.type(intm(me),2));
            incol = col(:,1);
            if ~handles.export
                plot(int(inte),0.5*ones(1,length(int(inte)))*handles.scale/50*(NbreChandisp),'^','Color',incol,'UIContextMenu',handles.DeletEvents,...
                    'LineWidth',2.5)
                for nse = 1 : length(inte)
                    plot(ones(1,2).*int(inte(nse)),[handles.scale/2 handles.scale*(NbreChandisp)+handles.scale],'Color',incol,'UIContextMenu',handles.DeletEvents,...
                        'LineWidth',0.5)
                end
            else
                plot(int(inte),0.5*ones(1,length(int(inte)))*handles.scale/50*(NbreChandisp),'^','Color',incol,...
                    'LineWidth',2.5)
                for nse = 1 : length(inte)
                    plot(ones(1,2).*int(inte(nse)),[handles.scale/2 handles.scale*(NbreChandisp)+handles.scale],'Color',incol,'LineWidth',0.5)
                end
            end
        end
        %chose between the display of type or of value
        fmric=0;
        if ~(max(size(handles.evt))==max(size(handles.chosevt)))
            disc=1;
            evtp=get(handles.popupmenu11,'Value');
            if evtp==1
                evtch=get(handles.popupmenu10,'String');
                evtp=get(handles.popupmenu10,'Value');
                chostype=1;
                typevt=evtch(evtp);
            elseif evtp==2
                fmric=1;
                evtch=get(handles.popupmenu10,'String');
                evtp=get(handles.popupmenu10,'Value');
                chostype=1;
                typevt=evtch(evtp);
            else
                evtch=get(handles.popupmenu11,'String');
                chostype=0;
                typevt=evtch(evtp);
            end
        else
            disc=0;
            if get(handles.popupmenu11,'Value')==2
                fmric=1;
            end
        end
        for jj = 1:Nev_dis
            %if types or values are chosen by the user, only display their
            %names
            if isempty(handles.base) || isempty(intersect(handles.base(:,1),{ev(indextrig).type})) % pas de selection faite ou pas de trigger correspondant à la selection faite
                if disc && chostype && strcmpi(etpv{jj},typevt) && ~fmric
                    msg = etpv{jj};
                elseif disc && chostype && strcmpi(etpv{jj},typevt) && fmric
                    msg = etype{jj};
                elseif disc && ~(strcmpi(etype{jj},typevt)) %(disc && chostype && strcmpi(etpv{jj},typevt)) || ...
                    msg='';
                elseif disc && ~chostype && strcmpi(etype{jj},typevt) || fmric
                    msg = etype{jj};
                elseif ~disc
                    msg = num2str(etpv{jj});
                end
            else
                %msg must contain value in base 10 :
                m = find([ev(itrigger(:)).time]==ev(indextrig(jj)).time);
                [evsel itri ibase]= intersect({ev(itrigger(m)).type},handles.base(:,1));
                if isempty(evsel)
                    sumtrig = 0;
                else
                    sumtrig = sum(2.^(ibase-1));
                end
                msg = num2str(sumtrig);
            end
            %Affichage
            if isempty(msg), msg = ''; end
            lgmsg = length(msg);
            if isfield(handles.Dmeg{1},'CRC') && ...
                    isfield(handles.Dmeg{1}.CRC,'goodevents') && ...
                    size(handles.Dmeg{1}.CRC.goodevents,2)>=indextrig(jj)&& ...
                    handles.Dmeg{1}.CRC.goodevents(indextrig(jj))==0
                b=[0.8 0.2 0.2];
            else
                b = 'k';
            end
            text(int(jj)-(0.4*lgmsg/(lgmsg+1))*handles.winsize/20, ...
                2.1*handles.scale/50*NbreChandisp,msg,'Color',b);
        end
    end
end
return

%--------------------------------------------------------------------------
%-------------------------- SUBFUNCTIONS ----------------------------------
%--------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function filterlowhigh(plt,ii,handles,frqcut,scale)

if nargin<4
    flc = handles.filter.other(1)/(fsample(handles.Dmeg{ii})/2);
    fhc = handles.filter.other(2)/(fsample(handles.Dmeg{ii})/2);
    if fsample(handles.Dmeg{ii})>1500
        forder = 1;
    else
        forder = 3;
    end
    [B,A] = butter(forder,[flc,fhc],'pass');
else
    B = frqcut(1,:);
    A = frqcut(2,:);
end

X = get(plt,'YData');
% Apply Butterworth filter
Y = filtfilt(B,A,X);
set(plt,'YData',Y+scale - mean(Y))%,'Color',Col{ii})

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Y = filterforspect(handles,X,frqcut,ii)
if fsample(handles.Dmeg{ii})>1500
    forder = 1;
else
    forder = 3;
end
[B,A] = butter(forder,[frqcut(1)/(fsample(handles.Dmeg{ii})/2),...
    frqcut(2)/(fsample(handles.Dmeg{ii})/2)],'pass');

% Apply Butterworth filter
Y = filtfilt(B,A,X);

return
%%%%%%%%%%%%%%%%%%%%%%%

function cleargraph(figure,axnum)

if nargin<2
    prop    =   'Type';
    compr   =   'axes';
else
    prop    =   'Tag';
    compr   =   axnum;
end

A       =   get(figure,'Children');
idx     =   find(strcmp(get(A,prop),compr)==1);
delete(get(A(idx),'Children'))

%new function : to calculate FFT in detailed (To be checked)
function counterspect_Callback(hObject, eventdata, handles)

if get(handles.counterspect,'Value')
    fprintf('FFT detailed if on \n')
    set(handles.counterspect,'BackgroundColor',[0.5 0.5 0.5])
    cla(handles.axes5)
    set(handles.figure1,'CurrentAxes',handles.axes1)
    set(handles.Cmp_Pwr_Sp,'Enable','off','Visible','off')
    set(handles.figure1, 'windowbuttonmotionfcn',@wbdcb)
else
    fprintf('FFT detailed if off \n')
    set(handles.figure1,'CurrentAxes',handles.axes1)
    delete(findobj('tag','O'))
    delete(findobj('tag','trc'))
    set(handles.counterspect,'BackgroundColor',[0.8 0.8 0.8])
    set(handles.figure1, 'windowbuttonmotionfcn', @update_powerspect)
    set(handles.figure1, 'WindowButtonUpFcn','')
    set(handles.Cmp_Pwr_Sp,'Enable','on','Visible','on')
end

guidata(hObject, handles)

function wbdcb(hObject, eventdata)

handles = guidata(hObject);
set(handles.figure1,'CurrentAxes',handles.axes1)
if strcmp(get(hObject,'SelectionType'),'open')
    delete(findobj('tag','O'))
    set(hObject,'pointer','circle')
    set(handles.figure1,'CurrentAxes',handles.axes1)
    cp = get(handles.axes1,'CurrentPoint');
    handles.xinit = cp(1,1); handles.yinit = cp(1,2);
    plot(cp(1,1),cp(1,2),'.','color','b','tag','O');
    set(hObject,'WindowButtonMotionFcn',@wbmcb)
    set(hObject,'WindowButtonUpFcn',@wbucb)
end
guidata(hObject,handles)


function wbmcb(hObject, eventdata)

handles = guidata(hObject);
cp = get(handles.axes1,'CurrentPoint');

xinit =   handles.xinit;
yinit =   handles.yinit;

delete(findobj('tag','trc'))

xdat1 = [xinit,cp(1,1)];
xdat2 = [xinit,xinit];
xdat3 = [cp(1,1), cp(1,1)];

ydat1 = [yinit, yinit];
ydat2 = [yinit, cp(1,2)];
ydat3 = [cp(1,2), cp(1,2)];

hold on
plot(xdat1,ydat1,'tag','trc');
plot(xdat2,ydat2,'tag','trc');
plot(xdat1,ydat3,'tag','trc');
plot(xdat3,ydat2,'tag','trc');

handles.coor = [xinit cp(1,1); yinit cp(1,2)];
guidata(hObject,handles)


function wbucb(hObject,eventdata)

handles = guidata(hObject);
if strcmp(get(hObject,'SelectionType'),'normal')
    
    set(hObject,'Pointer','arrow')
    set(hObject,'WindowButtonMotionFcn',@wbdcb)
    set(handles.figure1,'CurrentAxes',handles.axes5);
    
    delete(findobj('tag', 'powerspctrm'))   %Effacer le denier power spectrum
    delete(findobj('tag', 'error'))
    delete(findobj('tag', 'localizer'))
    
    NbreChandisp    =   str2double(get(handles.NbreChan,'String'));
    Chanslidval     =   get(handles.Chanslider,'Value');
    slidpos         =   Chanslidval-rem(Chanslidval,1);
    index           =   [handles.indnomeeg handles.indexMEEG];
    chan            =   index(slidpos : 1 : slidpos + NbreChandisp -1);
    
    chandeb = floor(min(handles.coor(2,:))/handles.scale)+1;
    chanfin = floor(max(handles.coor(2,:))/handles.scale);
    chan    =   chan(chandeb : chanfin);
    fs      =   fsample(handles.Dmeg{1});
    
    tt = sort(handles.coor(1,:));
    temps = tt(1) :1/fs: tt(end);
    toshow = ceil(temps*fs);
    
    set(handles.figure1,'CurrentAxes',handles.axes5);
    Xtot = 0;
    X = 0;
    for ichan = 1 : length(chan)
        hold on
        [dumb1,dumb2,index2] = ...
            intersect(upper(chanlabels(handles.Dmeg{1},chan(ichan))),handles.names);
        if abs(handles.crc_types(index2))>1
            if handles.crc_types(index2)>0
                [dumb1,index1,dumb2] = ...
                    intersect(upper(chanlabels(handles.Dmeg{1})), ...
                    upper(handles.names(handles.crc_types(index2))));
                try
                    X   =   handles.Dmeg{1}(chan(ichan),toshow) - ...
                        handles.Dmeg{1}(index1,toshow);
                catch
                    X   = 0;
                end
            else
                range   =   max(handles.Dmeg{1}(chan(ichan),toshow)) - ...
                    min(handles.Dmeg{1}(chan(ichan),toshow));
                try
                    X  = 	(handles.scale)*handles.Dmeg{1}(chan(ichan),toshow)/range;
                catch
                    X   =   0;
                end
            end
        else
            try
                X   =   handles.Dmeg{1}(chan(ichan),toshow);
                Col	=   3;
            catch
                X   =   0;
            end
        end
        Xtot = Xtot + X;
    end
    
    if length(X) == 1
        set(handles.figure1,'CurrentAxes',handles.axes5);
        text(0.75,1, 'No Signal here')
        xlim([0 2])
        ylim([0 2])
        grid off
    else
        X       =   filterforspect(handles,X,[0.001 fs/3],1);
        [P,F]   =   pwelch(X,[],[],[],fs);
        P       =   log(P);
        p_fft   =   plot(F,P,'Color','r');
        grid on
        axis auto
        xd  = str2double(get(handles.pwrblw,'String'));
        xf = str2double(get(handles.pwrabv,'String'));
        xlim([xd xf])      %(on zoom sur ce qui nous intéresse)
        set(p_fft,'tag', 'powerspctrm')
    end
    return;
    
elseif strcmp(get(hObject,'SelectionType'),'alt')
    set(hObject,'Pointer','arrow')
    delete(findobj('tag','trc'))
    delete(findobj('tag','O'))
    set(handles.figure1, 'windowbuttonmotionfcn', @wbdcb)
    set(hObject,'WindowButtonUpFcn','')
else
    
    return
    
end

guidata(hObject,handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%New funtion : To be checked %%%%%%%%%%%%%%%%%%%%

% --- Executes on selection change in fftchan.
function fftchan_Callback(hObject, eventdata, handles)
% hObject    handle to fftchan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns fftchan contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fftchan
if handles.gui_active ~= 1,
    return;
end
set(handles.fftchan,'visible','off');
%set(handles.figure1,'Parent',handles.axes5);
axes(handles.axes5);
Chan    =   get(handles.fftchan, 'Value');
delete(findobj('tag', 'powerspctrm'))   %Effacer le denier power spectrum
delete(findobj('tag', 'error'))
delete(findobj('tag', 'localizer'))

fs      =   fsample(handles.Dmeg{1});
slidval =   get(handles.slider1,'Value');    %Prendre le pos actuelle

if handles.multcomp % To be checked : I just take this part form old version but never tried to use it
    fil     =   min(max(1,Chan),length(handles.Dmeg));
    Ctodis  =	handles.Chantodis;
    [dumb1,dumb2,Chan]  =   intersect(handles.chanset{Ctodis}, ...
        upper(chanlabels(handles.Dmeg{fil})));
    start   =   datevec(handles.date(fil,1)-handles.mindate);
    start   =   start(4)*60^2+start(5)*60+start(6);
    beg     =   slidval - start;
    tdeb    =   round(beg*fsample(handles.Dmeg{fil}));
    temps   =   tdeb:1:min(tdeb+(fsample(handles.Dmeg{fil})*handles.winsize), ...
        nsamples(handles.Dmeg{fil}));
    toshow  =   temps;
    cmap 	=   hsv(length(handles.Dmeg));
    Col     =   fil;
else
    fil     =   1;
    Chan    =   min(max(1,Chan),length(handles.indexMEEG) + length(handles.indnomeeg));
    NbreChandisp    =   str2double(get(handles.NbreChan,'String'));
    Chanslidval     =   get(handles.Chanslider,'Value');
    slidpos         =   Chanslidval-rem(Chanslidval,1);
    index           =   [handles.indnomeeg handles.indexMEEG];
    handles.inddis  =   index(slidpos : 1 : slidpos + NbreChandisp -1);
    Chan            =   handles.inddis(Chan);
    chandis         =   chanlabels(handles.Dmeg{fil},Chan);
    chtype          =   chantype(handles.Dmeg{fil},Chan);
    fs              =   fsample(handles.Dmeg{fil});
    tdeb            =   round(slidval*fs);
    temps           =   tdeb:1:min(tdeb+(fs*handles.winsize), ...
        nsamples(handles.Dmeg{fil}));
    toshow          =   temps;
    cmap            =   [0 0 0; 1 0 0; 0 0 1];
end

tdeb_w = round(slidval*fs);
tend_w = min(tdeb+(fs*handles.winsize), ...
    nsamples(handles.Dmeg{fil}));
tohid_all = [];
if ~isempty(handles.score{5,handles.currentscore})
    tdebs   =   str2double(get(handles.currenttime,'String')) - handles.winsize/2;
    tfins   =   str2double(get(handles.currenttime,'String')) + handles.winsize/2;
    art     =   find(and((or(and(handles.score{5,handles.currentscore}(:,2)>tdebs,handles.score{5,handles.currentscore}(:,2)<tfins),...
        and(handles.score{5,handles.currentscore}(:,1)<tfins,handles.score{5,handles.currentscore}(:,1)>tdebs))),...
        or(handles.score{5,handles.currentscore}(:,3) == 0,handles.score{5,handles.currentscore}(:,3) == Chan)));
    art_concerned   =   handles.score{5,handles.currentscore}(art,1:2);
    if ~isempty(art_concerned)
        a=1;
        while a <= size(art_concerned,1)
            art_concerned(a);
            begart      =   max(tdeb_w,round(art_concerned(a,1)*fs));
            endart      =   min(tend_w,round(art_concerned(a,2)*fs));
            tohid2   	=   begart : endart;
            tohid_all   =   union(tohid_all,tohid2);
            tdeb_w      =   endart;
            a   =  a+1;
        end
        toshow 	=   setdiff(toshow,tohid_all);
    end
end

%set(handles.figure1,'Parent',handles.axes5);
axes(handles.axes5);
fs      =   fsample(handles.Dmeg{fil});
leg     =   cell(0);
hold on
[dumb1,dumb2,index2] = ...
    intersect(upper(chanlabels(handles.Dmeg{fil},Chan)),handles.names);
if abs(handles.crc_types(index2))>1
    if handles.crc_types(index2)>0
        [dumb1,index1,dumb2] = ...
            intersect(upper(chanlabels(handles.Dmeg{fil})), ...
            upper(handles.names(handles.crc_types(index2))));
        try
            X   =   handles.Dmeg{fil}(Chan,toshow) - ...
                handles.Dmeg{fil}(index1,toshow);
            Col	= 1;
        catch
            X   = 0;
        end
    else
        range   =   max(handles.Dmeg{fil}(Chan,toshow)) - ...
            min(handles.Dmeg{fil}(Chan,toshow));
        try
            X   = 	(handles.scale)*handles.Dmeg{fil}(Chan,toshow)/range;
            Col =   2;
        catch
            X   =   0;
        end
    end
else
    try
        X   =   handles.Dmeg{fil}(Chan,toshow);
        Col	=   3;
    catch
        X   =   0;
    end
end

if length(X) == 1
    %set(handles.figure1,'Parent',handles.axes5);
    axes(handles.axes5);
    text(0.75,1, 'No Signal here')
    xlim([0 2])
    ylim([0 2])
    grid off
else
    %set(handles.figure1,'Parent',handles.axes5)%,'PaperSize',[20.98 29.68]);
    axes(handles.axes5);
    reset(handles.axes5);
    X       =   filterforspect(handles,X,[0.001 fs/3],fil);
    [P,F]   =   pwelch(X,[],[],[],fs);
    P       =   log(P);
    p_fft   =   plot(handles.axes5,F,P,'Color',cmap(Col,:));
    set(p_fft,'tag', 'powerspctrm')
    titre   =   chandis;
    title(handles.axes5,titre,'FontSize',12,'FontWeight','demi','FontName','Consolas')
    ylabel(handles.axes5,'Log of (power/Hz) / ','FontSize',10,'FontName','Consolas')
    xlabel(handles.axes5,'Frequency (Hz)','FontSize',10,'FontName','Consolas')
    %     if P<0
    %         set(handles.axes5,'XTick',[0 5 10 15 20],'YTick',[-60 -40 -20 0],'ZTick',[],...
    %             'XTickLabel',{'0' '5' '10' '15' '20'},'YTickLabel',{'-60' '-40' '-20' '0'},'ZTickLabel',{})
    %         axis([0 20 -60 0])
    %     else
    %         set(handles.axes5,'XTick',[0 5 10 15 20],'YTick',[0 2 4 6 8],'ZTick',[],...
    %             'XTickLabel',{'0' '5' '10' '15' '20'},'YTickLabel',{'0' '2' '4' '6' '8'},'ZTickLabel',{})
    %          axis([0 20 0 8])
    %     end
    xd = str2double(get(handles.pwrblw,'String'));
    xf = str2double(get(handles.pwrabv,'String'));
    xlim(handles.axes5,[xd xf])
    grid on
end
return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function fftchan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fftchan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%New Function : To be checked %%%%%%%%%%%%%%%%%%%%%%
function update_powerspect(hObject, eventdata, handles)
handles = guidata(hObject);
% Get the current position of the mouse pointer
set(handles.figure1,'CurrentAxes',handles.axes1)
Mouse =     get(handles.axes1, 'CurrentPoint');
x = get(handles.axes1,'XTick');
y = get(handles.axes1,'YTick');
if and(and(Mouse(1,1)<x(end),Mouse(1,1)>x(1)),and(Mouse(1,2)<y(end),Mouse(1,2)>y(1)))
    chan    =   ceil((Mouse(1,2)-handles.scale/2)/handles.scale);
    propo 	=   cellstr(get(handles.fftchan,'String'));
    chanpop =   propo(chan);
    [chanpop in] =   intersect(propo,chanpop);
    set(handles.fftchan,'Value', in);
    fftchan_Callback(hObject,eventdata,handles);
end
guidata(hObject,handles)

% function Detection_Callback(hObject, eventdata, handles)
%
% set(handles.figure1, 'windowbuttonmotionfcn', '')
%
% flags.index	=   handles.index;
% flags.Dmeg  =   handles.Dmeg;
% flags.file  =   handles.file;
% flags.winsize = handles.winsize;
% flags.user = handles.currentscore;
%
% %faire passer la taille des époques choisies pour l'analyse des artefacts
% DC_detection(flags);

function saveart(handles)

D = handles.Dmeg{1};
D.CRC.artefacteeg = handles.artefacteeg;
save(D);
D.CRC

return

%%%%%%%%%%%%%%% Changement from base 2 to 10 for the triggers  %%%%%%%%%%

% --- Executes on selection Selection in menu bar.
function selection_trigger_Callback(hObject, eventdata, handles)
% hObject    handle to Selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.figure1, 'windowbuttonmotionfcn', '')

flags.index	=   handles.index;
flags.Dmeg  =   handles.Dmeg;
flags.file  =   handles.file;
% DC_selection(flags)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- power spectrum display
function pwrabv_Callback(hObject, eventdata, handles)
% hObject    handle to pwrabv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.counterspect,'Value')
    set(handles.figure1,'CurrentAxes',handles.axes5);
    
    delete(findobj('tag', 'powerspctrm'))   %Effacer le denier power spectrum
    delete(findobj('tag', 'error'))
    delete(findobj('tag', 'localizer'))
    
    NbreChandisp    =   str2double(get(handles.NbreChan,'String'));
    Chanslidval     =   get(handles.Chanslider,'Value');
    slidpos         =   Chanslidval-rem(Chanslidval,1);
    index           =   [handles.indnomeeg handles.indexMEEG];
    chan            =   index(slidpos : 1 : slidpos + NbreChandisp -1);
    
    chandeb = floor(min(handles.coor(2,:))/handles.scale)+1;
    chanfin = floor(max(handles.coor(2,:))/handles.scale);
    chan    = chan(chandeb : chanfin);
    fs      = fsample(handles.Dmeg{1});
    
    tt = sort(handles.coor(1,:));
    temps = tt(1) :1/fs: tt(end);
    toshow = ceil(temps*fs);
    
    set(handles.figure1,'CurrentAxes',handles.axes5);
    Xtot = 0;
    X = 0;
    for ichan = 1 : length(chan)
        hold on
        [dumb1,dumb2,index2] = ...
            intersect(upper(chanlabels(handles.Dmeg{1},chan(ichan))),handles.names);
        if abs(handles.crc_types(index2))>1
            if handles.crc_types(index2)>0
                [dumb1,index1,dumb2] = ...
                    intersect(upper(chanlabels(handles.Dmeg{1})), ...
                    upper(handles.names(handles.crc_types(index2))));
                try
                    X   =   handles.Dmeg{1}(chan(ichan),toshow) - ...
                        handles.Dmeg{1}(index1,toshow);
                catch
                    X   = 0;
                end
            else
                range   =   max(handles.Dmeg{1}(chan(ichan),toshow)) - ...
                    min(handles.Dmeg{1}(chan(ichan),toshow));
                try
                    X  = 	(handles.scale)*handles.Dmeg{1}(chan(ichan),toshow)/range;
                catch
                    X   =   0;
                end
            end
        else
            try
                X   =   handles.Dmeg{1}(chan(ichan),toshow);
                Col	=   3;
            catch
                X   =   0;
            end
        end
        Xtot = Xtot + X;
    end
    
    if length(X) == 1
        set(handles.figure1,'CurrentAxes',handles.axes5);
        text(0.75,1, 'No Signal here')
        xlim([0 2])
        ylim([0 2])
        grid off
    else
        X       =   filterforspect(handles,X,[0.001 fs/3],1);
        [P,F]   =   pwelch(X,[],[],[],fs);
        P       =   log(P);
        p_fft   =   plot(F,P,'Color','r');
        grid on
        axis auto
        xd  = str2double(get(handles.pwrblw,'String'));
        xf = str2double(get(handles.pwrabv,'String'));
        xlim([xd xf])      %(on zoom sur ce qui nous intéresse)
        set(p_fft,'tag', 'powerspctrm')
    end
    return;
end

function pwrblw_Callback(hObject, eventdata, handles)
% hObject    handle to pwrblw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.counterspect,'Value')
    set(handles.figure1,'CurrentAxes',handles.axes5);
    
    delete(findobj('tag', 'powerspctrm'))   %Effacer le denier power spectrum
    delete(findobj('tag', 'error'))
    delete(findobj('tag', 'localizer'))
    
    NbreChandisp    =   str2double(get(handles.NbreChan,'String'));
    Chanslidval     =   get(handles.Chanslider,'Value');
    slidpos         =   Chanslidval-rem(Chanslidval,1);
    index           =   [handles.indnomeeg handles.indexMEEG];
    chan            =   index(slidpos : 1 : slidpos + NbreChandisp -1);
    
    chandeb = floor(min(handles.coor(2,:))/handles.scale)+1;
    chanfin = floor(max(handles.coor(2,:))/handles.scale);
    chan    =   chan(chandeb : chanfin);
    fs      =   fsample(handles.Dmeg{1});
    
    tt = sort(handles.coor(1,:));
    temps = tt(1) :1/fs: tt(end);
    toshow = ceil(temps*fs);
    
    set(handles.figure1,'CurrentAxes',handles.axes5);
    Xtot = 0;
    X = 0;
    for ichan = 1 : length(chan)
        hold on
        [dumb1,dumb2,index2] = ...
            intersect(upper(chanlabels(handles.Dmeg{1},chan(ichan))),handles.names);
        if abs(handles.crc_types(index2))>1
            if handles.crc_types(index2)>0
                [dumb1,index1,dumb2] = ...
                    intersect(upper(chanlabels(handles.Dmeg{1})), ...
                    upper(handles.names(handles.crc_types(index2))));
                try
                    X   =   handles.Dmeg{1}(chan(ichan),toshow) - ...
                        handles.Dmeg{1}(index1,toshow);
                catch
                    X   = 0;
                end
            else
                range   =   max(handles.Dmeg{1}(chan(ichan),toshow)) - ...
                    min(handles.Dmeg{1}(chan(ichan),toshow));
                try
                    X  = 	(handles.scale)*handles.Dmeg{1}(chan(ichan),toshow)/range;
                catch
                    X   =   0;
                end
            end
        else
            try
                X   =   handles.Dmeg{1}(chan(ichan),toshow);
                Col	=   3;
            catch
                X   =   0;
            end
        end
        Xtot = Xtot + X;
    end
    
    if length(X) == 1
        set(handles.figure1,'CurrentAxes',handles.axes5);
        text(0.75,1, 'No Signal here')
        xlim([0 2])
        ylim([0 2])
        grid off
    else
        X       =   filterforspect(handles,X,[0.001 fs/3],1);
        [P,F]   =   pwelch(X,[],[],[],fs);
        P       =   log(P);
        p_fft   =   plot(F,P,'Color','r');
        grid on
        axis auto
        xd  = str2double(get(handles.pwrblw,'String'));
        xf = str2double(get(handles.pwrabv,'String'));
        xlim([xd xf])      %(on zoom sur ce qui nous intéresse)
        set(p_fft,'tag', 'powerspctrm')
    end
    return;
end
