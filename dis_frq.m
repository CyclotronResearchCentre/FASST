function varargout = dis_frq(varargin)
% DIS_FRQ M-file for dis_frq.fig
%      DIS_FRQ, by itself, creates a new DIS_FRQ or raises the existing
%      singleton*.
%
%      H = DIS_FRQ returns the handle to a new DIS_FRQ or the handle to
%      the existing singleton*.
%
%      DIS_FRQ('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DIS_FRQ.M with the given input arguments.
%
%      DIS_FRQ('Property','Value',...) creates a new DIS_FRQ or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dis_frq_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dis_frq_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%__________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre
 
% Written by Y. Leclercq & C. Phillips, 2008.
% Modified by J. Schrouff, 2011.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$
 
% Edit the above text to modify the response to help dis_frq
 
% Last Modified by GUIDE v2.5 08-Jun-2011 17:02:00
 
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @dis_frq_OpeningFcn, ...
    'gui_OutputFcn',  @dis_frq_OutputFcn, ...
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
 
 
% --- Executes just before dis_frq is made visible.
function dis_frq_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dis_frq (see VARARGIN)
 
% Choose default command line output for dis_frq
set(0,'CurrentFigure',handles.figure1);
handles.output = hObject;
 
crcdef = crc_get_defaults('disfrq');
warning off
try    
    handles.file = varargin{1}.file;
    D = crc_eeg_load(handles.file);
catch   
    D = crc_eeg_load;
    handles.file = fullfile(D.path,D.fname);
end
[pth na]= fileparts(handles.file);
handles.Dmeg = D;
namesscorer = D.CRC.score(2,:);
if length(namesscorer)>1
    [nscor,v] = listdlg('PromptString','Select a scorer:',...
                      'SelectionMode','single',...
                      'ListString',namesscorer);   
else 
    nscor = 1;
end
name = [na,'_',namesscorer{1,nscor}];

set(0,'CurrentFigure',handles.figure1);
%Check if the frequency information are available
if isfield(D, 'CRC')
    if isfield(D.CRC, 'pwrspect') && isfield(D.CRC.pwrspect, 'frqbins')
        flags.Dmeg = D;
        flags.file = handles.file;
    else
        beep
        disp('Frequency information not available')
        disp(' ')
        t = timer('StartDelay',3000);
        set(t,'StartFcn',{@dis_frqcomp,hObject,[],handles});
        set(t,'TimerFcn','jj=1');
        start(t)
        wait(t);
        D = crc_eeg_load(file);
        handles.file=file;
        handles.Dmeg=D;
        guidata(hObject, handles);
    end
else
    beep
    disp('Frequency information not available')
    disp(' ')
    t = timer('StartDelay',3000);
    set(t,'StartFcn',{@dis_frqcomp,hObject,[],handles});
    set(t,'TimerFcn','jj=1');
    start(t)
    wait(t);
    D = crc_eeg_load(file);
    handles.file=file;
    handles.Dmeg=D;
    guidata(hObject, handles);
end
 
namefile = fullfile(pth,filesep,[name,'.frq']);
file = fopen(namefile);
if file<0
    namefile = fullfile(pth,filesep,[na,'.frq']);
    file = fopen(namefile);
    if file<0
        h = msgbox('There is an error')
    else 
         handles.Dmeg.CRC.pwrspect.frqdata.fname = [pth ,filesep, na '.frq'];
    end
% else
%     handles.Dmeg.CRC.pwrspect.frqdata.fname = [pth ,filesep, name '.frq'];
end
        
if isfield(handles.Dmeg,'info')
    if isfield(handles.Dmeg.info,'hour')
        handles.offset = handles.Dmeg.info.hour(1) * 60^2 + ...
            handles.Dmeg.info.hour(2)*60 + ...
            handles.Dmeg.info.hour(3);
    end
end
 
if ~isfield(handles,'offset')
    handles.offset = 0;
end
%Set parametters of visualisation
handles.scale   = crcdef.scale(1);
handles.scales  = crcdef.scale;
handles.maxpix  = crcdef.maxpix;
 
%Load reference electrodes names and types
load('CRC_electrodes.mat');
handles.names=names;
 
handles.figz = 0;
set(handles.scorerpopmenu,'Value',nscor);
handles.currentscore=get(handles.scorerpopmenu,'Value');
handles.addpoi = 1;
if ~isfield(D.CRC.pwrspect,'frq_POI')
    D.CRC.pwrspect.frq_POI=[];
end
handles.frq_POI=D.CRC.pwrspect.frq_POI;
handles.export=0;
    
%Plot on the bottom window
set(handles.figure1,'CurrentAxes',handles.axes1)
cleargraph(handles.figure1)
 
% Determine which channel is displayed
cindex = meegchannels(handles.Dmeg);
handles.spectrochan = cindex(1);
cindex = cindex(1);
 
%Plot Blue line
fminblue=0.5;
fmaxblue=4;
handles.bluelines=[];
handles.redlines = [];
handles.fminblue=fminblue;
handles.fmaxblue=fmaxblue;
frqtoplot = find(handles.Dmeg.CRC.pwrspect.frqbins >= ...
    fminblue & handles.Dmeg.CRC.pwrspect.frqbins<=fmaxblue);
temps = handles.Dmeg.CRC.pwrspect.duration/2: ...
    handles.Dmeg.CRC.pwrspect.step: ...
    size(handles.Dmeg.CRC.pwrspect.frqdata,3)*handles.Dmeg.CRC.pwrspect.step;
cla(handles.axes1,'reset')
cla(handles.spectro,'reset')
cla(handles.Hypnogra,'reset')
 
 
handles.subsmpl = crcdef.subsmpl; %smoothing factor
 
%Plot data
for i=1:length(cindex)
    hold on
    if isempty(frqtoplot)
        vectortoplot=zeros(1,length(temps));
    else
        vectortoplot=mean(squeeze(handles.Dmeg.CRC.pwrspect.frqdata(cindex(i), ...
            frqtoplot,1:length(temps))));
    end
    vectortoplot=log(vectortoplot);
    if handles.subsmpl==0;
        %compte the subsampling factor according to the size of the file
        st=length(find(vectortoplot~=-Inf));
        if st<100
            handles.subsmpl=1;
        elseif st>=100 && st<500
            handles.subsmpl=4;
        elseif st>=500 && st<1000
            handles.subsmpl=8;
        else
            handles.subsmpl=16;
        end
    end
    subsmpl = handles.subsmpl;
    toadd = mod(length(vectortoplot),subsmpl);
    vectortoplot = [vectortoplot zeros(1,subsmpl-toadd)];
    vectortoplot = reshape(vectortoplot,subsmpl,length(vectortoplot)/subsmpl);
    vectortoplot = mean(vectortoplot);
    temps = temps(1:subsmpl:length(temps));
    if length(vectortoplot)>length(temps)
        tps = temps(length(temps))-temps(length(temps)-1);
        temps = [temps (temps(length(temps))+tps)];
    end
    numline=plot(temps,handles.scale*i+vectortoplot);
    handles.bluelines=[handles.bluelines numline];
end
ylim([handles.scale/2 handles.scale*(i+1)])
xlim([handles.Dmeg.CRC.pwrspect.duration/2-1  ...
    nsamples(handles.Dmeg)/(fsample(handles.Dmeg)+1)])
set(handles.bluelines,'Visible','off')
 
%Define the label/ticks of the axes
set(handles.axes1,'YTick',[handles.scale:handles.scale/2:i*handles.scale+handles.scale/2]);
ylabels=[];
for j=1:length(cindex)
    q=cell(1);
    q{1}='0';
    ylabels=[ylabels q];
    ylabels=[ylabels num2str(round(handles.scale/2))];
end
set(handles.axes1,'YTickLabel',ylabels);
 
%Plot the spectrogram
set(handles.figure1,'CurrentAxes',handles.spectro)
picstep = ceil(size(handles.Dmeg.CRC.pwrspect.frqdata,3)/handles.maxpix);
toshowpic = 1:picstep:size(handles.Dmeg.CRC.pwrspect.frqdata,3);
Spectr(:,:)=handles.Dmeg.CRC.pwrspect.frqdata(cindex,2:end-1,toshowpic);
hspec=imagesc(log(Spectr));
set(hspec,'UIContextMenu',handles.POI_def); %à voir
colormap(jet)
axis xy
 
% Ticks
tim = ceil(nsamples(handles.Dmeg)/(fsample(handles.Dmeg)));
Xtic=(sort((tim):-(tim)/10:0));
[dumb, times]=crc_time_converts(mod(Xtic+handles.offset,24*60^2));
tim = ceil(nsamples(handles.Dmeg)/(fsample(handles.Dmeg)*picstep));
Xtic=(sort((tim):-(tim)/10:0));
set(handles.spectro,'XTick',Xtic/handles.Dmeg.CRC.pwrspect.step);
set(handles.spectro,'XTickLabel',times)
 
List=get(handles.spectro,'YTick');
try
    set(handles.spectro,'YTickLabel',round(handles.Dmeg.CRC.pwrspect.frqbins(List)));
catch
    % Sometimes (Matlab version issue ?) 'List' is badly initialized...
    % It consists in a list 0 .1 .2 .3 .4 .5 .6 .7 .8 9 1.0
    % so initialising it as it should reasonably be
    Nfreqbins = length(handles.Dmeg.CRC.pwrspect.frqbins);
    List = (1:(round(Nfreqbins/10)-1))*10;
    set(handles.spectro,'YTickLabel',round(handles.Dmeg.CRC.pwrspect.frqbins(List)));
end
 
%Plot Periods of Interest
a=get(gca,'Xlim');
b=get(handles.axes1,'Xlim');
scl=(b(2)-b(1))/(a(2)-a(1));
a=get(gca,'Ylim');
hold on
handles.numpoi=[];
if ~isempty(handles.frq_POI)
    npoi=size(handles.frq_POI,1);
    for i=1:npoi
        numpoi=plot(ones(1,2)*(handles.frq_POI(i,1)/scl),[0 a(2)], ...
             'UIContextMenu',handles.Deletemenu,'Color',[0 0 0]);
        text(handles.frq_POI(i,1)/scl,a(2)*7/8, ...
             'S.POI','Color',[0 0 0],'FontSize',14)
         handles.numpoi=[handles.numpoi, numpoi];
    end
    for i=1:npoi        
        numpoi=plot(ones(1,2)*(handles.frq_POI(i,2)/scl),[0 a(2)], ...
             'UIContextMenu',handles.Deletemenu,'Color',[0 0 0]);
        text(handles.frq_POI(i,2)/scl,a(2)*7/8, ...
             'E.POI','Color',[0 0 0],'FontSize',14)
         handles.numpoi=[handles.numpoi, numpoi];
    end
end
 
%Plot the hypnogram
 
if isfield(handles.Dmeg,'CRC')
    if isfield(handles.Dmeg.CRC,'score')
        def=crc_get_defaults('score');
        strsw=def.stnames_S(1:6);
        handles.score=handles.Dmeg.CRC.score;
        handles.scoring = true;
        strsw=[strsw,{'Slow Wave Sleep (3-4)'},{'NREM Sleep (2-4)'}];
        set(handles.figure1,'CurrentAxes',handles.Hypnogra)
        crc_hypnoplot(handles.Hypnogra,...
            handles,...
            handles.Dmeg.CRC.score{3,nscor},...
            handles.Dmeg.CRC.score{1,nscor},...
            handles.Dmeg.CRC.score{2,nscor});
        set(handles.scorerpopmenu,'String',handles.Dmeg.CRC.score(2,:))
        set(handles.scorerpopmenu,'Value',nscor)
        set(handles.popstage,'String',strsw,'Value',1)
    else
        set(handles.figure1,'CurrentAxes',handles.Hypnogra)
        text(0.45,0.5,'NO HYPNOGRAM')
        set(handles.Hypnogra,'YTick',[]);
        set(handles.Hypnogra,'XTick',[]);
        set(handles.scorerpopmenu,'Visible','off')
        set(handles.scorertxt,'Visible','off')
        set(handles.meanpwrspect,'Visible','off')
        set(handles.popstage,'Visible','off')
        set(handles.pwrstage,'Visible','off')
    end
else
    set(handles.figure1,'CurrentAxes',handles.Hypnogra)
    text(0.45,0.5,'NO HYPNOGRAM')
    set(handles.Hypnogra,'YTick',[]);
    set(handles.Hypnogra,'XTick',[]);
    set(handles.scorerpopmenu,'Visible','off')
    set(handles.scorertxt,'Visible','off')
    set(handles.meanpwrspect,'Visible','off')
    set(handles.popstage,'Visible','off')
    set(handles.pwrstage,'Visible','off')
end
 
% Ticks
tim = ceil(nsamples(handles.Dmeg)/fsample(handles.Dmeg));
Xtic=(sort((tim):-(tim)/10:0));
Xtick = mod(Xtic + handles.offset,24*60^2);
[dumb, times]=crc_time_converts(Xtick);
set(handles.axes1,'XTick',Xtic);
set(handles.axes1,'XTickLabel',times)
 
 
% Ticks on spectrogram
tim = ceil(nsamples(handles.Dmeg)/fsample(handles.Dmeg));
Xtic=(sort((tim):-(tim)/10:0));
[dumb, times]=crc_time_converts(mod(Xtic+handles.offset,24*60^2));
% tim = ceil(nsamples(handles.Dmeg)/(fsample(handles.Dmeg)*picstep));
% Xtic=(sort((tim):-(tim)/10:0));
Xtic=get(handles.axes1,'XTick');
a=get(handles.spectro,'Xlim');
b=get(handles.axes1,'Xlim');
scl=abs((b(2)-b(1))/(a(2)-a(1)));
set(handles.spectro,'XTick',Xtic/scl);
% set(handles.spectro,'XTick',Xtic/handles.Dmeg.CRC.pwrspect.step);
set(handles.spectro,'XTickLabel',times)
 
 
%Set the other elements
set(handles.selchanspectro,'String',upper(chanlabels(handles.Dmeg)));
scdmenu = ['None' upper(chanlabels(handles.Dmeg))];
set(handles.secondchan,'String',scdmenu);
set(handles.filename,'String',handles.Dmeg.fname);
handles.spectrochan=1;
 
%Set the visibility of elements
set(handles.axes1,'Visible','off')
set(handles.Blue,'Visible','off')
set(handles.fmin_blue,'Visible','off')
set(handles.blue_fmax,'Visible','off')
set(handles.text10,'Visible','off')
set(handles.text14,'Visible','off')
set(handles.text2,'Visible','off')
set(handles.Scale,'Visible','off')
set(handles.secondchan,'Visible','off')
set(handles.txt2ndchan,'Visible','off')
set(handles.bluelines,'Visible','off')
set(handles.redlines,'Visible','off')
set(handles.Pwr_RA,'Visible','off')
set(handles.numpoi,'Visible','on');
 
 
% Set the value of the radiobutton
set(handles.dis_fband,'Value',0);
set(handles.dis_spectr,'Value',1);
 
% Allow Mongrain Representation for Frequency band representation
try
    handles.Dmeg.CRC.score;
catch
    labels=get(handles.Pwr_RA,'String');
    set(handles.Pwr_RA,'String',labels(1:2));
end
 
% Update handles structure
guidata(hObject, handles);
 
 
% --- Outputs from this function are returned to the command line.
function varargout = dis_frq_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Get default command line output from handles structure
varargout{1} = handles.output;
 
 
 
function Scale_Callback(hObject, eventdata, handles)
% hObject    handle to Scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: get(hObject,'String') returns contents of Scale as text
%        str2double(get(hObject,'String')) returns contents of Scale as a double
 
if not(str2double(get(hObject,'String'))>0)&& not(str2double(get(hObject,'String'))<0)
    set(handles.Scale,'String',num2str(handles.scale));
else
    handles.scale=str2double(get(hObject,'String'));
end
 
handles = mainplot_frq(handles);
 
% Update handles structure
guidata(hObject, handles);
 
 
% --- Executes during object creation, after setting all properties.
function Scale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Scale (see GCBO)
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
 
D = crc_eeg_load;
file = fullfile(D.path,D.fname);
handles.file=file;
dis_frq(handles);
 
function fmin_blue_Callback(hObject, eventdata, handles)
% hObject    handle to fmin_blue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: get(hObject,'String') returns contents of fmin_blue as text
%        str2double(get(hObject,'String')) returns contents of fmin_blue as a double
 
handles.fminblue=str2double(get(hObject,'String'));
 
handles = mainplot_frq(handles);
 
% Update handles structure
guidata(hObject, handles);
 
% --- Executes during object creation, after setting all properties.
function fmin_blue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fmin_blue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
function blue_fmax_Callback(hObject, eventdata, handles)
% hObject    handle to blue_fmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: get(hObject,'String') returns contents of blue_fmax as text
%        str2double(get(hObject,'String')) returns contents of blue_fmax as a double
 
handles.fmaxblue=str2double(get(hObject,'String'));
 
handles = mainplot_frq(handles);
 
% Update handles structure
guidata(hObject, handles);
 
 
% --- Executes during object creation, after setting all properties.
function blue_fmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blue_fmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
 
% --- Executes on selection change in selchanspectro.
function selchanspectro_Callback(hObject, eventdata, handles)
% hObject    handle to selchanspectro (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: contents = get(hObject,'String') returns selchanspectro contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selchanspectro
 
contents = get(hObject,'String');
selchan=upper(contents{get(hObject,'Value')});
[dumb1,index]=intersect(upper(chanlabels(handles.Dmeg)),selchan);
handles.spectrochan=index;
 
% mainplot_frq(handles);
 
toshow=1:size(handles.Dmeg.CRC.pwrspect.frqdata,3);
temps = handles.Dmeg.CRC.pwrspect.duration/2: ...
    handles.Dmeg.CRC.pwrspect.step: ...
    size(handles.Dmeg.CRC.pwrspect.frqdata,3)*handles.Dmeg.CRC.pwrspect.step;
 
%Plot the spectrogram
if get(handles.dis_spectr,'Value')
 
    set(handles.figure1,'CurrentAxes',handles.spectro)
    toshow=toshow(1:length(toshow)-1);
    picstep = ceil(size(handles.Dmeg.CRC.pwrspect.frqdata,3)/handles.maxpix);
    toshowpic = 1:picstep:size(handles.Dmeg.CRC.pwrspect.frqdata,3);
    Spectr(:,:)=handles.Dmeg.CRC.pwrspect.frqdata(handles.spectrochan,2:end-1,toshowpic);
    hspec=imagesc(log(Spectr));
    set(hspec,'UIContextMenu',handles.POI_def);
    colormap(jet)
    axis xy
    step=round(size(handles.Dmeg.CRC.pwrspect.frqdata,3)/10);
 
    % Ticks
    tim = ceil(nsamples(handles.Dmeg)/(fsample(handles.Dmeg)));
    Xtic=(sort((tim):-(tim)/10:0));
    [dumb, times]=crc_time_converts(mod(Xtic+handles.offset,24*60^2));
%     tim = ceil(nsamples(handles.Dmeg)/(fsample(handles.Dmeg)*picstep));
%     Xtic=(sort((tim):-(tim)/10:0));
    Xtic=get(handles.axes1,'XTick');
    a=get(handles.spectro,'Xlim');
    b=get(handles.axes1,'Xlim');
    scl=abs((b(2)-b(1))/(a(2)-a(1)));
    a=get(handles.spectro,'Ylim');
    set(handles.spectro,'XTick',Xtic/scl);
    % set(handles.spectro,'XTick',Xtic/handles.Dmeg.CRC.pwrspect.step);
    set(handles.spectro,'XTickLabel',times)
    List=get(handles.spectro,'YTick');
    set(handles.spectro,'YTickLabel',round(handles.Dmeg.CRC.pwrspect.frqbins(List)));
    % % %Plot Periods of Interest
    if ~isempty(handles.numpoi)
        hold on
        npoi=size(handles.frq_POI,1);
        for i=1:npoi
            if handles.export
                plot(ones(1,2)*(handles.frq_POI(i,1)/scl),[0 a(2)], ...
                    'Color',[0 0 0]);
            else
                plot(ones(1,2)*(handles.frq_POI(i,1)/scl),[0 a(2)], ...
                    'UIContextMenu',handles.Deletemenu,'Color',[0 0 0]);
            end
            text((handles.frq_POI(i,1)/scl),a(2)*(7/8), ...
                 'S.POI','Color',[0 0 0],'FontSize',14)
        end
        for i=1:npoi
            if handles.export
                plot(ones(1,2)*(handles.frq_POI(i,2)/scl),[0 a(2)], ...
                    'Color',[0 0 0]);
            else
               plot(ones(1,2)*(handles.frq_POI(i,2)/scl),[0 a(2)], ...
                    'UIContextMenu',handles.Deletemenu,'Color',[0 0 0]);
            end
             text(handles.frq_POI(i,2)/scl,a(2)*(7/8), ...
                 'E.POI','Color',[0 0 0],'FontSize',14)
        end
    end
 
end
 
if get(handles.dis_fband,'Value')
    handles = mainplot_frq(handles);
end
 
guidata(hObject, handles);
 
% --- Executes during object creation, after setting all properties.
function selchanspectro_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selchanspectro (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
 
 
 
% --- Executes on selection change in secondchan.
function secondchan_Callback(hObject, eventdata, handles)
% hObject    handle to secondchan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: contents = get(hObject,'String') returns secondchan contents as cell array
%        contents{get(hObject,'Value')} returns selected item from secondchan
 
handles = mainplot_frq(handles);
 
% Update handles structure
guidata(hObject, handles);
 
 
 
% --- Executes during object creation, after setting all properties.
function secondchan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to secondchan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
function dis_spectr_Callback(hObject, eventdata, handles)
 
% Set the value of the radiobutton
set(handles.dis_spectr,'Value',1);
set(handles.dis_fband,'Value',0);
 
%set the visibility of the elements
 
set(handles.spectro,'Visible','on')
set(get(handles.spectro,'Children'),'Visible','on')
set(handles.axes1,'Visible','off')
set(handles.Blue,'Visible','off')
set(handles.text2,'Visible','off')
set(handles.Scale,'Visible','off')
set(handles.secondchan,'Visible','off')
set(handles.txt2ndchan,'Visible','off')
try
    for i=1:length(handles.numpoi)
        set(handles.numpoi(i),'Visible','off')
    end
end
try
    set(handles.bluelines,'Visible','off')
    set(handles.redlines,'Visible','off')   
end
 
set(handles.Blue,'Visible','off')
set(handles.fmin_blue,'Visible','off')
set(handles.blue_fmax,'Visible','off')
set(handles.text10,'Visible','off')
set(handles.text14,'Visible','off')
set(handles.Pwr_RA,'Visible','off')
 
contents = get(handles.selchanspectro,'String');
selchan=upper(contents{get(handles.selchanspectro,'Value')});
[dumb1,index]=intersect(upper(chanlabels(handles.Dmeg)),selchan);
handles.spectrochan=index;
 
if ~handles.export
    set(handles.figure1,'CurrentAxes',handles.spectro)
    cla(handles.spectro)
    cla(handles.axes1)
end

toshow=1:size(handles.Dmeg.CRC.pwrspect.frqdata,3);
temps = handles.Dmeg.CRC.pwrspect.duration/2: ...
    handles.Dmeg.CRC.pwrspect.step: ...
    size(handles.Dmeg.CRC.pwrspect.frqdata,3)*handles.Dmeg.CRC.pwrspect.step;
 
%Plot the spectrogram

toshow=toshow(1:length(toshow)-1);
picstep = ceil(size(handles.Dmeg.CRC.pwrspect.frqdata,3)/handles.maxpix);
toshowpic = 1:picstep:size(handles.Dmeg.CRC.pwrspect.frqdata,3);
Spectr(:,:)=handles.Dmeg.CRC.pwrspect.frqdata(handles.spectrochan,2:end-1,toshowpic);
hspec=imagesc(log(Spectr));
if ~handles.export
   set(hspec,'UIContextMenu',handles.POI_def);
end
colormap(jet)
axis xy
step=round(size(handles.Dmeg.CRC.pwrspect.frqdata,3)/10);
 
% Ticks
tim =ceil(nsamples(handles.Dmeg)/(fsample(handles.Dmeg)));
Xtic=(sort((tim):-(tim)/10:0));
[dumb, times]=crc_time_converts(mod(Xtic+handles.offset,24*60^2));
% tim = ceil(nsamples(handles.Dmeg)/(fsample(handles.Dmeg)*picstep));
% Xtic=(sort((tim):-(tim)/10:0));
Xtic=get(handles.axes1,'XTick');
a=get(gca,'Xlim');
b=get(handles.axes1,'Xlim');
scl=(b(2)-b(1))/(a(2)-a(1));
set(gca,'XTick',Xtic/scl);
% set(gca,'XTick',Xtic/handles.Dmeg.CRC.pwrspect.step);
set(gca,'XTickLabel',times)
List=get(gca,'YTick');
set(gca,'YTickLabel',round(handles.Dmeg.CRC.pwrspect.frqbins(List)));
a=get(gca,'Ylim');
% 
% % %Plot Periods of Interest
if ~isempty(handles.numpoi)
    hold on
    npoi=size(handles.frq_POI,1);
    for i=1:npoi
        if handles.export
            plot(ones(1,2)*(handles.frq_POI(i,1)/scl),[0 a(2)], ...
                'Color',[0 0 0]);
        else
            plot(ones(1,2)*(handles.frq_POI(i,1)/scl),[0 a(2)], ...
                'UIContextMenu',handles.Deletemenu,'Color',[0 0 0]);
        end
        text((handles.frq_POI(i,1)/scl),a(2)*(7/8), ...
             'S.POI','Color',[0 0 0],'FontSize',14)
    end
    for i=1:npoi
        if handles.export
            plot(ones(1,2)*(handles.frq_POI(i,2)/scl),[0 a(2)], ...
                'Color',[0 0 0]);
        else
           plot(ones(1,2)*(handles.frq_POI(i,2)/scl),[0 a(2)], ...
                'UIContextMenu',handles.Deletemenu,'Color',[0 0 0]);
        end
         text(handles.frq_POI(i,2)/scl,a(2)*(7/8), ...
             'E.POI','Color',[0 0 0],'FontSize',14)
    end
end
 
% Update handles structure
guidata(hObject, handles);
 
return
 
 
function dis_fband_Callback(hObject, eventdata, handles)
 
% Correct radiobutton values
set(handles.dis_spectr,'Value',0);
set(handles.dis_fband,'Value',1);
 
 
%set visibility
set(handles.spectro,'Visible','off')
set(get(handles.spectro,'Children'),'Visible','off')
set(handles.axes1,'Visible','on')
set(handles.Blue,'Visible','on')
set(handles.text2,'Visible','on')
set(handles.Scale,'Visible','on')
set(handles.secondchan,'Visible','on')
set(handles.txt2ndchan,'Visible','on')
try
    set(handles.redlines,'Visible','on')
end
set(handles.Blue,'Visible','on')
set(handles.fmin_blue,'Visible','on')
set(handles.blue_fmax,'Visible','on')
set(handles.text10,'Visible','on')
set(handles.text14,'Visible','on')
set(handles.Pwr_RA,'Visible','on')
 
%Plot Elements
handles = mainplot_frq(handles);
 
% Update handles structure
guidata(hObject, handles);
 
return
 
 
% --- Executes on selection change in scorerpopmenu.
function scorerpopmenu_Callback(hObject, eventdata, handles)
% hObject    handle to scorerpopmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: contents = get(hObject,'String') returns scorerpopmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from scorerpopmenu
 
Val=get(handles.scorerpopmenu,'Value');
Str=get(handles.scorerpopmenu,'String');
 
handles.currentscore=Val;
 
if get(handles.Pwr_RA,'Value')==3
    mainplot_frq(handles);
end
 
guidata(hObject, handles);
set(handles.figure1,'CurrentAxes',handles.Hypnogra)
 
crc_hypnoplot(handles.Hypnogra,...
    handles,...
    handles.Dmeg.CRC.score{3,handles.currentscore},...
    handles.Dmeg.CRC.score{1,handles.currentscore},...
    handles.Dmeg.CRC.score{2,handles.currentscore});
 
[pth na]= fileparts(handles.file); 
name = [na,'_',char(Str(Val))];
namefile = fullfile(pth,filesep,[name,'.frq']);
file = fopen(namefile);
if file<0
    handles.Dmeg.CRC.pwrspect.frqdata.fname = spm_select(Inf, 'any', 'Select imported frq file associated','' ...
            ,pwd,'\.[fF][rR][qQ]');
else
    handles.Dmeg.CRC.pwrspect.frqdata.fname = [pth ,filesep, name '.frq'];%C'est ici qu'on vient chercher le fichier frq associé
end
if get(handles.dis_spectr,'Value')
    dis_spectr_Callback(hObject, eventdata, handles)
else
    dis_fband_Callback(hObject, eventdata, handles)
end
 
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
 
% --- Executes on selection change in popstage.
function popstage_Callback(hObject, eventdata, handles)
% hObject    handle to popstage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: contents = get(hObject,'String') returns popstage contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popstage
 
 
% --- Executes during object creation, after setting all properties.
function popstage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popstage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
 
% --- Executes on button press in pwrstage.
function pwrstage_Callback(hObject, eventdata, handles)
% hObject    handle to pwrstage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
contents = get(handles.popstage,'String');
stagestr = contents{get(handles.popstage,'Value')};
 
contents = get(handles.selchanspectro,'String');
selchan=upper(contents{get(handles.selchanspectro,'Value')});
[dumb1,channumb]=intersect(upper(chanlabels(handles.Dmeg)),selchan);
 
numscore = get(handles.scorerpopmenu,'Value');
stage = get(handles.popstage,'Value');
 
% Find the moment of interest from the score
switch  stage
    case {1,2,3,4,5,6}
        score = find(handles.Dmeg.CRC.score{1,numscore} == stage-1); 
        % Difference of one between the stage nr and place in the popupmenu
    case {7}
        % between 3 and 4
        score = find(and(handles.Dmeg.CRC.score{1,numscore} > 2, ... 
            handles.Dmeg.CRC.score{1,numscore} < 5));
    case {8}
        % between 2 and 4
        score = find(and(handles.Dmeg.CRC.score{1,numscore} >1, ...  
            handles.Dmeg.CRC.score{1,numscore} < 5 ));
end
 
% Computing the time index of all frequency information available
timevect = handles.Dmeg.CRC.pwrspect.step:handles.Dmeg.CRC.pwrspect.step: ...
    size(handles.Dmeg.CRC.pwrspect.frqdata,3)*handles.Dmeg.CRC.pwrspect.step;

% Computing the time index of the window of interest
timewin = [(score*20 -20) ; score*20 ];
 
% Find the time index of the frequency info that matches the score of
% interest
timofint=[];
for i = 1:size(timewin, 2)
    tmp=find(timevect > timewin(1,i) & timevect <= timewin(2,i));
    timofint = [timofint tmp];
end
 
% Find the time index of the frequency info between FPL and OPL
% Use by default the 1st user's score !!!
try
    betpl = find(timevect > handles.Dmeg.CRC.score{4,1}(1) & ...
        timevect < handles.Dmeg.CRC.score{4,1}(2));
catch
    betpl = find(timevect > handles.Dmeg.CRC.pl(1) & ...
        timevect < handles.Dmeg.CRC.pl(2));
end
% Intersection of the two
inter = intersect(betpl,timofint);
 
if handles.figz~=0
    z=handles.figz;
else
    z=figure;
end
 
figure(z)
set(z,'NumberTitle','off')
set(z,'Name','Mean Power Spectrum')
subplot(1,1,1)
try
    semilogy(handles.Dmeg.CRC.pwrspect.frqbins, ...
        mean(squeeze(handles.Dmeg.CRC.pwrspect.frqdata(channumb,:,inter))'))
 
    grid on
    titre=['Power on ' chanlabels(handles.Dmeg,channumb) ' during ' stagestr];
    title(titre)
    ylabel('Log of power')
    xlabel('Frequency in Hz')
    xlim([0 20])
end
handles.figz=z;
 
% Update handles structure
guidata(hObject, handles);
 
 
% --- Executes on selection change in Pwr_RA.
function Pwr_RA_Callback(hObject, eventdata, handles)
% hObject    handle to Pwr_RA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
 
if get(handles.Pwr_RA,'Value')==1
    set(handles.text2,'String','Scale (dB)')
    set(handles.Scale,'String',num2str(handles.scales(1)))
elseif get(handles.Pwr_RA,'Value')==2
        set(handles.text2,'String','Scale (no unit)')
        set(handles.Scale,'String',num2str(handles.scales(2)))
else
    set(handles.text2,'String','Scale (no unit)')
    set(handles.Scale,'String',num2str(handles.scales(3)))
end
 
handles.scale = str2double(get(handles.Scale,'String'));
handles = mainplot_frq(handles);
 
% Update handles structure
guidata(hObject, handles);
 
 
% --- Executes during object creation, after setting all properties.
function Pwr_RA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pwr_RA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
% --------------------------------------------------------------------
%Menu to export the main axes in a Matlab figure
function export_figure_Callback(hObject, eventdata, handles)
% hObject    handle to export_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.export=1;
z=figure;
axs=axes;
set(z,'CurrentAxes',axs)
handles.currentfigure=z;
if get(handles.dis_spectr,'Value')==1
    dis_spectr_Callback(hObject, eventdata, handles)
elseif get(handles.dis_fband,'Value')==1
    dis_fband_Callback(hObject, eventdata, handles)
end
ha=get(handles.currentfigure,'CurrentAxes');
index=max(length(handles.spectrochan),1);
 
%Define the label/ticks of the axes
%display y-labels
set(ha,'YTick',[handles.scale:handles.scale/2:index*handles.scale+handles.scale/2]);
ylabels=[];
for j=1:index
    q=cell(1);
    q{1}='0';
    ylabels=[ylabels q];
    ylabels=[ylabels num2str(round(handles.scale*100/2)/100)];
end
set(ha,'YTickLabel',ylabels);
 
%display x-labels
Xtick=get(handles.axes1,'XTick');
a=get(ha,'Xlim');
b=get(handles.axes1,'Xlim');
scl=abs((b(2)-b(1))/(a(2)-a(1)));
set(ha,'XTick',Xtick/scl);
% xtick = get(handles.axes1,'XTick');
if isfield(handles,'offset')
    Xtick = mod(Xtick + handles.offset,24*60^2);
end
[time string] = crc_time_converts(Xtick);
set(ha,'XTickLabel',string)
 
handles.figz=z;
handles.export=0;
 
% Update handles structure
guidata(hObject, handles);
 
return
 
 
% --------------------------------------------------------------------
%Right-click menu to mark a period of interest
function StartPOI_Callback(hObject, eventdata, handles)
% hObject    handle to POI_def (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Mouse=get(handles.axes1,'CurrentPoint');
timedbtart=Mouse(1,1);
 
if get(handles.dis_spectr,'Value')==1
    a=get(gca,'Xlim');
    b=get(handles.axes1,'Xlim');
    scl=(b(2)-b(1))/(a(2)-a(1));
elseif get(handles.dis_fband,'Value')==1
    scl=1;
end
a=get(gca,'Ylim');
 
if handles.addpoi== 1
    if ~isempty(handles.frq_POI)
        if sum(and(timedbtart>handles.frq_POI(:,1),...
                timedbtart<handles.frq_POI(:,2)))
            beep
            disp('Invalid "start Period of Interest" point')
            disp(' ')
 
        else
            handles.frq_POI= ...
                [handles.frq_POI; timedbtart timedbtart];
            set(handles.StartPOI,'Label','Add "End POI point"');
            handles.addpoi= 0;
            plot(ones(1,2)*(timedbtart/scl),[0 a(2)], ...
                'UIContextMenu',handles.Deletemenu,'Color',[0 0 0])
            text(timedbtart/scl,a(2)*7/8, ...
                'S.POI','Color',[0 0 0],'FontSize',14)
        end
    else
        handles.frq_POI=...
            [handles.frq_POI; timedbtart timedbtart];
        set(handles.StartPOI,'Label','Add "End POI point"');
        handles.addpoi = 0;
        plot(ones(1,2)*(timedbtart/scl),[0 a(2)], ...
            'UIContextMenu',handles.Deletemenu,'Color',[0 0 0])
        text(timedbtart/scl,a(2)*7/8, ...
            'S.POI','Color',[0 0 0],'FontSize',14)
    end
else
    Mouse=get(handles.axes1,'CurrentPoint');
    timefinart=Mouse(1,1);
    [row,col]=find(and(timefinart>handles.frq_POI(:,1),...
        timefinart>handles.frq_POI(:,2)));
    test3=sum(and(handles.frq_POI...
        (size(handles.frq_POI,1),1)<...
        handles.frq_POI(row,1), ...
        handles.frq_POI...
        (size(handles.frq_POI,1),2)<...
        handles.frq_POI(row,2)));
    if or(or(sum(and(timefinart>handles.frq_POI(:,1),...
            timefinart<handles.frq_POI(:,2)))>0, ...
            timefinart<handles.frq_POI...
            (size(handles.frq_POI,1),2)),test3)
        beep
        disp('Invalid "end Period of Interest" point')
        disp(' ')
    else
        set(handles.StartPOI,'Label','Add "Start POI point"');
        handles.frq_POI...
            (size(handles.frq_POI,1),2)= timefinart;
        handles.addpoi = 1;
        plot(ones(1,2)*(timedbtart/scl),[0 a(2)], ...
            'UIContextMenu',handles.Deletemenu,'Color',[0 0 0])
        text(timedbtart/scl,a(2)*7/8, ...
            'E.POI','Color',[0 0 0],'FontSize',14)
    end
end
handles.Dmeg.CRC.pwrspect.frq_POI=handles.frq_POI;
D=handles.Dmeg;
save(D);
 
% Update handles structure
guidata(hObject, handles);
 
% --------------------------------------------------------------------
function Deletemenu_Callback(hObject, eventdata, handles)
% hObject    handle to Deletemenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% --------------------------------------------------------------------
function Delpoi_Callback(hObject, eventdata, handles)
% hObject    handle to Delart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Mouse=get(handles.axes1,'CurrentPoint');
minimum=min(min(abs(handles.frq_POI-Mouse(1,1))));
[row,col]=find((abs(handles.frq_POI-Mouse(1,1))-minimum)==0);
 
if length(row)==2
    handles.addpoi = 1;
    set(handles.StartPOI,'Label','Add "Start POI point"');
end
 
handles.frq_POI = ...
    [handles.frq_POI(1:row-1,:) ; ...
    handles.frq_POI(row+1:size(handles.frq_POI,1),:)];
 
%Save the changes
handles.Dmeg.CRC.pwrspect.frq_POI=handles.frq_POI;
D=handles.Dmeg;
save(D);
 
if get(handles.dis_spectr,'Value')==1
    dis_spectr_Callback(hObject, eventdata, handles)
elseif get(handles.dis_fband,'Value')==1
    mainplot_frq(handles)
end
 
 
% Update handles structure
guidata(hObject, handles);
 
% --------------------------------------------------------------------
function POI_def_Callback(hObject, eventdata, handles)
% hObject    handle to infoPOI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
 
% --------------------------------------------------------------------
%right click menu to obtain more info on the selected POI
function infoPOI_Callback(hObject, eventdata, handles)
% hObject    handle to InfoPOI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
npoi=size(handles.frq_POI,1);
tpoi=mod(handles.frq_POI(:,1) + handles.offset,24*60^2);
[spoid, spois]=crc_time_converts(tpoi);
tpoi=mod(handles.frq_POI(:,2) + handles.offset,24*60^2);
[epoid, epois]=crc_time_converts(tpoi);
durpoi=handles.frq_POI(:,2)-handles.frq_POI(:,1);
 
if handles.figz~=0
    z=handles.figz;
else
    z=figure;
end
figure(z)
close(z)
figure(z)
 
set(z,'NumberTitle','off')
set(z,'Name','Periods of Interest')
handles.figz = z;
 
axis off
set(z,'PaperOrientation','landscape');
y=0.9:-0.05:0.1;
text(0,1,'POI n°')
text(0.25,1,'Start time')
text(0.5,1,'End time')
text(0.75,1,'Duration')
if ~isempty(handles.frq_POI)
    for i=1:size(handles.frq_POI,1)
        text(0 ,y(i),num2str(i))
        text(0.25, y(i),spois(i,:))
        text(0.5, y(i),epois(i,:))
        text(0.75, y(i),num2str(durpoi(i)))
    end
    text(0,0.05,['Number of POIs: ',num2str(npoi)])
else
    text(0,1,'No Periods of Interest in this file')
end
 
filename=[pth(handles.Dmeg),filesep,'Periods_of_interest.txt'];
fid = fopen(filename,'w+');
A=[(1:npoi)',spoid,epoid,durpoi];
fprintf(fid,'POI n° %2.0f starts at %2.0f hh- %2.0f mm- %2.0f ss and stops at %2.0f hh- %2.0f mm- %2.0f ss, with a duration of %5.2f s\n',A');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function handles = mainplot_frq(handles)
 
%see in which axes to plot (if export, then use the new figure)
if handles.export
    axs=get(handles.currentfigure,'CurrentAxes');
    h=axs;
else
    set(handles.figure1,'CurrentAxes',handles.axes1)
    cla(handles.axes1)
    h=handles.axes1;
end
axes(h)
 
% subsmpl = handles.subsmpl;
toshow=1:size(handles.Dmeg.CRC.pwrspect.frqdata,3);
temps = handles.Dmeg.CRC.pwrspect.duration/2: ...
    handles.Dmeg.CRC.pwrspect.step: ...
    size(handles.Dmeg.CRC.pwrspect.frqdata,3)*handles.Dmeg.CRC.pwrspect.step;
contents = get(handles.selchanspectro,'String');
selchan=upper(contents{get(handles.selchanspectro,'Value')});
[dumb1,index]=intersect(upper(chanlabels(handles.Dmeg)),selchan);
handles.spectrochan=index;
 
 
%Plot blue line
handles.bluelines=[];
frqtoplot=find(handles.Dmeg.CRC.pwrspect.frqbins>=handles.fminblue & ...
    handles.Dmeg.CRC.pwrspect.frqbins<=handles.fmaxblue);
for i=1:length(index)
    hold on
    if isempty(frqtoplot)
        vectortoplot=zeros(1,length(temps));
    else
        % Classic case
        vectortoplot=mean(squeeze(handles.Dmeg.CRC.pwrspect.frqdata(index(i), ...
            frqtoplot,1:length(temps))));
 
        % Mongrain case
        if get(handles.Pwr_RA,'Value')==3
            scorer = get(handles.scorerpopmenu,'Value');
            SWS=find(handles.Dmeg.CRC.score{scorer}==3|handles.Dmeg.CRC.score{scorer}==4);
            timSWS=[];
            for ii=1:length(SWS)
                tmp=(SWS(ii)-1)*handles.Dmeg.CRC.pwrspect.step+1:SWS(ii)*handles.Dmeg.CRC.pwrspect.step;
                timSWS=[timSWS tmp];
            end
            try
                dv=mean(mean(squeeze(handles.Dmeg.CRC.pwrspect.frqdata(index(i),frqtoplot,timSWS))));
            catch
                beep
                disp('No Sleep Stage 3 or 4 scored: dividing by the mean over the whole recording')
                dv=mean(mean(squeeze(handles.Dmeg.CRC.pwrspect.frqdata(index(i),frqtoplot,:))));
            end
            vectortoplot=vectortoplot/dv;            
        elseif get(handles.Pwr_RA,'Value')==2         
            dv=sum(squeeze(handles.Dmeg.CRC.pwrspect.frqdata(index(i), ...
            :,1:length(temps))),1);        
            vectortoplot=vectortoplot./dv;
        else    
            vectortoplot=20*log(vectortoplot);
        end
    end
    subsmpl = handles.subsmpl;
    toadd = mod(length(vectortoplot),subsmpl);
    vectortoplot = [vectortoplot zeros(1,subsmpl-toadd)];
    vectortoplot = reshape(vectortoplot,subsmpl,length(vectortoplot)/subsmpl);
    vectortoplot = mean(vectortoplot);
    temps = temps(1:subsmpl:length(temps));
    if length(vectortoplot)>length(temps)
        tps = temps(length(temps))-temps(length(temps)-1);
        temps = [temps (temps(length(temps))+tps)];
    end
    numline=plot(temps,handles.scale*i+vectortoplot);
    handles.bluelines=[handles.bluelines numline];
    vectortoplot=[];
 
end
 
ylim([handles.scale/2 handles.scale*(i+1)])
xlim([handles.Dmeg.CRC.pwrspect.duration/2-1  ...
    nsamples(handles.Dmeg)/(fsample(handles.Dmeg)+1)])
 
%Define the label/ticks of the axes
set(handles.axes1,'YTick',[handles.scale:handles.scale/2:i*handles.scale+handles.scale/2]);
ylabels=[];
for j=1:length(index)
    q=cell(1);
    q{1}='0';
    ylabels=[ylabels q];
    ylabels=[ylabels num2str(round(handles.scale*100/2)/100)];
end
set(handles.axes1,'YTickLabel',ylabels);
 
% Red lines def
contents = get(handles.secondchan,'String');
selchan=upper(contents{get(handles.secondchan,'Value')});
[dumb1,index]=intersect(upper(chanlabels(handles.Dmeg)),selchan);
handles.spectrochan=index;
 
toshow=1:size(handles.Dmeg.CRC.pwrspect.frqdata,3);
temps = handles.Dmeg.CRC.pwrspect.duration/2: ...
    handles.Dmeg.CRC.pwrspect.step:...
    size(handles.Dmeg.CRC.pwrspect.frqdata,3)*handles.Dmeg.CRC.pwrspect.step;
set(handles.figure1,'CurrentAxes',handles.axes1)
handles.redlines=[];
frqtoplot=find(handles.Dmeg.CRC.pwrspect.frqbins>=handles.fminblue & ...
    handles.Dmeg.CRC.pwrspect.frqbins<=handles.fmaxblue);
 
cindex=index;
for i=1:length(cindex)
    hold on
    if isempty(frqtoplot)
        vectortoplot=zeros(1,length(temps));
    else
        % Classic case
        vectortoplot=mean(squeeze(handles.Dmeg.CRC.pwrspect.frqdata(index(i), ...
            frqtoplot,1:length(temps))));
 
        % Mongrain case
        if get(handles.Pwr_RA,'Value')==3
            scorer = get(handles.scorerpopmenu,'Value');
            SWS=find(handles.Dmeg.CRC.score{scorer}==3|handles.Dmeg.CRC.score{scorer}==4);
            timSWS=[];
            for ii=1:length(SWS)
                tmp=(SWS(ii)-1)*handles.Dmeg.CRC.pwrspect.step+1:SWS(ii)*handles.Dmeg.CRC.pwrspect.step;
                timSWS=[timSWS tmp];
            end
            try
                dv=mean(mean(squeeze(handles.Dmeg.CRC.pwrspect.frqdata(index(i),frqtoplot,timSWS))));
            catch
                beep
                disp('No Sleep Stage 3 or 4 scored: dividing by the mean over the whole recording')
                dv=mean(mean(squeeze(handles.Dmeg.CRC.pwrspect.frqdata(index(i),frqtoplot,:))));
            end
            vectortoplot=vectortoplot/dv;            
        elseif get(handles.Pwr_RA,'Value')==2            
            dv=sum(squeeze(handles.Dmeg.CRC.pwrspect.frqdata(index(i), ...
            :,1:length(temps))),1);        
            vectortoplot=vectortoplot./dv;
        else    
            vectortoplot=20*log(vectortoplot);
        end
    end
    subsmpl = handles.subsmpl;
    toadd = mod(length(vectortoplot),subsmpl);
    vectortoplot = [vectortoplot zeros(1,subsmpl-toadd)];
    vectortoplot = reshape(vectortoplot,subsmpl,length(vectortoplot)/subsmpl);
    vectortoplot = mean(vectortoplot);
    temps = temps(1:subsmpl:length(temps));
    if length(vectortoplot)>length(temps)
        tps = temps(length(temps))-temps(length(temps)-1);
        temps = [temps (temps(length(temps))+tps)];
    end
    numline=plot(temps,handles.scale*i+vectortoplot,'r');
    handles.redlines=[handles.redlines numline];
    vectortoplot=[];
end
 
 
%Plot Periods of Interest
handles.numpoi=[];
if ~isempty(handles.frq_POI)
    npoi=size(handles.frq_POI,1);
    fact=max(length(handles.spectrochan),1);
    for i=1:npoi
        if handles.export
            numpoi=plot(ones(1,2)*handles.frq_POI(i,1),[0 handles.scale*(fact+1)], ...
                'Color',[0 0 0]);
        else
            numpoi=plot(ones(1,2)*handles.frq_POI(i,1),[0 handles.scale*(fact+1)], ...
                'UIContextMenu',handles.Deletemenu,'Color',[0 0 0]);
        end
        text(handles.frq_POI(i,1),handles.scale*(fact+7/8), ...
             'S.POI','Color',[0 0 0],'FontSize',14)
         handles.numpoi=[handles.numpoi, numpoi];
    end
    for i=1:npoi
        if handles.export
            numpoi=plot(ones(1,2)*handles.frq_POI(i,2),[0 handles.scale*(fact+1)], ...
                'Color',[0 0 0]);
        else
           numpoi=plot(ones(1,2)*handles.frq_POI(i,2),[0 handles.scale*(fact+1)], ...
                'UIContextMenu',handles.Deletemenu,'Color',[0 0 0]);
        end
         text(handles.frq_POI(i,2),handles.scale*(fact+7/8), ...
             'E.POI','Color',[0 0 0],'FontSize',14)
         handles.numpoi=[handles.numpoi, numpoi];
    end
end
 
 
return
 
%%%%
 
function cleargraph(fig)
 
A=get(fig,'Children');
idx=find(strcmp(get(A,'Type'),'axes')==1);
try
    delete(get(A(idx),'Children'))
end