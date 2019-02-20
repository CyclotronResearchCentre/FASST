function varargout = crc_z3score_genreport(varargin)
% CRC_Z3SCORE_GENREPORT MATLAB code for crc_z3score_genreport.fig
%      CRC_Z3SCORE_GENREPORT, by itself, creates a new CRC_Z3SCORE_GENREPORT or raises the existing
%      singleton*.
%
%      H = CRC_Z3SCORE_GENREPORT returns the handle to a new CRC_Z3SCORE_GENREPORT or the handle to
%      the existing singleton*.
%
%      CRC_Z3SCORE_GENREPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CRC_Z3SCORE_GENREPORT.M with the given input arguments.
%
%      CRC_Z3SCORE_GENREPORT('Property','Value',...) creates a new CRC_Z3SCORE_GENREPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before crc_z3score_genreport_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to crc_z3score_genreport_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help crc_z3score_genreport

% Last Modified by GUIDE v2.5 19-Jun-2018 02:11:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @crc_z3score_genreport_OpeningFcn, ...
                   'gui_OutputFcn',  @crc_z3score_genreport_OutputFcn, ...
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


% --- Executes just before crc_z3score_genreport is made visible.
function crc_z3score_genreport_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to crc_z3score_genreport (see VARARGIN)

% Choose default command line output for crc_z3score_genreport
handles.output = hObject;
handles.Dmeg = varargin{1};
handles.flags = varargin{2};
handles.score = handles.Dmeg.CRC.score;
current = handles.flags.currentscore;
handles.currentscore = current;
handles.scoring = 1;
handles.offset  =   0;

scores = handles.score{1,handles.currentscore};
[~,name,ext] = fileparts(fullfile(handles.Dmeg));
start = datetime([handles.Dmeg.info.date handles.Dmeg.info.hour]);
stop = datetime([handles.Dmeg.info.date handles.Dmeg.info.hour]) + seconds(nsamples(handles.Dmeg)/fsample(handles.Dmeg));
lightsoff = datetime(handles.Dmeg.lightsOff);
lightson = datetime(handles.Dmeg.lightsOn);

%everything outside lights off and on are counted as wake
epoch_size = handles.score{3,handles.currentscore};
off_epoch = ceil(seconds(lightsoff-start)/epoch_size);
on_epoch = min([ceil(seconds(lightson-start)/epoch_size), numel(scores)]);
scores(1:off_epoch-1) = 0;
scores(on_epoch:end) = 0;

% first onset of sleep
onset = find(scores~= 0 & scores <= 5,1,'first');
filename = [name ext];
scoredby = handles.score{2,handles.currentscore};

TIB = (on_epoch - off_epoch + 1)*epoch_size/60; 
latency = (onset - off_epoch)*epoch_size/60; 
REMonset = (find(scores== 5,1,'first') - off_epoch)*epoch_size/60;
REMtime = sum(scores == 5)*epoch_size/60;
NREMtime = sum(scores == 1 | scores == 2 | scores == 3)*epoch_size/60;
S1time = sum(scores == 1)*epoch_size/60;
S2time = sum(scores == 2)*epoch_size/60;
S3time = sum(scores == 3)*epoch_size/60;
TST = REMtime + NREMtime;
efficiency = TST*100/TIB;
s1percent = S1time*100/TST;
s2percent = S2time*100/TST;
s3percent = S3time*100/TST;
REMpercent = REMtime*100/TST;
NREMpercent = NREMtime*100/TST;
scr = scores(onset:on_epoch); %sleep onset to wake
WASO = sum(scr==0)*epoch_size/60;
no_of_awakenings = sum(diff(scr == 0)==1);

summary{1,1} = 'File name:'; summary{1,2} = filename;
summary{2,1} = 'Scored by:'; summary{2,2} = scoredby;
summary{3,1} = 'Epoch size (sec):'; summary{3,2} = epoch_size;
summary{4,1} = 'Recording started at:'; summary{4,2} = datestr(start);
summary{5,1} = 'Lights off time:'; summary{5,2} = datestr(lightsoff);
summary{6,1} = 'Lights on time:'; summary{6,2} = datestr(lightson);
summary{7,1} = 'Recording stopped at:'; summary{7,2} = datestr(stop);
summary{8,1} = 'Time in bed (TIB) mins:'; summary{8,2} = sprintf('%.2f',TIB);
summary{9,1} = 'Total sleep time (TST) mins:'; summary{9,2} = sprintf('%.2f',TST);
summary{10,1} = 'Wake after sleep onset (WASO) mins:'; summary{10,2} = sprintf('%.2f',WASO); 
summary{11,1} = 'Sleep onset latency (mins):'; summary{11,2} = sprintf('%.2f',latency);
summary{12,1} = 'REM latency (mins):'; summary{12,2} = sprintf('%.2f',REMonset);
summary{13,1} = 'Number of awakenings:'; summary{13,2} = no_of_awakenings;
summary{14,1} = 'Efficiency (%):'; summary{14,2} = sprintf('%.2f',efficiency);
summary{15,1} = 'Stage N1 sleep (mins):'; summary{15,2} = sprintf('%.2f',S1time);
summary{16,1} = 'Stage N2 sleep (mins):'; summary{16,2} = sprintf('%.2f',S2time);
summary{17,1} = 'Stage N3 sleep (mins):'; summary{17,2} = sprintf('%.2f',S3time);
summary{18,1} = 'Total NREM sleep (mins):'; summary{18,2} = sprintf('%.2f',NREMtime);
summary{19,1} = 'Stage REM sleep (mins):'; summary{19,2} = sprintf('%.2f',REMtime);
summary{20,1} = 'Stage N1 sleep (%TST):'; summary{20,2} = sprintf('%.2f',s1percent);
summary{21,1} = 'Stage N2 sleep (%TST):'; summary{21,2} = sprintf('%.2f',s2percent);
summary{22,1} = 'Stage N3 sleep (%TST):'; summary{22,2} = sprintf('%.2f',s3percent);
summary{23,1} = 'Total NREM sleep (%TST):'; summary{23,2} = sprintf('%.2f',NREMpercent);
summary{24,1} = 'Stage REM sleep (%TST):'; summary{24,2} = sprintf('%.2f',REMpercent);

set(handles.summarytable, 'data', summary);

set(handles.figure1,'CurrentAxes',handles.hypno)
crc_hypnoplot(handles.hypno, handles,handles.score{3,handles.currentscore},scores,scoredby);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes crc_z3score_genreport wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = crc_z3score_genreport_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on button press in savecsv.
function savecsv_Callback(hObject, eventdata, handles)
% hObject    handle to savecsv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[~,name,~] = fileparts(fullfile(handles.Dmeg));
[file, path, indx] = uiputfile([name '_' handles.Dmeg.CRC.score{2,handles.currentscore} '_sleepsummary.csv']);
if indx,
    s=get(handles.summarytable,'data');
    cell2csv([path file], s);
end


% --- Executes on button press in copytable.
function copytable_Callback(hObject, eventdata, handles)
% hObject    handle to copytable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
s=get(handles.summarytable,'data');
mat2clip(s);
msgbox('Table copied, use CTRL+V to paste.','Success');
