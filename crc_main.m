function varargout = crc_main(varargin)
%__________________________________________________________________________
%     ________  _____________ ___
%    / __ _/  |/ ___/_/_  __/|__ \
%   / /_ / /| |\__ \ \ / /   __/ /
%  / __// ___ |__/ / // / _ / __/
% /_/  /_/  |_|___/_//_/ (_)____/
%
% fMRI Artefact removal and Sleep Scoring Toolbox, FASST.2
% http://www.montefiore.ulg.ac.be/~phillips/FASST.html
%__________________________________________________________________________
%
% CRC_MAIN M-file for crc_main.fig
%      CRC_MAIN, by itself, creates a new CRC_MAIN or raises the existing
%      singleton.
%
%      H = CRC_MAIN returns the handle to a new CRC_MAIN or the handle to
%      the existing singleton*.
%
%      CRC_MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CRC_MAIN.M with the given input arguments.
%
%      CRC_MAIN('Property','Value',...) creates a new CRC_MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before crc_main_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to crc_main_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help crc_main

% Last Modified by GUIDE v2.5 16-Jun-2018 15:00:03

% Display ASCII Welcome
try
    version = load('z3score-version.txt');
    fprintf(' FASST-Z3Score version no: %d.%02d.%02d \n',version);
catch
    
end
disp('                                                            ');
disp('     ________  _____________ ___                            ');
disp('    / ____/  |/ ___/_/_  __/|__ \                           ');
disp('   / /_ / /| |\__ \ \ / /   __/ /                           ');
disp('  / __// ___ |__/ / // / _ / __/                            ');
disp(' /_/  /_/  |_|___/_//_/ (_)____/                            ');
disp('                                                            ');
disp(' fMRI Artefact removal and Sleep Scoring Toolbox, FASST.2   ');
disp(' http://www.montefiore.ulg.ac.be/~phillips/FASST.html       ');
disp(' An SPM8/12-compatible toolbox.                             ');
disp('															  ');
disp('															  ');
disp('         ****                                               ');
disp('        */// *                                              ');
disp(' ******/    /*  ******  *****   ******  ******  *****       ');
disp('////**    ***  **////  **///** **////**//**//* **///**      ');
disp('   **    /// *//***** /**  // /**   /** /** / /*******      ');
disp('  **    *   /* /////**/**   **/**   /** /**   /**////       ');
disp(' ******/ ****  ****** //***** //****** /***   //******      ');
disp('//////  ////  //////   /////   //////  ///     //////       ');
disp('															  ');
disp(' Z3Score.com - Automatic Sleep Scoring & Artifact Detection ');
disp(' Z3Score is powered by Neurobit Technologies Pte. Ltd.      ');
disp(' More details at: https://z3score.com          			  ');
disp(' For API documentation and sample code visit:               ');
disp(' https://github.com/neurobittechnologies/z3score-api        ');
disp(' this version now uses Z3Score-V2 (NEO) sleep staging system.');
disp(' Z3Score V3 shows 35 to 40% lower error rates compared to V1.');
disp('															  ');
disp('															  ');
disp(' Cite this version as:                                       ');
disp(' Z3Score-Cloud Version 3.0, Neurobit Technologies, Singapore.');
disp(' and Patanaik et al., Sleep 41(5), 2018                       ');


%% Checking the installation and defaults
persistent flag_filter flag_initialize

if isempty(flag_initialize)
    % Check if SPM is available
    [ok,spm_v] = check_installation;
    if ~ok
        beep
        fprintf('INSTALLATION PROBLEM!');
        return
    end
    
    % Add the fieldtrip toolbox from SPM, if necessary
    if ~exist('ft_defaults','file')
        addpath(fullfile(spm('Dir'),'external','fieldtrip'));
    end
    ft_defaults;
    flag_initialize = true;
end

% Check for Signal Processing Toolbox & SPM12 fixes
if isempty(flag_filter)
    flag_TBX = license('checkout','signal_toolbox');
    if ~flag_TBX && spm_v ~= 12
        pth = fullfile(spm_str_manip(mfilename('fullpath'),'h'),'SPTfunctions');
        addpath(pth)
        disp(['warning: using freely distributed equivalent to filtering functions ', ...
            'as Signal Processing Toolbox is not available.']);
    end
    if spm_v == 12 && flag_TBX
        disp(['warning: using FiledTrip distribution (SPM12) of filtering functions ', ...
            'instead of Signal Processing Toolbox.']);
        
    end
    if spm_v == 12
        pth = fullfile(spm_str_manip(mfilename('fullpath'),'h'),'spm12_utils');
        addpath(pth)
    end
    flag_filter = true;
end

if ~nargin || (ischar(varargin{1}) && ~strcmp(varargin{1},'SetDefs'))
    %% Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
        'gui_Singleton',  gui_Singleton, ...
        'gui_OpeningFcn', @crc_main_OpeningFcn, ...
        'gui_OutputFcn',  @crc_main_OutputFcn, ...
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
else
    varargout{1} = true; % Defaults were set
end

% --- Executes just before crc_main is made visible.
function crc_main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to crc_main (see VARARGIN)


[A] = imread('LOGO_Simple.png','BackgroundColor',0.94*[1 1 1]);
image(A)
axis off

% Choose default command line output for crc_main
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes crc_main wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = crc_main_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in push_dis_main.
function push_dis_main_Callback(hObject, eventdata, handles)
% hObject    handle to push_dis_main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%delete(handles.figure1)
dis_selchan;

% --- Executes on button press in push_dis_cmp.
function push_dis_cmp_Callback(hObject, eventdata, handles)
% hObject    handle to push_dis_cmp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% flags.multcomp=1;
% setappdata(hObject,'multcomp',1);
flags.multcomp=1;
dis_selchan(flags);

% --- Executes on button press in push_credits.
function push_credits_Callback(hObject, eventdata, handles)
% hObject    handle to push_credits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fig = figure;
set(fig,'Position',[73   145   398   503])
set(fig,'NumberTitle','off')
set(fig,'Name','License & Copyright')
h = axes('Position',[0 0 1 1],'Visible','off');

str{1}= 'Please refer to this version as "FASST, fMRI Artefact Removal and';
str{2}= 'Sleep scoring Toolbox" in papers and communication';
str{3}= ' ';
str{4}= '_____________________________________________________________________';
str{5}= ' ';
str{6}= '';
str{7}= 'For bug reports, contact directly';
str{8}= '';
str{9}= 'Jessica Schrouff, jschrouff@doct.ulg.ac.be';
str{10}='Dorothée Coppieters, d.coppieters@ulg.ac.be';
str{11}='Christophe Phillips, c.phillips@ulg.ac.be';
str{12}=' ';
str{13}='Feel free to submit add-ons to this toolbox';
str{14}='but please note that support will not be provided for';
str{15}='"home-made" changes of the distributed code.';
str{16}=' ';
str{17}='More details are available on';
str{18}='             http://www.montefiore.ulg.ac.be/~phillips/FASST.html';
str{19}=' ';
str{20}='Reference:';
str{21}='   Y. Leclercq, J. Schrouff, Q. Noirhomme, P. Maquet, C. Phillips,';
str{22}='   "fMRI Artefact rejection and Sleep Scoring Toolbox",';
str{23}='   Computational Intelligence and Neuroscience, 2011,';
str{24}='   http://www.hindawi.com/journals/cin/2011/598206/';
str{25}='___________________________________________________________________';
str{26}=' ';
str{27}='FASST is developed by the Cyclotron Research Centre, part of';
str{28}='the University of Liege (ULg), BE. This work is supported by the ';
str{29}='FRS-FNRS, the Queen Elizabeth''s funding and the University of Liege.';
str{30}=' ';
str{31}='___________________________________________________________________';
str{32}=' ';
str{33}='FASST (being the collection of files given in the manifest in the';
str{34}='Contents.m file) is free but copyright software, distributed under the';
str{35}='terms of the GNU General Public Licence as published by the Free Software';
str{36}='Fundation (either version 2, as given in file CRC_LICENCE.man or at your';
str{37}='option, any later version). Further details on "copyleft" can be found at';
str{38}='http://www.gnu.org/copyleft/';
str{39}=' ';
str{40}='___________________________________________________________________';
str{41}='Copyright (C) 2010 Cyclotron Research Centre, University of Liege';
str{42}=' ';
text(.025,.5,str,'FontSize',8)

% --- Executes on button press in push_concatenate.
function push_concatenate_Callback(hObject, eventdata, handles)
% hObject    handle to push_concatenate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
crc_concatenate;

% --- Executes on button press in push_disfrqcomp.
function push_disfrqcomp_Callback(hObject, eventdata, handles)
% hObject    handle to push_disfrqcomp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dis_frqcomp;

% --- Executes on button press in push_freqplot.
function push_freqplot_Callback(hObject, eventdata, handles)
% hObject    handle to push_freqplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dis_frq;

% --- Executes on button press in push_crcgar.
function push_crcgar_Callback(hObject, eventdata, handles)
% hObject    handle to push_crcgar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
crc_gar;

% --- Executes on button press in push_crcpar.
function push_crcpar_Callback(hObject, eventdata, handles)
% hObject    handle to push_crcpar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
crc_par;

% --- Executes on button press in push_score.
function push_score_Callback(hObject, eventdata, handles)
% hObject    handle to push_score (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dis_selchan;

% --- Executes on button press in push_freqplotstat.
function push_freqplotstat_Callback(hObject, eventdata, handles)
% hObject    handle to push_freqplotstat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dis_frq;

% --- Executes on button press in push_chunk.
function push_chunk_Callback(hObject, eventdata, handles)
% hObject    handle to push_chunk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
crc_chunks;

% --- Executes on button press in sws.
function sws_Callback(hObject, eventdata, handles)
% hObject    handle to sws (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
crc_wave_detection;

function push_batch_Callback(hObject, eventdata, handles)
% hObject    handle to push_chunk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
crc_batch;

%% SUBFUNCTION

function [ok,spm_v] = check_installation
% function to check installation state of toolbox,
% particullarly the SPM path setup
%
% OUTPUT
% ok    : logical flag
% spm_v : spm version number

ok = true;

% Check SPM installation
if exist('spm.m','file')
    [SPMver, SPMrel] = spm('Ver');
    if ~(strcmpi(SPMver,'spm8') && str2double(SPMrel)>8.5) && ...
            ~strcmpi(SPMver,'spm12')
        beep
        fprintf('\nERROR:\n')
        fprintf('\tThe *latest* version of SPM8 or SPM12 should be installed on your computer,\n')
        fprintf('\tand be available on MATLABPATH!\n\n')
        ok = false;
    end
    spm_v = str2double(SPMver(4:end));
else
    beep
    fprintf('\nERROR:\n')
    fprintf('\tThe *latest* version of SPM8/12 should be installed on your computer,\n')
    fprintf('\tand be available on MATLABPATH!\n\n')
    ok = false;
end

return
