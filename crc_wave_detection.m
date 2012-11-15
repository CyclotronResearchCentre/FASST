function varargout = crc_wave_detection(varargin)
% CRC_WAVE_DETECTION M-file for crc_wave_detection.fig
%      CRC_WAVE_DETECTION, by itself, creates a new CRC_WAVE_DETECTION or raises the existing
%      singleton*.
%
%      H = CRC_WAVE_DETECTION returns the handle to a new CRC_WAVE_DETECTION or the handle to
%      the existing singleton*.
%
%      CRC_WAVE_DETECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CRC_WAVE_DETECTION.M with the given input arguments.
%
%      CRC_WAVE_DETECTION('Property','Value',...) creates a new CRC_WAVE_DETECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before crc_wave_detection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to crc_wave_detection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%__________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by J. Schrouff & C. Phillips, 2009.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id:$

% Edit the above text to modify the response to help crc_wave_detection

% Last Modified by GUIDE v2.5 21-Sep-2010 15:10:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @crc_wave_detection_OpeningFcn, ...
                   'gui_OutputFcn',  @crc_wave_detection_OutputFcn, ...
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


% --- Executes just before crc_wave_detection is made visible.
function crc_wave_detection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to crc_wave_detection (see VARARGIN)

% Choose default command line output for crc_wave_detection
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes crc_wave_detection wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = crc_wave_detection_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pb_spindles.
function pb_spindles_Callback(hObject, eventdata, handles)
% hObject    handle to pb_spindles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
crc_SP_detection

% --- Executes on button press in pb_sws.
function pb_sws_Callback(hObject, eventdata, handles)
% hObject    handle to pb_sws (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
crc_sws_detection

