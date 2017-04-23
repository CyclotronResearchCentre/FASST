function varargout = crc_z3score_settings(varargin)
% CRC_Z3SCORE_SETTINGS MATLAB code for crc_z3score_settings.fig
%      CRC_Z3SCORE_SETTINGS, by itself, creates a new CRC_Z3SCORE_SETTINGS or raises the existing
%      singleton*.
%
%      H = CRC_Z3SCORE_SETTINGS returns the handle to a new CRC_Z3SCORE_SETTINGS or the handle to
%      the existing singleton*.
%
%      CRC_Z3SCORE_SETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CRC_Z3SCORE_SETTINGS.M with the given input arguments.
%
%      CRC_Z3SCORE_SETTINGS('Property','Value',...) creates a new CRC_Z3SCORE_SETTINGS or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before crc_z3score_settings_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to crc_z3score_settings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help crc_z3score_settings

% Last Modified by GUIDE v2.5 14-Jan-2017 18:25:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @crc_z3score_settings_OpeningFcn, ...
                   'gui_OutputFcn',  @crc_z3score_settings_OutputFcn, ...
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

% --- Executes just before crc_z3score_settings is made visible.
function crc_z3score_settings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to crc_z3score_settings (see VARARGIN)

% Choose default command line output for crc_z3score_settings
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

initialize_gui(hObject, handles, false);

% UIWAIT makes crc_z3score_settings wait for user response (see UIRESUME)
% uiwait(handles.settings_gui);


% --- Outputs from this function are returned to the command line.
function varargout = crc_z3score_settings_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function initialize_gui(fig_handle, handles, isreset)
% If the metricdata field is present and the reset flag is false, it means
% we are we are just re-initializing a GUI by calling it from the cmd line
% while it is up. So, bail out as we dont want to reset the data.
if isfield(handles, 'metricdata') && ~isreset
    return;
end

settings_path = fullfile(getuserdir,'/z3license.mat');

if exist(settings_path, 'file') == 2,
    load(settings_path);
    server = settings.serverURL;
    email = settings.email;
    key = settings.key;
else
    server = 'http://z3score.com/api/v1';
    email = 'name@domain.com';
    key = 'yourAPIKey';
end

set(handles.server, 'String', server);
set(handles.email,  'String', email);
set(handles.key,  'String', key);


% Update handles structure
guidata(handles.settings_gui, handles);



function setLicense_Callback(hObject, eventdata, handles)
% hObject    handle to server (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of server as text
%        str2double(get(hObject,'String')) returns contents of server as a double
%check license
serverURL = get(handles.server,'String');
email = get(handles.email,'String');
key = get(handles.key,'String');

settings_path = fullfile(getuserdir,'/z3license.mat');

h = waitbar(0,'Please wait, communicating with server...');

try
    response = loadjson(urlreadpost([serverURL '/check'],...
                                        {'email',email,'key',key}));
catch
    close(h);
    errordlg('Error connecting server. Please check server address and internet connection.','Error Connecting Server');
    return
end

if response.status == 0,
    close(h);
    errordlg([' Error message: ' response.message],'License check failed.');
    return
end

close(h);
settings.serverURL = serverURL;
settings.email = email;
settings.key = key;
settings.call_limit = response.call_limit;
settings.epoch_limit = response.epoch_limit;
r = regexp(response.message,'\d*-\w*-\d*','match');
%settings.validity = datetime(r{1},'InputFormat','dd-MMMM-yyyy','TimeZone','UTC');
settings.validity = r{1};

save(settings_path, 'settings');

h = msgbox([response.message ' Call limit (hourly): ' num2str(response.call_limit) ' Epoch limit (daily): ' num2str(response.epoch_limit) '.'],'License Validated');
waitfor(h);
close(handles.settings_gui);
