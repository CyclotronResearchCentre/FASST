function newD = crc_interpol(Dbadc,config)

% csg_interpol interpolate bad channels detected from the csg_badchannels detection method (for each 30s window).
% if varargin is not empty, it is a struct containing:
% ° .filename*   : filename of a meeg file
% ° .method      : 'spline' or 'average' interpolation method
% ° .lambda      : smoothing (?)
% ° .order       : m-constant for the spline interpolation
%  This code needs functions fromt the fieldtrip toolbox (i.e. ft_channelrepair) available in SPM12

% parameters from the data file loaded
fs      =   fsample(Dbadc);         % sampling frequency
Lw      =   config.winsize;        % fix time window used for bad channels detection
nspl    =   nsamples(Dbadc);        % number of samples
NofW    =   ceil(nspl/(fs*Lw));     % number of 15s window
channels    =   config.channels;
newdata.label   =   chanlabels(Dbadc,channels);

for ibc = 1 : NofW
    badchannels{ibc} = union(Dbadc.CRC.DC.badchannels.chan_defaillant{ibc},Dbadc.CRC.DC.badchannels.chan_incoherent{ibc});
end

data        =   spm2fieldtrip(Dbadc); % transform meeg object from spm into a raw data of fieldtrip
newdata     =   data;             % new raw data with interpolated values 
filename    =   fnamedat(Dbadc);
[pth, fname]    =   fileparts(filename);
newname     =   ['I' fname];
newfilename =   fullfile(pth,newname);

% clone the data 
newD    =   clone(Dbadc,newfilename);
coord   =   coor2D(Dbadc)';

% interpolation: configuration (to be checked, still in progess)
elec.chanpos    =   [coord(channels,:) zeros(numel(channels),1)];
elec.elecpos    =   [coord(channels,:) zeros(numel(channels),1)];
elec.label      =   chanlabels(Dbadc,channels); 
cfglayout.elec  =	elec;
cfglayout.layout	= 'ordered';
lay     =   ft_prepare_layout(cfglayout, data);

cfgneigh.elec   =   elec;
cfgneigh.method	=   'triangulation';% or 'template' (default = 'distance')
cfgneigh.layout	=   lay;
cfg.neighbours  =   ft_prepare_neighbours(cfgneigh, data);

cfg.elec    =   elec; % channels position and labels given from the elec_EGI256.mat file according to the fieldtrip toolbox
cfg.method  =   ft_getopt(cfg, 'method','spline');
cfg.lambda  =   ft_getopt(cfg, 'lambda',[]); % subfunction will handle this
cfg.order   =   ft_getopt(cfg, 'order',[]); % subfunction will handle this
newD.CRC.interpolation.cfg  =   cfg;
newD(:,:,1) =   Dbadc(:,:,1);

% loop to search bad channels in each 30s-window and interpolate bad
% channels from the method chosen
for iw = 1 : NofW
    fprintf(1,' \n');
    tempo   =   max(1,(iw-1)*Lw*fs):min(nspl,iw*fs*Lw);
    if ~isempty(badchannels{iw})
        cfg.badchannel      =   chanlabels(Dbadc,badchannels{iw});  % bad channel to interpolate
        newdata.trial{1}    =   Dbadc(channels,tempo);
        newdata.time{1}     =   data.time{1}(:,tempo);
        inter               =   ft_channelrepair(cfg,newdata);
        newD(channels,tempo)       =   inter.trial{1};
    end
end
fprintf(1,' -------------- Interpolation processed ----------------- \n ');
save(newD);

