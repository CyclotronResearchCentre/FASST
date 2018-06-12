function crc_ASDM(config,iid)

% load data 
try 
    D = spm_eeg_load(char(deblank(config.filename(iid,:))));
catch
    error('unable to load ',deblank(config.filename(iid,:)))
end

% ----- load parameters values -------
name = fname(D);
fs = fsample(D);
nspl = nsamples(D);
Time = floor(nspl/fs);
tr_tsup = 5;
pathname = path(D);
filename = fname(D);

% electrods selected 
if ischar(cfg.channels) && strcmp(cfg.channels,'ALL')
    channels = meegchannels(D,'EEG');
elseif ischar(cfg.channels)
    channels = find(strcmp(chantype(D),cfg.channels));
    if isempty(channels) || ~isempty(intersect(eogchannels(D),cfg.channels)) || ~isempty(intersect(emgchannels(D),cfg.channels)) || ~isempty(intersect(ecgchannels(D),cfg.channels))
        error('Error with channels selected')
    end
else
    % check indexes
    if max(cfg.channels)>size(D,1)
        error('At least one index exceeds the number of channels available')
    elseif ~isempty(intersect(eogchannels(D),cfg.channels)) || ~isempty(intersect(emgchannels(D),cfg.channels)) || ~isempty(intersect(ecgchannels(D),cfg.channels))
        error('Error with channels selected')
    end  
    
end 
 
ne = numel(channels);
uni = units(D,channels);
scal = zeros(1,size(D,1));
for iu = 1 : length(channels)
    if strcmp(uni(iu),'T/m')
        scal(channels(iu)) = 5e-13;
    elseif strcmp(uni(iu),'T')
        scal(channels(iu)) = 5e-14;
    elseif strcmpi(uni(iu),'V')
        scal(channels(iu)) = 1e-6;
    elseif strcmpi(uni(iu),'UV')
        scal(channels(iu)) = 1;
    else 
        scal(channels(iu)) = 1;
    end
end 

% I ran it already (FR)
% % make the preprocessing including filtering and artefact and arousal
% % detection
% DC_artifact(D);

% initialization of structure for sleep spindle characteristics
spindle_features = struct('time',zeros(ne,Time,2),'type',zeros(ne,Time,2),'frequency',struct('mode_mean',nan(ne,Time),'mode_arma',nan(ne,Time),'intrafreq',nan(ne,Time,tr_tsup*fs)),'shape',struct('skewness',nan(ne,Time),'kurtosis',nan(ne,Time),'magnitude',nan(ne,Time,2),'energy',nan(ne,Time)),'power',struct('max',nan(ne,Time,1),'intrapower',nan(ne,Time,tr_tsup*fs)),'chan',nan(Time,ne));
spindle_features.chan = chanlabels(D,channels);

% NREM PS computation for each channel selected at a time 
addpath 'C:\Users\Franziska\Documents\FRANZISKA\SiRENE\FHR\Cooperation\Liege_Toolbox\ToolboxES\AASM_toolbox\ASDM'
NREMcell = crc_NREMclassifier(D,channels); 
% % % % 
for ie = 1 : numel(channels) 
    
    nrem = NREMcell{ie};
    % estimate the sigma energy to detect sigma burst
    % data    =   spm2fieldtrip(D);
    data = D.ftraw;
    cfg.lpfilter    =   'yes';
    cfg.hpfilter    =   'yes';
    cfg.lpfreq      =   16;    % Low pass filter
    cfg.hpfreq      =   11;    % High pass filter
    cfg.datafile    =   fnamedat(D);
    cfg.channel     =   chanlabels(D,ie);
    newfile     =   ft_preprocessing(cfg,data);
    sigmact     =   newfile.trial{1};
    [burst teo] =   crc_energydecomposition(sigmact,nrem,fs,cfg); 
    for ib = 1 : size(burst,1)
         st     =   burst(ib,1)-fs;
         et     =   burst(ib,2)+fs;
         event  =   D(ie,st:et);  
         % TF and AR decomposition => + check of the consistency between
         % both frequency computations
         [flags_sst tempo fmode_mean fmode_arma intrafreq Pmax intrapow] = crc_TFdecomposition(event,fs,teo(st:et));
         for ifl = 1 : flags_sst
             indx  =  find(~spindle_features.time(ie,:,1),1);
             timin =  st+tempo(ifl,:);
             spindle_features.type(ie,indx,:) =  [cfg.hpfreq  cfg.lpfreq ];
             spindle_features.time(ie,indx,:) =  timin./fs;
             spindle_features.frequency.mode_mean(ie,indx) =  fmode_mean(ifl);
             spindle_features.frequency.mode_arma(ie,indx) =  fmode_arma(ifl);
             spindle_features.frequency.intrafreq(ie,indx,1:numel(intrafreq{ifl})) = intrafreq{ifl};
             spindle_features.power.max(ie,indx)  =  Pmax(ifl);
             spindle_features.power.intrapower(ie,indx,1:numel(intrafreq{ifl}))  =  intrapow{ifl};
             spindle_features.shape.skewness(ie,indx) = skewness(teo(timin(1):timin(2)));                 
             spindle_features.shape.kurtosis(ie,indx) = kurtosis(teo(timin(1):timin(2)));
             [magnitude imx] = max(teo(timin(1):timin(2)));
             imx  =  (imx + timin(1))/fs;
             spindle_features.shape.magnitude(ie,indx,:)  =  [magnitude imx];
             spindle_features.shape.energy(ie,indx)       =  sum(teo(timin(1):timin(2)));
         end
    end
end

 indx_sst = zeros(size(channels));
 for ie = 1 : ne
     indx_sst(ie) = find(~spindle_features.time(ie,:,1),1);
 end
 nmaxspi = max(indx_sst);
 spindle_features.type  =  spindle_features.type(:,1:nmaxspi,:);
 spindle_features.time  =  spindle_features.time(:,1:nmaxspi,:);
 spindle_features.frequency.mode_mean  =  spindle_features.frequency.mode_mean(:,1:nmaxspi);
 spindle_features.frequency.mode_arma  =  spindle_features.frequency.mode_arma(:,1:nmaxspi);
 spindle_features.power.max  =  spindle_features.power.max(:,1:nmaxspi);
 spindle_features.shape.skewness  =  spindle_features.shape.skewness(:,1:nmaxspi);                 
 spindle_features.shape.kurtosis  =  spindle_features.shape.kurtosis(:,1:nmaxspi);
 spindle_features.shape.magnitude =  spindle_features.shape.magnitude(:,1:nmaxspi,:);
 spindle_features.shape.energy    =  spindle_features.shape.energy(:,1:nmaxspi);
 spindle_features.frequency.intrafreq  =  spindle_features.frequency.intrafreq(ie,1:nmaxspi,:);
 spindle_features.power.intrapower = spindle_features.power.intrapower(ie,1:nmaxspi,:);
% % % % 
save(fullfile(pathname,['spindle_features_' filename]),'spindle_features','-v7.3')
% Save spindles
D.CRC.DC.spindle = spindle_features;