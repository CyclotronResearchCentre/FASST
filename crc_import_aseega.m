function Do = crc_import_aseega(dat, flags)
% Do = crc_import_aseega(dat, flags)
%
% Importing ASEEGA score and structure into the corresponding @meeg object.
% Specifically a few operations can be specified into the flags structure 
% to control what is included as an event.
% Otherwise the whole structure is saved into D.aseega substructure
%
% INPUT
%   dat     : filename to or assega structure itself
%   flags   : structure with options
%       .addspin : add spindles as events in the object [1, def]
%       .addab   : add alpha bursts as events in the object [1, def]
%       .save    : save or not the output updatred structure [1, def]
%
% OUTPUT
%   Do       : updated meeg object(s) (also saved on disk if requested)
%__________________________________________________________________________
% Copyright (C) 2013 Cyclotron Research Centre

% Written by C. Phillips, 2013.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

def_flags = struct(...
    'addspin', true,...
    'addab', true,...
    'save', true,...
    'pathname',false...
    );
if nargin<2
    flags = def_flags;
else
    flags = crc_check_flag(def_flags,flags);
end

if nargin<1 || isempty(dat)
    fn = spm_select(Inf,'mat','Select ASEEGA score mat file(s)', ...
        '',pwd,'.*AseegaAnalysis.A4R.mat');
    [dat,Ndat] = subfct_loaddata(fn);
else
    if ischar(dat)
        fn = dat;
        [dat,Ndat] = subfct_loaddata(fn);
    elseif isstruct(dat)
        % assume it's the aseega structure (possibly an array of it)
        fn = [];
        Ndat = numel(dat);
    else
        error('fasst:aseega','Wrong input data format')
    end
end
% check the path for the fasst files 
if numel(flags.pathname)>1
    pth = flags.pathname;
else 
    pth = spm_select(1,'dir','Select folder with files with Fasst');
end
%% Loop through all the aseega data
Do = cell(Ndat,1); % initialize output
for ii=1:Ndat
    % find the EEGdata .mat file
    fn_mati = dat(ii).info.recording.filename;
    if ~isempty(fn)
        % try adding path to the .mat file,
        ini = spm_str_manip(fullfile(pth,['M',fn_mati]),'r');
        fn_mat= [ini '.mat'];
    else
        % or assume it's in the current dir
        ini= spm_str_manip(fullfile(pwd,fn_mati),'r'); 
        fn_mat= [ini '.mat'];
    end
    %essayer d'ajouter le 'M' devant le nom pour dire que ça a été
    %rereference
    if ~exist(fn_mat,'file')
        ini = spm_str_manip(fullfile(pth,fn_mati),'r');
        fn_mat= [ini '.mat'];
    end
    if ~exist(fn_mat,'file')
        beep
        fprintf('\n\n WARNING !!!\n')
        fprintf('\tCan''t find .mat file "%s".\n',fn_mat)
        fprintf('\tSkipping it and moving on (if any).\n')
    else
        % load the meeg data
        D = crc_eeg_load(fn_mat);
        % include aseega info into object ("D.other" substructure)
        D.aseega = dat(ii);
        
        % proceed with specific events
        if flags.addspin % Add spindles
            Nsp = dat(ii).micro.spindles.number;
            ev = D.events;
            clear ev_new
            ev_new(Nsp) = struct('type','','value',[],'time',[], ...
                'duration',[],'offset',[]); %#ok<*AGROW>
            for jj=1:Nsp
                info = [
                    dat(ii).micro.spindles.power(jj) ...
                    dat(ii).micro.spindles.frequency(jj) ...
                    dat(ii).micro.spindles.frequential_instability(jj) ...
                    dat(ii).micro.spindles.frequential_purity(jj) ...
                    dat(ii).micro.spindles.temporal_instability(jj) ...
                    dat(ii).micro.spindles.maximum_amplitude(jj)];
                ev_new(jj).type     = 'ASEEGA_sp';
%                 ev_new(jj).value    = dat(ii).micro.spindles.frequency(jj);
                ev_new(jj).value    = info;
                ev_new(jj).time     = dat(ii).micro.spindles.position(jj);
                ev_new(jj).duration = dat(ii).micro.spindles.duration(jj);
                ev_new(jj).offset   = 0;
            end
            ev = [ev ev_new]; %#ok<AGROW>
            D = events(D, 1, ev);
        end
        if flags.addab % Add alpha_bursts
            Nab = dat(ii).micro.alpha_bursts.number;
            ev = D.events;
            clear ev_new
            ev_new(Nab) = struct('type','','value',[],'time',[], ...
                'duration',[],'offset',[]); %#ok<*AGROW>
            for jj=1:Nab
                info = [
                    dat(ii).micro.alpha_bursts.power(jj) ...
                    dat(ii).micro.alpha_bursts.frequency(jj) ...
                    dat(ii).micro.alpha_bursts.frequential_instability(jj) ...
                    dat(ii).micro.alpha_bursts.frequential_purity(jj) ...
                    dat(ii).micro.alpha_bursts.temporal_instability(jj) ...
                    dat(ii).micro.alpha_bursts.maximum_amplitude(jj)];
                ev_new(jj).type     = 'ASEEGA_ab';
%                 ev_new(jj).value    = dat(ii).micro.alpha_bursts.frequency(jj);
                ev_new(jj).value    = info;
                ev_new(jj).time     = dat(ii).micro.alpha_bursts.position(jj);
                ev_new(jj).duration = dat(ii).micro.alpha_bursts.duration(jj);
                ev_new(jj).offset   = 0;
            end
            ev = [ev ev_new]; %#ok<AGROW>
            D = events(D, 1, ev);
        end
        
        % Add Aseega scores as FASST scores
        sc = D.CRC.score;
        n = D.CRC.score(2,:);
        %if the scores have already imported, the news events are included
        %in the Aseega scorer previously added
        if any(strcmp(n,'Aseega'))
            isc = find(strcmp(n,'Aseega'));
        else
            isc = size(sc,2)+1;
        end
        sc_aseega = crc_convert_AseegaScore(dat(ii).scoring.hypno6);
        sc{1,isc} = sc_aseega; %#ok<*NASGU>
        sc{2,isc} = 'Aseega';
        sc{3,isc} = 30;
        sc{4,isc} = sc{4,1};             
        D.CRC.score = sc;
        
        Do{ii} = D; %#ok<AGROW>     
        % saving, only if requested
        if flags.save
            save(D);
        end
    end
end
end

%% SUBFUNCTIONS
function [dat,Ndat] = subfct_loaddata(fn)

Ndat = size(fn,1);
for ii=1:Ndat
    tmp = load(deblank(fn(ii,:)));
    dat(ii) = tmp.aseega_analysis;
end

end
