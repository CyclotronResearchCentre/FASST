function out = crc_run_chunking(varargin)
% FASST job execution function
% takes a harvested job data structure and call FASST functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%_______________________________________________________________________
% Copyright (C) 2010 Cyclotron Research Centre

% Christophe Phillips
% $Id$

job = varargin{1};

Dmeg = crc_eeg_load(job.data{1}); % Load data
Nchunks = numel(job.chunk); % number of chunks to create
% % identify the index for chk#_fname
% numchk = 1;
% while exist(fullfile(path(Dmeg),['chk' num2str(numchk) '_' fname(Dmeg)]),'file')
%     numchk=numchk+1;
% end

for ii=1:Nchunks
    % Beginning
    %==========
    field_n = fieldnames(job.chunk(ii).chunk_beg);
    switch field_n{1}
        % Time defined as clock time
        %---------------------------
        case 't_abs'
            % Is clock time available ???
            if ~isfield(Dmeg,'info') || isempty(Dmeg.info)
                error('Missing actual recording time of data.')
            end
            t_beg =  job.chunk(ii).chunk_beg.t_abs;
            % Is time properly defined ?
            if any(~isfinite(t_beg)) || any(t_beg<0) || any(t_beg>=[24 60 60])
                error('The parameters to determine the beginning of the chunk are not correct.')
            end

            ch_sec_rounding = round(Dmeg.info.hour(3)) - ...
                            Dmeg.info.hour(3);
            if abs(ch_sec_rounding)<1e-4
                % Rounding off if "seconds" saved differs from rounded
                % "seconds" by less than .1ms.
                Dmeg.info.hour(3) = round(Dmeg.info.hour(3));
            end
            
            datenum_start=datenum(round([0 0 0 Dmeg.info.hour]));
            datenum_stch=datenum([0 0 0 t_beg]);
            if datenum_stch<datenum_start % assume no recording is >24hours !!!
                datenum_stch=datenum([0 0 1 t_beg]);
            end
            datediff = datevec(datenum_stch-datenum_start); %#ok<NASGU>
            Begpts = max(1,round(...
                (datediff(4)*60^2 + datediff(5)*60 + datediff(6))*fsample(Dmeg)));

            if Begpts>=nsamples(Dmeg)
                error('Beginning time is over the end of recording!');
            end

        % Time defined relative to beginning of data
        %-------------------------------------------
        case 't_rel'
            t_beg =  job.chunk(ii).chunk_beg.t_rel;
            % Is time properly defined ? Assuming no data is >24h
            if any(~isfinite(t_beg)) || any(t_beg<0) || any(t_beg>=[24 60 60])
                error('The parameters to determine the beginning of the chunk are not correct.')
            end
            Begpts = max(1,round( ...
                (t_beg(1)*60^2+t_beg(2)*60+t_beg(3)) * fsample(Dmeg)));
            if Begpts>=nsamples(Dmeg)
                error('Beginning time is over the end of recording!');
            end

        % Time defined by an event.
        %--------------------------
        case 't_mark'
            % This bit should be seriously checked!
            % Any one using this way of chunking could report to me please?
            t_mark = job.chunk(ii).chunk_beg.t_mark;
            evts = events(Dmeg);
            marker_list = [evts(:).value];
            marker_match = find(marker_list==t_mark.m_type);
            % check if user's selected marker is ok
            if isempty(marker_match) % No match
                error('No matching event/marker value.')
            elseif length(marker_match)<t_mark.m_ind % Index too big
                error('Selected index beyong number of such marker/event.')
            else % ok, so pick time of event
                Begpts = max(1,round( ...
                    evts(marker_match(t_mark.m_ind)).time*fsample(Dmeg)));
            end
    end

    % End
    %====
    field_n = fieldnames(job.chunk(ii).chunk_end);
    switch field_n{1}
        % Time defined as clock time
        %---------------------------
        case 't_abs'
            % Is clock time available ???
            if ~isfield(Dmeg,'info') || isempty(Dmeg.info)
                error('Missing actual recording time of data.')
            end
            t_end =  job.chunk(ii).chunk_end.t_abs;
            % Is time properly defined ?
            if any(~isfinite(t_end)) || any(t_end<0) || any(t_end>=[24 60 60])
                error('The parameters to determine the beginning of the chunk are not correct.')
            end

            datenum_start=datenum([0 0 0 Dmeg.info.hour]);
            datenum_ench=datenum([0 0 0 t_end]);
            if datenum_ench<datenum_start % assume no recording is >24hours !!!
                datenum_ench=datenum([0 0 1 t_end]);
            end
            datediff = datevec(datenum_ench-datenum_start); %#ok<NASGU>
            Endpts = min(nsamples(Dmeg),round(...
                (datediff(4)*60^2 + datediff(5)*60 + datediff(6))*fsample(Dmeg)));

        % Time defined relative to beginning of data
        %-------------------------------------------
        case 't_rel'
            t_end =  job.chunk(ii).chunk_end.t_rel;
            % Is time properly defined ? Assuming no data is >24h
            if any(~isfinite(t_end)) || any(t_end<0) || any(t_end>=[24 60 60])
                error('The parameters to determine the beginning of the chunk are not correct.')
            end
            Endpts = min(nsamples(Dmeg),round( ...
                (t_end(1)*60^2+t_end(2)*60+t_end(3)) * fsample(Dmeg)));

        % Time defined by an event.
        %--------------------------
        case 't_mark'
            % This bit should be seriously checked!
            % Any one using this way of chunking could report to me please?
            t_mark = job.chunk(ii).chunk_end.t_mark;
            evts = events(Dmeg);
            marker_list = [evts(:).value];
            marker_match = find(marker_list==t_mark.m_type);
            % check if user's selected marker is ok
            if isempty(marker_match) % No match
                error('No matching event/marker value.')
            elseif length(marker_match)<t_mark.m_ind % Index too big
                error('Selected index beyong number of such marker/event.')
            else % ok, so pick time of event
                Endpts = max(1,round( ...
                    evts(marker_match(t_mark.m_ind)).time*fsample(Dmeg)));
            end
    end
    if Endpts<Begpts
        error('Ending time is set before the beginning time!');
    end

    % Check various flags
    flags = struct( ...
        'overwr',job.options.overwr, ... % overwriting previous chunks, or not
        'prefix',job.options.fn_prefix); % define chunk filename prefix
    if isfield(Dmeg,'info') && ~isempty(Dmeg.info) && ...
            isfield(Dmeg.info,'hour') && isfield(Dmeg.info,'date')
        flags.clockt = 1; % use clock time or not
    else
        flags.clockt = 0;
    end
    if flags.overwr
        flags.numchunk = ii+job.options.numchunk-1; % fixing the chunk number
    else
        flags.numchunk = job.options.numchunk; % or letting the routine finding it
    end
    
    % Do it!
    out.chunk(ii).cfiles{1} = crc_process_chunk(Dmeg,Begpts,Endpts, flags);
end

return;
