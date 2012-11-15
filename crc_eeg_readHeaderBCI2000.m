function [fdata,ev,param,states] = crc_eeg_readHeaderBCI2000(fname)
% Reads in EEG data from a BCI2000 data file and returns:
% - the path to a binary data file containing only the data
% - a structure array with all the "events" as extracted from the 'states'
%   variable of BCI2000
% - a structure with all the parameters from the BCI2000 header
% - the states variables themselves
%
% Data are saved in 'float' format, so no need anymore for an offset or
% scale factor.
%__________________________________________________________________
% Copyright (C) 2012 Cyclotron Research Centre

% Written by C. Phillips, 2012.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

if nargin<1
    fname = spm_select(1,'^.*\.dat$','Select BCI2000 data file');
end
[pth,fn_orig,ext] = fileparts(fname);
fn_info = dir(fname);

%% Read in text part
ii_line = 1;
go_on = true;
fid = fopen(fname);
while go_on
    tmp = textscan(fid,'%s',1,'delimiter','\n','whitespace','');
    file{ii_line,1} = tmp{1}{1};
    if isempty(file{ii_line})
        go_on = false;
    end
    ii_line = ii_line+1;
end

%% Check first line
if strfind(file{1},'BCI2000V')
    % Version 1.1 or above
    p_DForm = strfind(file{1},'DataFormat=');
    DataFormat = file{1}(p_DForm+12:end); % to be read from 1st line
    switch DataFormat
        case 'int16'
            DataNbytes = 2;
        case 'float32'
            DataNbytes = 4;
        otherwise
            error('Data format unknown -> update routine!');
    end
    p_SvLen = strfind(file{1},'StatevectorLen=');
    StatevectorLen = str2num(file{1}(p_SvLen+15:p_DForm-2));
    
else
    % Version 1.0
    DataFormat = 'int16';
    DataNbytes = 2;
    p_SvLen = strfind(file{1},'StatevectorLen=');
    StatevectorLen = str2num(file{1}(p_SvLen+15:end));
end

p_HLen = strfind(file{1},'HeaderLen=');
p_SCh = strfind(file{1},'SourceCh=');
param.HeaderLen = str2num(file{1}(p_HLen+10:p_SCh-1));
param.NSourceCh = str2num(file{1}(p_SCh+9:p_SvLen-1));
param.StatevectorLen = StatevectorLen;
param.Nsamples = (fn_info.bytes - param.HeaderLen) / ...
    (DataNbytes*param.NSourceCh + param.StatevectorLen);
param.DataFormat = DataFormat;
param.DataNbytes = DataNbytes;

%% Find delimiters for State and Parameters definition
i_StateDef = find( strncmp('[ State Vector Definition ]',file, ...
    length('[ State Vector Definition ]')));
i_ParamDef = find( strncmp('[ Parameter Definition ] ',file, ...
    length('[ Parameter Definition ] ')));

%% Extract the information from State/Parameters into structures.
States_info(i_ParamDef-1-i_StateDef) = struct('name','','type','', ...
    'Nbytes',[]);
is = 1;
for ii = (i_StateDef+1):(i_ParamDef-1)
    l_space = strfind(file{ii},' ');
    States_info(is).name = file{ii}(1:l_space(1)-1);
    Nbits = str2num(file{ii}(l_space(1)+1:l_space(2)-1));
    States_info(is).Nbits = Nbits;
    States_info(is).type = ['ubit',num2str(Nbits)];
    States_info(is).Nbytes = Nbits/8;
    is = is+1;
end
SV_BitLen = sum([States_info(:).Nbits]);
dNbits = param.StatevectorLen*8-SV_BitLen;
if dNbits>0
    % Still read in the last few bits
    States_info(is).name = 'UnknownLastOne';
    States_info(is).type = ['ubit',num2str(dNbits)];
    States_info(is).Nbytes = dNbits/8;
    States_info(is).Nbits = dNbits;
elseif dNbits<0
    error('Looks like #bits of of States > than planned!');
end

E_para = struct('Section','','Type','', 'DefaultValue','', ...
    'LowRange','','HighRange','','Comment','', ...
    'Value','','NumericValue',[]);
is = 1;
param.ErrorListReadHeader = [];
param.ErrorNameReadHeader = [];
for ii = (i_ParamDef+1):(ii_line-2)
    l_space = strfind(file{ii},' ');
    l_bb = strfind(file{ii},'//');
    S_name = file{ii}(l_space(2)+1:l_space(3)-2);
    ii_para = E_para;
    ii_para.Section = file{ii}(1:l_space(1)-1);
    ii_para.Type = file{ii}(l_space(1)+1:l_space(2)-1);
    switch ii_para.Type
        
        case {'int','float'}
            val = file{ii}(l_space(3)+1:l_space(4)-1);
            ii_para.Value = {val};
            ii_para.NumericValue = str2num(val);
            if numel(l_space)>4
                ii_para.DefaultValue = file{ii}(l_space(4)+1:l_space(5)-1) ;
                ii_para.LowRange = file{ii}(l_space(5)+1:l_space(6)-1) ;
                ii_para.HighRange = file{ii}(l_space(6)+1:l_space(7)-1) ;
            end
            
        case 'string'
            ii_para.DefaultValue = file{ii}(l_space(4)+1:l_space(5)-1) ;
            ii_para.LowRange = file{ii}(l_space(5)+1:l_space(6)-1) ;
            ii_para.HighRange = file{ii}(l_space(6)+1:l_space(7)-1) ;
            ii_para.Value = {file{ii}(l_space(3)+1:l_space(4)-1)};
            ii_para.NumericValue = NaN;
            
        case 'matrix'
            try
                sz = [str2num(file{ii}(l_space(3)+1:l_space(5)-1))];
                p_sh = 5 + prod(sz);
                ii_para.DefaultValue = file{ii}(l_space(p_sh)+1:l_space(p_sh+1)-1) ;
                ii_para.LowRange = file{ii}(l_space(p_sh+1)+1:l_space(p_sh+2)-1) ;
                ii_para.HighRange = file{ii}(l_space(p_sh+2)+1:l_space(p_sh+3)-1) ;
                ijm = 1;
                val = cell(sz);
                valn = zeros(sz);
                for im = 1:sz(1)
                    for jm = 1:sz(2)
                        val{im,jm} = file{ii}(l_space(4+ijm)+1:l_space(4+ijm+1)-1) ;
                        vn = str2num(val{im,jm});
                        if numel(vn)==0
                            valn(im,jm) = NaN ;
                        else
                            valn(im,jm) = vn ;
                        end
                        ijm = ijm+1;
                    end
                end
                ii_para.Value = val;
                ii_para.NumericValue = valn;
            catch
                warning(['Sorry couldn''t read line ',num2str(ii), ...
                    ' of header!']);
                ii_para.Value = file{ii};
                ii_para.NumericValue = 'Error reading in!';
                param.ErrorListReadHeader = [param.ErrorListReadHeader ii];
                param.ErrorNameReadHeader = strvcat( ...
                    param.ErrorNameReadHeader,S_name);
            end
            
        case {'floatlist','intlist'}
            Nel = str2num(file{ii}(l_space(3)+1:l_space(4)-1));
            ii_para.DefaultValue = file{ii}(l_space(4+Nel)+1:l_space(5+Nel)-1) ;
            ii_para.LowRange = file{ii}(l_space(5+Nel)+1:l_space(6+Nel)-1) ;
            ii_para.HighRange = file{ii}(l_space(5+Nel)+1:l_space(7+Nel)-1) ;
            val = cell(Nel,1);
            valn = zeros(Nel,1);
            for im = 1:Nel
                val{im} = file{ii}(l_space(3+im)+1:l_space(3+im+1)-1) ;
                valn(im) = str2num(val{im}) ;
            end
            ii_para.Value = val;
            ii_para.NumericValue = valn;
            
        case 'list'
            Nel = str2num(file{ii}(l_space(3)+1:l_space(4)-1));
            val = cell(Nel,1);
            for im = 1:Nel
                val{im} = file{ii}(l_space(3+im)+1:l_space(3+im+1)-1) ;
            end
            ii_para.Value = val;
            
        otherwise
            error('Parameter type not accounted for -> update routine');
    end
    if ~isempty(l_bb)
        ii_para.Comment = file{ii}(l_bb+3:end);
    end
    param.(S_name) = ii_para;
    is = is+1;
end

%% Read in the data and 'states' and write data out cleanly
fdata = fullfile(pth,['spm8_',fn_orig,'.dat']);
fo = fopen(fdata,'w','ieee-le');
for ii=1:numel(States_info)
    states.(States_info(ii).name) = zeros(param.Nsamples,1);
end
for ii = 1:param.Nsamples
    tmp_d = fread(fid,param.NSourceCh,DataFormat) ;
    tmp_d = (tmp_d + param.SourceChOffset.NumericValue) ... % Raw offset
        .* param.SourceChGain.NumericValue; % gain conversion to µV
    % NOTE : note 100% sure about the offset, it could be '-'
    % instead of '+'. Should check this with BCI2000 display.
    fwrite(fo,tmp_d,'float32');
    
    for jj=1:numel(States_info)
        tmp_e = fread(fid,1,States_info(jj).type);
        states.(States_info(jj).name)(ii) = tmp_e;
    end
end
fclose(fid);
fclose(fo);

%% Deal with the stimuli, in order no to return states time series.
Srate = param.SamplingRate.NumericValue;
ev = struct('type','','value',[],'time',[],'duration',[],'offset',0);
i_ev = 1;
field_n = fieldnames(states);
for ii=1:numel(field_n)
    if isempty(strfind(field_n{ii},'Time')) ... % not continuous time signal
            && ~strcmp(field_n{ii},'StimulusCode') % not stim code
        % look for signal change!
        sig = states.(field_n{ii});
        d_sig = diff(sig);
        lpos = find(d_sig>0);
        lneg = find(d_sig<0);
        if isempty(lpos) && isempty(lneg)
            ev(i_ev).type = field_n{ii};
            ev(i_ev).value = ii;
            ev(i_ev).time = 1/Srate;
            ev(i_ev).duration = param.Nsamples/Srate;
            ev(i_ev).offset = 0;
            i_ev = i_ev+1;
        elseif strfind(field_n{ii},'Stimulus') % deal with stim + code
            % NOTE: from the example dataset, it looks like a stimulus is
            % linked to a 1->0 change in signal, hence list of 'lneg'
            for jj=1:numel(lneg)
                code = states.StimulusCode(lneg(jj)+1);
                ev(i_ev).type = ['Stimulus ',num2str(code)];
                ev(i_ev).value = 100 + code;
                ev(i_ev).time = (lneg(jj)+1)/Srate;
                if ~isempty(lpos) && numel(lpos)>=jj
                    ev(i_ev).duration = (lpos(jj)-lneg(jj))/Srate;
                else
                    ev(i_ev).duration = (param.Nsamples-lneg(jj))/Srate;
                end
                ev(i_ev).offset = 0;
                i_ev = i_ev+1;
            end
        else % deal with others
            % NOTE: from the example dataset, it looks like a other events
            % are linked to a 0->1 change in signal, so use 'lpos' list
            for jj=1:numel(lpos)
                ev(i_ev).type = field_n{ii};
                ev(i_ev).value = ii;
                ev(i_ev).time = (lpos(jj)+1)/Srate;
                if ~isempty(lneg) && numel(lneg)>=jj
                    ev(i_ev).duration = (lneg(jj)-lpos(jj))/Srate;
                else
                    ev(i_ev).duration = (param.Nsamples-lpos(jj))/Srate;
                end
                ev(i_ev).offset = 0;
                i_ev = i_ev+1;
            end
        end
    end
end

return

% % display stims
% figure
% for ii=1:numel(field_n)
%     subplot(ceil(numel(field_n)/4),4,ii)
%     plot(single(states.(field_n{ii})))
%     title(field_n{ii})
% end

% Check file reading
% ftell(fid)
% ftell(fid)-param.HeaderLen
