function [header,marker] = crc_eeg_readHeaderBrainP(fname)

%_________________________________________________________________________________________
%
% FORMAT [header,marker] = crc_eeg_read_HeaderBrainP(fname)
%
% 'crc_eeg_read_HeaderBrainP' reads the header of BrainProduct data.
% It produces 2 separate structures: header and marker
%
% Header has the following fields: 
% DataFile, MarkerFile, DataFormat, DataOrientation, DataType, 
% NumberOfChannels, SamplingInterval, Averaged, AveragedSegments, 
% SegmentDataPoints, SegmentationType, DataPoints, binary/ascii, channels, 
% info.
% All fields are explained in the recorder.pdf manual.
%
% Marker has the following fields: 
% DataFile, nMrk, Type, Descript, Position, Points, ChanNumb, Date.
%
% Interstingly, ALL data have a marker file, with at least the marker type 
% 'New Segment', with the date and time of the beginning of recording.
%__________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by Y. Leclercq & C. Phillips, 2008.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

% Originally written by C. Phillips, 2004.10.20

if nargin<1
    fname = spm_select(1,'^.*\.vhdr$','EEG data header file');
end

path_file = spm_str_manip(fname,'H');

fid = fopen(fname);
tmp = textscan(fid,'%s','delimiter','\n','whitespace','');
fclose(fid);
file = tmp{1};

%#ok<*ST2NM>
%%
%========================
%
% READ HEADER FILE
%
% Start with the main fields of header
%===========================
l_fields = {'DataFile','MarkerFile','DataFormat','DataOrientation', ...
            'DataType','NumberOfChannels','SamplingInterval','Averaged', ...
            'AveragedSegments','SegmentDataPoints','SegmentationType', ...
            'DataPoints'};
fieldval = {'','','ASCII','Multiplexed','TimeDomain',[],[],'No',0,0, ...
            'NotSegmented',0};

for ii=1:length(l_fields)
    pfield = find(strncmp(file,l_fields{ii},length(l_fields{ii})));
    if ~isempty(pfield)
        if ischar(fieldval{ii})
            fieldval{ii} = file{pfield}((findstr(file{pfield},'=')+1):end);
        else
            fieldval{ii} = str2num(file{pfield}((findstr(file{pfield},'=')+1):end));
        end
    end
end

% Just check if there are any placeholder, stated by '$b', and add path info
for ii=1:2
    p_holder = findstr(fieldval{ii},'$b');
    if ~isempty(p_holder)
        fieldval{ii} = [spm_str_manip(fname,'r'),fieldval{ii}((p_holder+2):end)];
    end
    if strcmp(spm_str_manip(fieldval{ii},'H'),'.')
        fieldval{ii} = fullfile(path_file,fieldval{ii});
    end
end
header = cell2struct(fieldval,l_fields,2);

% Deal with extra ASCII/BINARY-format information.
%=================================================
if strcmp(header.DataFormat,'ASCII')
    ascii_f = {'DecimalsSymbol','SkipLines','SkipColumns'};
    ascii_v = {'.',0,0};
    for ii=1:length(ascii_v)
        pfield = find(strncmp(file,ascii_f{ii},length(ascii_f{ii})));
        if ~isempty(pfield)
            if ischar(ascii_v{ii})
                ascii_v{ii} = file{pfield}((findstr(file{pfield},'=')+1):end);
            else
                ascii_v{ii} = str2num(file{pfield}((findstr(file{pfield},'=')+1):end));
            end
        end
    end
    header.ascii = cell2struct(ascii_v,ascii_f,2);
elseif strcmp(header.DataFormat,'BINARY')
    bin_f = {'BinaryFormat','ChannelOffset','DataOffset', ...
            'SegmentHeaderSize','TrailerSize','UseBigEndianOrder'};
    bin_v = {'INT_16',0,0,0,0,'No'};
    for ii=1:length(bin_v)
        pfield = find(strncmp(file,bin_f{ii},length(bin_f{ii})));
        if ~isempty(pfield)
            if ischar(bin_v{ii})
                bin_v{ii} = file{pfield}((findstr(file{pfield},'=')+1):end);
            else
                bin_v{ii} = str2num(file{pfield}((findstr(file{pfield},'=')+1):end));
            end
        end
    end
    header.binary = cell2struct(bin_v,bin_f,2);
end

% Deal with the channels information
%===================================
pChanInfo = find(strcmp(file,'[Channel Infos]'));
p1stChan = find(strncmp(file,'Ch',2));
p1stChan = min(p1stChan(find(p1stChan>pChanInfo)))-1;

Cnames = cell(1,header.NumberOfChannels);
RefChanName = cell(1,header.NumberOfChannels);
Resolution = zeros(header.NumberOfChannels,1);
CUnit = cell(1,header.NumberOfChannels);
for ii=1:header.NumberOfChannels
    pequal = findstr(file{ii+p1stChan},'=');
    pcomma = findstr(file{ii+p1stChan},',');
    Cnames{ii} = file{ii+p1stChan}((pequal+1):(pcomma(1)-1));
    RefChanName{ii} = file{ii+p1stChan}((pcomma(1)+1):(pcomma(2)-1));
    if length(pcomma)<3
        % old vhdr, only 2 commas !
        % Leave Resolution empty and assume units is µV
        res_tmp = str2num(file{ii+p1stChan}((pcomma(2)+1):end));
        if ~isempty(res_tmp)
            Resolution(ii) = res_tmp;
        end
        CUnit{ii} = 'µV';
    else
        res_tmp = str2num(file{ii+p1stChan}((pcomma(2)+1):pcomma(3)-1));
        if ~isempty(res_tmp)
            Resolution(ii) = res_tmp;
        end
        CUnit{ii} = file{ii+p1stChan}((pcomma(3)+2):end);
        % Somehow 1st character after last comma is 'weird', just drop it.
    end
end

header.channels.Cnames = Cnames ;
header.channels.RefChanName = RefChanName;
header.channels.Resolution = Resolution;
header.channels.CUnit = CUnit;

% Sometimes electr coordinates are provided as well
%=================================================
pCoord = find(strcmp(file,'[Coordinates]'));
if pCoord
    p1stChan = find(strncmp(file,'Ch',2));
    p1stChan = min(p1stChan(find(p1stChan>pCoord)))-1;
    r_th_ph = zeros(header.NumberOfChannels,3);
    d2r = pi/180;
    for ii=1:header.NumberOfChannels
        pequal = findstr(file{ii+p1stChan},'=');
        pcomma = findstr(file{ii+p1stChan},',');
        r_th_ph(ii,1) = str2num(file{ii+p1stChan}((pequal+1):(pcomma(1)-1)));
        r_th_ph(ii,2) = str2num(file{ii+p1stChan}((pcomma(1)+1):(pcomma(2)-1)))*d2r;
        r_th_ph(ii,3) = str2num(file{ii+p1stChan}((pcomma(2)+1):end))*d2r;
    end
    header.channels.r_th_ph = r_th_ph;
end

% Select the rest of the file as comment
%=======================================
pComment = find(strcmp(file,'[Comment]'));
header.info.comment = char(file{(pComment+2):end});
[v,r] = crc_fasst_utils('Ver',mfilename);
header.info.ver = struct('v_nr',v,'rel',r);

% Save file name in .info
header.info.fname = fname;

%%
%========================
%
% READ MARKER FILE
%
%========================
fid = fopen(header.MarkerFile);
tmp = textscan(fid,'%s','delimiter','\n','whitespace','');
fclose(fid);
file = tmp{1};

% Deal with the common information
%===================================
pDataFile = find(strncmp(file,'DataFile',8));
marker.DataFile = file{pDataFile}((findstr(file{pDataFile},'=')+1):end);
if strcmp(spm_str_manip(marker.DataFile,'H'),'.')
    marker.DataFile = fullfile(path_file,marker.DataFile);
end

% Deal with the markers information
%===================================
p1stMrk = find(strncmp(file,'Mk',2));
marker.nMrk = length(p1stMrk);
p1stMrk = min(p1stMrk)-1;

Type = cell(1,marker.nMrk);
Descript = cell(1,marker.nMrk);
Position = zeros(marker.nMrk,1);
Points = zeros(marker.nMrk,1);
ChanNumb = zeros(marker.nMrk,1);
Date = cell(1,marker.nMrk);
for ii=1:marker.nMrk 
    pequal = findstr(file{ii+p1stMrk},'=');
    pcomma = findstr(file{ii+p1stMrk},',');
    Type{ii} = file{ii+p1stMrk}((pequal+1):(pcomma(1)-1));
    Descript{ii} = file{ii+p1stMrk}((pcomma(1)+1):(pcomma(2)-1));
    Position(ii) = str2num(file{ii+p1stMrk}((pcomma(2)+1):(pcomma(3)-1)));
    Points(ii) = str2num(file{ii+p1stMrk}((pcomma(3)+1):(pcomma(4)-1)));
    if length(pcomma)==5
        ChanNumb(ii) = str2num(file{ii+p1stMrk}((pcomma(4)+1):(pcomma(5)-1)));
        Date{ii} = file{ii+p1stMrk}((pcomma(5)+1):end);
    else
        ChanNumb(ii) = str2num(file{ii+p1stMrk}((pcomma(4)+1):end));
    end
end

marker.Type = Type;
marker.Descript = Descript;
marker.Position = Position;
marker.Points = Points;
marker.ChanNumb = ChanNumb;
marker.Date = Date;

return