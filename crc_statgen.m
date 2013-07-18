function [out_V,fig_z] = crc_statgen(D,i_scorer,fig_z)
% FORMAT [out_V,fig_z] = crc_statgen(D,i_scorer,fig_z)
%
% Build sleep statistics for a (series of) scored data file(s).
% NOTE:
% - By defaults this stat program use the score of the *first* scorer.
% - If a single data object is passed then it generates a plot (as when
%   called from the main disply GUI), otherwise it will generate a text
%   file on disk with the output stats.
%
% INPUT:
%   D        : (array of) data object(s). If empty, select them manually.
%   i_scorer : index of the scorer [1 by default]. A list can be passed to
%              specify a different scorer for each data set.
%   fig_z    : figure index to use [create new one by default].
%
% OUTPUT:
%   out_V    : vector with all the sleep stats.
%   fig_z    : figure index used
%
% Output vector, in this order:
%   W       : time awake
%   S1      : time in stage 1
%   S2      : time in stage 2
%   S3      : time in stage 3
%   S4      : time in stage 4
%   REM     : time in REM
%   MT      : time in "movement time"
%   Unsc    : time in "unscorable"
%   TNotsc  : time not scored
%   LatS1   : latency of stage 1
%   LatS2   : latency of stage 2
%   LatREM  : latency of REM
%   TRS     : time allowed to sleep
%   TPS     : time of sleeping period
%   TST     : total sleep time
%   SEff    : sleep efficiency
%   S1Eff   : stage 1 efficiency
%   S2Eff   : stage 2 efficiency
%   S3Eff   : stage 3 efficiency
%   S4Eff   : stage 4 efficiency
%   REMEff  : REM efficiency
%   Arhour  : arousals per hour
%
% TODO: account for artefacted periods in the estimated stats!
%__________________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by Y. Leclercq & C. Phillips, 2008.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

% Input/load
crcdef = crc_get_defaults('score');

if nargin<2
    i_scorer = 1;
end
save_file = 0;
disp_fig = 0;
if nargin<1 || isempty(D)
    files    = spm_select(Inf,'mat','Select the EEG file(s)');
    for ii=1:size(files,1)
        D{ii} = crc_eeg_load(deblank(files(ii,:)));
    end
    save_file = 1;
    filename = 'Sleep_stat.txt';
    fid      = fopen(filename,'w+');
elseif isa(D,'meeg')
    D = {D};
    disp_fig = 1;
end

if (nargin<3 | fig_z==0) & disp_fig
    fig_z = figure;
else
    try
        close(fig_z)
        figure(fig_z)
    catch
        fig_z = figure;
    end
end
N_data = numel(D);

% Check score and scorer index
if length(i_scorer)==1
    i_scorer = ones(1,N_data)*i_scorer;
elseif length(i_scorer)>1 & length(i_scorer)~=N_data
    error('Mismatch between #data sets and #i_scorer');
end
for ii=1:N_data
    if ~isfield(D{ii},'CRC') | ~isfield(D{ii}.CRC,'score')
        error(['No score available for data : ',D{ii}.fname]);
    end
    if i_scorer(ii)<1 | i_scorer(ii)>size(D{ii}.CRC.score,2)
        error(['Wrong i_scorer for data : ',D{ii}.fname]);
    end
end

% Proceed through the data files.
out_V   = [];
for ii=1:N_data
    iisc = D{ii}.CRC.score{1,i_scorer(ii)};
    try
        FPL = D{ii}.CRC.score{4,i_scorer(ii)}(1);
        OPL = D{ii}.CRC.score{4,i_scorer(ii)}(2);
        Winsize = D{ii}.CRC.score{3,i_scorer(ii)};
    catch
        FPL = D{ii}.CRC.pl(1);
        OPL = D{ii}.CRC.pl(2);
        Winsize = 20;
    end
    
    % Time allowed to sleep
    TRS = OPL - FPL;
    [TRStime TRSstring] = crc_time_converts(TRS);
    
    % Bits not useful
    adapted = 1:length(iisc);
    nottobescored = find(adapted < FPL/Winsize | adapted > OPL/Winsize);
    iisc(nottobescored)=-1;
    
    % Score indexes
    Zero = find(iisc==0);   % awake
    One = find(iisc==1);    % Stage 1
    Two = find(iisc==2);    % Stage 2
    Three = find(iisc==3);  % Stage 3
    Four = find(iisc==4);   % Stage 4
    Five = find(iisc==5);   % REM
    Six = find(iisc==6);    % MT
    Seven = find(iisc==7);  % unscorable
    NotSc = find(isnan(iisc)); % Bits that are not scored or not merged
    
    % Duration of the sleeping period, begin_sleep -> end_sleep
    TPS = (max([Two Three Four Five])-1)*Winsize - ...
        (min([Two Three Four Five])-1)*Winsize;
    if isempty(TPS)
        TPS = nan;
        TPStime = [Nan Nan Nan];
        TPSstring = '';
    else
        [TPStime,TPSstring] = crc_time_converts(TPS);
    end
    
    % Total Sleep Time
    TST = length([Two Three Four Five])*Winsize;
    if isempty(TST)
        TST = nan;
        TSTtime = [Nan Nan Nan];
        TSTstring = '';
    else
        [TSTtime,TSTstring] = crc_time_converts(TST);
    end
    
    % Lantency St1
    if isempty(One)
        LatS1 = nan;
        LatS1string = ['No ',crcdef.stnames_L{2},' Scored'];
    else
        LatS1 = (One(1)-1)*Winsize - FPL;
        [LatS1time,LatS1string] = crc_time_converts(LatS1);
    end
    
    % Lantency St2
    if isempty(Two)
        LatS2 = nan;
        LatS2string = ['No ',crcdef.stnames_L{3},' Scored'];
    else
        LatS2 = (Two(1)-1)*Winsize - FPL;
        [LatS2time,LatS2string] = crc_time_converts(LatS2);
    end
    
    % Lantency REM
    if isempty(Five)
        LatREM = nan;
        LatREMstring = ['No ',crcdef.stnames_L{6},' Scored'];
    else
        LatREM = (Five(1)-1)*Winsize - FPL;
        [LatREMtime,LatREMstring] = crc_time_converts(LatREM);
    end
    
    
    % Min & Percentage of W
    W = length(Zero)*Winsize;
    PW = W/TST;
    [Wtime,Wstring] = crc_time_converts(W);
    
    % Min & Percentage of S1
    S1 = length(One)*Winsize;
    PS1 = S1/TST;
    [S1time,S1string] = crc_time_converts(S1);
    
    % Min & Percentage of S2
    S2 = length(Two)*Winsize;
    PS2 = S2/TST;
    [S2time,S2string] = crc_time_converts(S2);
    
    % Min & Percentage of S3
    S3 = length(Three)*Winsize;
    PS3 = S3/TST;
    [S3time,S3string] = crc_time_converts(S3);
    
    % Min & Percentage of S4
    S4 = length(Four)*Winsize;
    PS4 = S4/TST;
    [S4time,S4string] = crc_time_converts(S4);
    
    % Min & Percentage of SWS
    SWS = length([Three Four])*Winsize;
    [SWStime,SWSstring] = crc_time_converts(SWS);
    PSWS = SWS/TST;
    
    % Min & Percentage of REM
    REM = length(Five)*Winsize;
    PREM = REM/TST;
    [REMtime,REMstring] = crc_time_converts(REM);
    
    % Min & Percentage of MT
    MT = length(Six)*Winsize;
    PMT = MT/TST;
    [MTtime,MTstring] = crc_time_converts(MT);
    
    %Min & Percentage of unscorable
    Unsc = length(Seven)*Winsize;
    PUnsc = Unsc/TST;
    [Unsctime Unscstring] = crc_time_converts(Unsc);
    %warning if 'unscorable' is present in 1/4 of the file
    if Unsc>=(nsamples(D{ii})/fsample(D{ii}))/4
        beep
        warning(['At least 25% of the file is marked as "Unscorable".',...
            '  Please be careful when analysing the results and consider re-scoring'])
    end
    
    %Min & Percentage of 'not scored time'
    TNotSc = length(NotSc)*Winsize;
    PNotSc = TNotSc/TST;
    [Notsctime Notscstring] = crc_time_converts(TNotSc);
    
    
    %Sleep Efficiency
    SEff = TST/TRS;
    
    %S1 Efficiency
    S1Eff = S1/TRS;
    
    %S2 Efficiency
    S2Eff = S2/TRS;
    
    %S3 Efficiency
    S3Eff = S3/TRS;
    
    %S4 Efficiency
    S4Eff = S4/TRS;
    
    %REM Efficiency
    REMEff = REM/TRS;
    
    %Arousal/hour
    nb_ar = size(D{ii}.CRC.score{6,i_scorer(ii)},1);
    Arhour = nb_ar/(TST/60);
    
    A = [W S1 S2 S3 S4 REM MT Unsc TNotSc...        % 9
        round(LatS1) round(LatS2) round(LatREM)...  % 3
        round(TRS) round(TPS) round(TST)...         % 3
        SEff S1Eff S2Eff S3Eff...                   % 4
        S4Eff REMEff Arhour ];                      % 3
    out_V = [out_V ; A]; %#ok<AGROW>
    
    if save_file
        fprintf( fid, [repmat('%s ',1,23) '\n'], 'Name',...
            'W','S1' ,'S2', 'S3', 'S4', 'REM' ,'MT' ,'Unsc', 'Notsc', ...
            'LatS1' ,'LatS2' ,'LatREM' ,'TRS' ,'TPS' ,'TST',...
            'SEff', 'S1Eff', 'S2Eff', 'S3Eff', 'S4Eff', 'REMEff', 'Arhour');
        fprintf(fid,['%s ' repmat('%6.3f ',1,22) '\n'],D{ii}.fname , A);
        fclose(fid);
    end
    if disp_fig
        cmap = [ ...
            0.2 0.75 0.6; ...
            0 0.8 1; ...
            0.1 0.5 0.9; ...
            0.1 0.2 0.8; ...
            0.1 0.15 0.5; ...
            0.5 0.5 0.9; ...
            0.9 0.4 0.4; ...
            0.9 0.6 0.3; ...
            0.9 0.9 0.9];
        legendm = crcdef.stnames_S;
        legendm{end+1} = 'NotSco';
        set(fig_z,'NumberTitle','off')
        set(fig_z,'Name','Scoring Statistics')
        
        try
            subplot(121)
            x = [PW PS1 PS2 PS3 PS4 PREM PMT PUnsc PNotSc];
            explode = [1 1 0 0 0 0 1 1 1];
            x(find(x==0))=eps;
            pie(x,explode)
            legend(legendm,'location','SouthOutside')
            colormap(cmap)
        end
        
        a = subplot(122);
        set(a,'Color',[0.8,0.8,0.8])
        axis off
        text(0 ,1,['Time ',crcdef.stnames_sS{1},': '])
        text(0.6, 1, Wstring )
        text(0 ,0.95,['Time ',crcdef.stnames_sS{2},': '])
        text(0.6, 0.95, S1string )
        text(0 ,0.9,['Time ',crcdef.stnames_sS{3},': '])
        text(0.6, 0.9, S2string)
        text(0 ,0.85,['Time ',crcdef.stnames_sS{4},': '])
        text(0.6, 0.85,S3string)
        text(0 ,0.8,['Time ',crcdef.stnames_sS{5},': '])
        text(0.6, 0.8, S4string )
        text(0 ,0.75,['Time ',crcdef.stnames_sS{6},': '])
        text(0.6, 0.75,REMstring )
        text(0 ,0.7,['Time ',crcdef.stnames_sS{7},': '])
        text(0.6, 0.7,MTstring )
        text(0 ,0.65,['Time ',crcdef.stnames_sS{8},': '])
        text(0.6, 0.65,Unscstring )
        if TNotSc
            text(0 ,0.60,'Time NotScored: ')
            text(0.6, 0.60,Notscstring )
            vshift = .05;
        else
            vshift = 0;
        end
        
        text(0 ,0.55-vshift,[crcdef.stnames_sS{2},' Latency: '])
        text(0.6, 0.55-vshift,LatS1string )
        text(0 ,0.50-vshift,[crcdef.stnames_sS{3},' Latency: '])
        text(0.6, 0.50-vshift,LatS2string )
        text(0 ,0.45-vshift,[crcdef.stnames_sS{6},' Latency: '])
        text(0.6, 0.45-vshift,LatREMstring )
        
        text(0 ,0.35-vshift,'Time Allowed to Sleep: ')
        text(0.8, 0.35-vshift,TRSstring )
        text(0 ,0.3-vshift,'Sleeping Period: ')
        text(0.8, 0.3-vshift,TPSstring )
        text(0 ,0.25-vshift,'Total Sleep Time: ')
        text(0.8, 0.25-vshift,TSTstring )
        
        text(0 ,0.15-vshift,'Sleep efficiency: ')
        text(0.6, 0.15-vshift,[num2str(round(SEff*100)) ' %'])
    end
end

return