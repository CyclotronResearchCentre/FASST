function C = DC_arousal(C)
 
% FORMAT DC_emg(args)
% It detects three kinds of artifact:
%     * movement (C.CRC.DC.shortartf.movement)
%     * arousal in non-REM sleep (C.CRC.DC.shortartf.arousal)
%     * arousal in REM sleep (C.CRC.DC.shortartf.arousal_nrem)
%           
% INPUT
%       .file   - data file (.mat files)
% 
% NB: it has to be processed after DC_preprocessing, DC_badchannels and DC_pwr
% detection) in every 1-s epoch, in the following
% _________________________________________________________________
% Copyright (C) 2014 Cyclotron Research Centre

% Written by D. Coppieters 't Wallant, 2014.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$
% ----------------------

% -------------------------------------------------------------------------
%                               Paramètres
% -------------------------------------------------------------------------
% Paramètres généraux
% if ~isempty(varargin)
%     C = spm_eeg_load(varargin{1});
% else 
%     C = spm_eeg_load;
% end
fs      =   fsample(C);    
nspl    =   nsamples(C);
Time    =   ceil(nspl/fs);
% Initialization
mvt_emg     =   C.CRC.DC.artemg;
mt          =   C.CRC.DC.MT;
arouf       =   [C.CRC.DC.arousalr;C.CRC.DC.arousal];
beta_s      =   C.CRC.DC.power.beta_s;
theta_s     =   C.CRC.DC.power.theta_s;
alpha_s     =   C.CRC.DC.power.alpha_s;
spindle_s   =   C.CRC.DC.power.spindle_s;
sp          =   spindle_s./(spindle_s+beta_s+alpha_s);
msp         =   max(sp);
rspindle    =   find(sp/msp>0.85);
arouf       =   setdiff(arouf,rspindle); % returns the data in A that is not in B, with no repetitions. C is in sorted order.

% % INDIVIDUAL AROUSAL
% arousal_r = [];
% ind_events = find(diff(arouf)~=1);
% fin = sort([arouf(ind_events) arouf(end)]); 
% deb = sort([arouf(1) arouf(ind_events+1)]);
% for ic = 1:length(ind_events)
%     arour_ev = deb(ic):fin(ic); % jeweilige Arousalsekunden
%     if length(arour_ev)>3 && isempty(intersect(arour_ev,rspindle)) % wenn länger als 3 sec und KEINE Spindel da, dann Arousal
%         arousal_r = [arousal_r arour_ev];
%     end
% end

% OTHER MOVEMENT ARTEFACTS (according to Brunner et al., 1997) and MT 
artifacts = [];
ind_events = find(diff(mvt_emg)~=1);
fin = sort([mvt_emg(ind_events) mvt_emg(end)]);
deb = sort([mvt_emg(1) mvt_emg(ind_events+1)]);
for ic = 1:length(ind_events)
    ii = 3;
    c = 1;
    bef = [];
    while deb(ic) - ii >0 && c < 11;
        if isempty(intersect(deb(ic)-ii, mvt_emg))
            bef = union(bef,ii);
            c = c+1;
        end 
        ii =ii +1;
    end
    ii = 3;
    c = 1;
    aft = [];
    while fin(ic) + ii <Time && c < 11;
        if isempty(intersect(fin(ic)+ii, mvt_emg))
            aft = union(aft,ii);
            c = c+1;
        end 
        ii =ii +1;
    end
    er_epoch = max(deb(ic)-3,1):min(fin(ic)+3,Time);
    sl = [deb(ic)-bef fin(ic)+aft];  
    act = [];
    % check, if median beta > 3
    act1 = er_epoch((beta_s(er_epoch)/median(beta_s(sl)))>3);
    if ~isempty(act1)
        ind_ev = find(diff(act1)~=1);
        if ~isempty(ind_ev)
            f = sort([act1(ind_ev) act1(end)]);
            d = sort([act1(1) act1(ind_ev+1)]);
            for ie = 1 : length(d)
                evm = d(ie) : f(ie);
                sp = intersect(evm,rspindle);
                if isempty(sp)&&length(evm)>1
                    act = union(act,evm);
                end
            end
        end
    end
    % check, if median theta > 3
    act2 = er_epoch((theta_s(er_epoch)/median(theta_s(sl)))>3);
    if ~isempty(act2)
        ind_ev = find(diff(act2)~=1);
        if ~isempty(ind_ev)
            f = sort([act2(ind_ev) act2(end)]);
            d = sort([act2(1) act2(ind_ev+1)]);
            for ie = 1 : length(d)
                evm = d(ie) : f(ie);
                sp = intersect(evm,rspindle);
                if isempty(sp)&& length(evm)>1
                    act = union(act,evm);
                end
            end
        end
    end
    if ~isempty(act)&& ~isempty(intersect(act,mvt_emg))
        artifacts = [artifacts(:); act(:)];
    end
end

artifacts   =   union(artifacts,mt);
% on va aussi calculer les mauvaises fenêtres via DC_badwindows
C.CRC.DC.shortartf.movement     =   artifacts;
C.CRC.DC.shortartf.arousal      =   C.CRC.DC.arousalr;

%% Write Activation output table (FR)
allT = union(C.CRC.DC.arousal,C.CRC.DC.arousalr);
ind_events = find(diff(allT)~=1);
fin = sort([allT(ind_events); allT(end)]);
deb = sort([allT(1); allT(ind_events+1)]);
    
% Duration longer than 3 seconds?
artal = [];
ita = find((fin-deb+1)>3);
for tt = 1:length(ita) 
    artal = [artal deb(ita(tt)):fin(ita(tt))];
end

% More than 10 seconds apart?
ind_events = find(diff(artal)>10);
tEnd_Sec_a = transpose([artal(ind_events), artal(length(artal))]);
tBeg_Sec_a = transpose([artal(1), artal(ind_events+1)]); 
duration = tEnd_Sec_a - tBeg_Sec_a;
 
% alpha arousal
% Duration longer than 3 seconds?
artalpha = C.CRC.DC.artalpha;
ind_ev = find(diff(artalpha)~=1);
finn = sort([allT(ind_ev); allT(end)]);
debb = sort([allT(1); allT(ind_ev+1)]);
artalpha = [];
ita = find((finn-debb+1)>3);
for tt = 1:length(ita) 
    artalpha = [artalpha debb(ita(tt)):finn(ita(tt))];
end
% More than 10 seconds apart?
ind_ev = find(diff(artalpha)>10);
tEnd_Sec_al = transpose([artalpha(ind_ev), artalpha(length(artalpha))]);
tBeg_Sec_al = transpose([artalpha(1), artalpha(ind_ev+1)]); 

% Add indication for EMG-increase and/or movement relation
% Get info about the epoch
Epoch = floor(tBeg_Sec_a/30)+1;
%Epoch = floor(tBeg_Sec_a-1/30)+1;

EMG = [];
EMG_bef = [];
MT_beta = [];
REM_est = [];
SL_est = [];
SL_est_1 = [];
ALPH = [];
for ii = 1:length(tBeg_Sec_a)
    % is Arousal associated with alpha
    xx = intersect(tBeg_Sec_a(ii):tEnd_Sec_a(ii),tBeg_Sec_al(ii):tEnd_Sec_al(ii));
    if ~isempty(xx)
        ALPH = [ALPH;1];
    else
        ALPH = [ALPH;0];
    end
    % is Arousal associated with concomitant EMG increase
    xx = intersect(tBeg_Sec_a(ii):tEnd_Sec_a(ii),C.CRC.DC.arousalr);
    if ~isempty(xx)
        EMG = [EMG;1];
    else
        EMG = [EMG;0];
    end
    % is Arousal associated with concomitant EMG increase (evaluated
    % with alternative: only tonus before potential arousal)
    xy = intersect(tBeg_Sec_a(ii):tEnd_Sec_a(ii),C.CRC.DC.arousalrALT);
    if ~isempty(xy)
        EMG_bef = [EMG_bef;1];
    else
        EMG_bef = [EMG_bef;0];
    end
    % is Arousal associated with a movement (as evaluated using
    % criteria from Brunner et al., 1997)
    yy = intersect(tBeg_Sec_a(ii):tEnd_Sec_a(ii),C.CRC.DC.shortartf.movement);
    if ~isempty(yy)
        MT_beta = [MT_beta;1];
    else
        MT_beta = [MT_beta;0];
    end
    % occuring during REM or not
    zz = intersect(Epoch(ii),C.CRC.DC.REM);
    if ~isempty(zz)
        REM_est = [REM_est;1];
    else
        REM_est = [REM_est;0];
    end
    % occuring during Sleep or not, excluding stage 1
    zzz = intersect(Epoch(ii),C.CRC.DC.SL);
    if ~isempty(zzz)
        SL_est = [SL_est;1];
    else
        SL_est = [SL_est;0];
    end    
    % occuring during Sleep or not, INCLUDING stage 1
    zzz = intersect(Epoch(ii),C.CRC.DC.SL_incl_1);
    if ~isempty(zzz)
        SL_est_1 = [SL_est_1;1];
    else
        SL_est_1 = [SL_est_1;0];
    end    

end
save(C);
    
% Write Output    
ALL = [tBeg_Sec_a,tEnd_Sec_a,duration,Epoch,REM_est,EMG,EMG_bef,MT_beta,SL_est,SL_est_1];
% ALL = [tBeg_Sec_a,tEnd_Sec_a,duration,Epoch,REM_est,EMG,EMG_bef,MT_beta,SL_est,SL_est_1,ALPH];
ind = findstr('COF',fname(C));
tmp = C.fname;
% tBeg_Sec_a = tBeg_Sec_a-1;%to realign the 1s "delay" between detector and eeg file (and get a second 0)
% tEnd_Sec_a = tEnd_Sec_a-1;%to realign the 1s "delay" between detector and eeg file (and get a second 0)
colonnes = [{'tBeg_Sec_a'} {'tEnd_Sec_a'} {'duration'} {'Epoch'} {'REM_est'} {'EMG'} {'EMG_bef'} {'MT_beta'} {'SL_est'} {'SL_est_1'} ];
% colonnes = [{'tBeg_Sec_a'} {'tEnd_Sec_a'} {'duration'} {'Epoch'} {'REM_est'} {'EMG'} {'EMG_bef'} {'MT_beta'} {'SL_est'} {'SL_est_1'} {'ALPH'} ];
tmp = tmp(ind:ind+5);
xlswrite(strcat(tmp,'_AROUSAL.xls'),colonnes,'Sheet1','A1');
xlswrite(strcat(tmp,'_AROUSAL.xls'),ALL,'Sheet1','A2');

% fprintf('* movement artifacts detected per 1-s epoch: %d \n',length(unique(artifacts)))
% fprintf('* arousal artifacts detected in REM per 1-s epoch: %d \n',length( unique(arousal_r)))
end