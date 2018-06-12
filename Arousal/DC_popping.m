function C = DC_popping(C)

% FORMAT DC_popping(args)
% It detects rapid transitions also called poppin artifacts over all good EEG channels  
% The rapid transitions are evaluated by 0.5-s epochs.
%        popping artifacts are saved in C.CRC.DC.shortartf.popping
% INPUT
%       .file   - data file (.mat files)
% 
% NB: it has to be processed after DC_preprocessing, DC_badchannels and DC_pwr
% detection) in every 1-s epoch
% _________________________________________________________________
% Copyright (C) 2014 Cyclotron Research Centre

% Written by D. Coppieters 't Wallant, 2014.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$
% ----------------------

%loading
% if ~isempty(varargin)
%     C = spm_eeg_load(varargin{1});
% else 
%     C = spm_eeg_load;
% end
% parameters depending on sleep recordings
winsize =   30;  
fs      =   fsample(C);
nspl    =   nsamples(C);
Time    =   ceil(nspl/fs);
NofW    =   ceil(nspl/(fs*winsize));

% channels loading (meegchannels do not work because of rereferencing)
channels = meegchannels(C);

ref = find(strncmp(chanlabels(C),'A1',2)); % ref = meegchannels(C,'REF')
popping = [];

for ibc = 1 : NofW
        badchannels{ibc} = union(C.CRC.DC.badchannels.chan_defaillant{ibc},C.CRC.DC.badchannels.chan_incoherent{ibc});
end

% ===========================
% Don't do this @ home :)
fprintf(1,'No worries, everything will be OK: get a coffee instead :D \n');
% ===========================

% Initialization
beta_s = C.CRC.DC.power.beta_s;
spindle = C.CRC.DC.power.spindle_s;
alpha = C.CRC.DC.power.alpha_s;
sp = spindle./(alpha+spindle+beta_s);
msp = max(sp);
rspindle = find(sp/msp>0.85);
h = waitbar(0,'Progression for short artefact detection...');
str = find(beta_s>prctile(beta_s,90));
Rr = zeros(1,Time);
ar = zeros(1,Time);
tr = zeros(1,Time);
locmax = zeros(1,Time);
locmin = zeros(1,Time);  C.CRC.DC.popcar.Rr = Rr;
if ~isempty(ref)
    for is = 1 : Time
        [vmax locmax(is)] = max(C(ref(1),(is-1)*fs+1 : min(nspl,is*fs)));
        [vmin locmin(is)] = min(C(ref(1),(is-1)*fs+1 : min(nspl,is*fs)));
        ar(is) = vmax - vmin;
        tr(is)  = abs(locmax(is)-locmin(is))/fs; 
        Rr(is) = ar(is)/tr(is);
        String  =  ['Progress for popping artefact on active mastoid channel... : ' num2str(100*is/Time) ' %'];
        waitbar(is/Time,h,String);
    end
end
ae = zeros(1,2*Time);
te = zeros(1,2*Time);
Re = zeros(1,2*Time);
for epoch = str
    win = ceil(epoch/winsize);
    gd = setdiff(channels,badchannels{win});
    if ~isempty(gd)
        for it = (epoch-1)*2+1 : min(ceil(2*nspl/fs),epoch*2)
            treeg = C(gd,((it-1)*fs/2+1: min(nspl,it*fs/2)));
            [vmax lmax] = max(treeg,[],2);
            [vmin lmin] = min(treeg,[],2);
            [Re(it) c] = max((vmax-vmin)./(abs(lmax-lmin)/fs));
            ae(it) = vmax(c)-vmin(c);
            te(it) = (abs(lmax(c)-lmin(c))/fs);
            String  =  ['Progress for popping artefact on all EEG channels... : ' num2str(100*epoch/(Time)) ' %'];
            waitbar(epoch/(Time),h,String);
        end   
    end
end  
poppingr = find(or(and(Rr>1e3,ar>20),and(ar>80,Rr>2e2)));
poppinge = unique(ceil(find(and(Re>3e3,or(ae<120,te<0.03)))/2));
poppinge = setdiff(poppinge,rspindle);
popping = union(poppingr,poppinge);
close(h) %enregistrement------------------------------------------------------------
C.CRC.DC.popcar.Rr = Rr;
C.CRC.DC.popcar.Re = Re;
C.CRC.DC.popcar.ar = ar;
C.CRC.DC.popcar.ae = ae;
C.CRC.DC.popcar.te = te;
C.CRC.DC.shortartf.popping = popping;
save(C);
fprintf('* popping artifacts detected per 1-s epoch: %d \n',length(popping))