function Dmeg = crc_spectcompute(args)

% FORMAT crc_spectcompute(args)
% Compute spectrogram for data using pwelsh and sliding window.
% Power is estimated between a lower (fmin) and upper (fmax) frequency
% bound. 
% To avoid low (DC) drifts and high frequency noise, data are 1st filtered
% between .01 and fmax+1 Hz, with a high-pass followed by a low-pass
% Butter (3rd order) filter.
%
% INPUT
%   args : structure with the following fields
%       .file   - data filename
%       .D      - data object (skipping any provided filename then)
%       .fmax   - max frequency to consider [25Hz, def.]
%       .fmin   - min frequency to consider [.5Hz, def.]
%       .dur    - duration of time window [4s, def.]
%       .step   - time step to use [2s, def.]
%       .scorer - id of night scorer [1, def.]
%__________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by Y. Leclercq & C. Phillips, 2008.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

% get the data
if isfield(args,'D')
    Dmeg = args.D;
    file = fullfile(Dmeg.path,Dmeg.fname);
else
    if isfield(args,'file')
        file = args.file ;
    else
        file = spm_select(1,'mat', 'Select imported EEG file','' ,pwd) ;
    end
    Dmeg = crc_eeg_load(file);
end

% get the parameters & extract.
crcdef = crc_get_defaults('cps');
args = crc_check_flag(crcdef,args);

uplimit = args.fmax ;
downlimit = args.fmin ;
duration = args.dur ;
step = args.step ;
scorer = args.scorer;
reference = args.ref;

% cut off frequencies for pre-filter
fs = fsample(Dmeg);
frqcut(1)= 0.1; % Low cutfrq to suppress DC components
frqcut(2) = uplimit+1;
[Bh,Ah] = butter(3, frqcut(1)/(fs/2), 'high');
[Bl,Al] = butter(3, frqcut(2)/(fs/2), 'low');


% Determine reference
% Coded reference :
% 1 Means original reference
% x (above 1 and below length(D.channels) + 2)
% -1 Mean of Ref and Ref 2
% -2 Mean of M1 & M2
ref1 = [];
ref2 = [];
if (reference == nchannels(Dmeg)+2 && ...
        ~ismember('REF2',upper(chanlabels(Dmeg)))) || ...
        reference == nchannels(Dmeg)+3
    [dumb,ref1]=ismember('M1',upper(chanlabels(Dmeg)));
    [dumb,ref2]=ismember('M2',upper(chanlabels(Dmeg)));
    reference = -2;
elseif reference == nchannels(Dmeg)+2 && ismember('REF2',upper(chanlabels(Dmeg)))
    [dumb,ref2]=ismember('REF2',upper(chanlabels(Dmeg)));
    reference = -1;
end

% List of windows to reject, i.e. with "zero power"
zerolist = [-1 -1];

if isfield(Dmeg,'CRC')
    if isfield(Dmeg.CRC,'score')
        % Remove the artefact & arousal from the computation
        if ~isempty(Dmeg.CRC.score{5,scorer})
            zerolist = [zerolist ; Dmeg.CRC.score{5,scorer}(:,1:2)];
        end
        zerolist=[zerolist ; Dmeg.CRC.score{6,scorer}];
        % Remove the movement time and the 'unscorable' pages from computation
        a = find(Dmeg.CRC.score{1,scorer} ~= 0 & ...
            Dmeg.CRC.score{1,scorer} ~= 1 & ...
            Dmeg.CRC.score{1,scorer} ~= 2 &...
            Dmeg.CRC.score{1,scorer} ~= 3 & ...
            Dmeg.CRC.score{1,scorer} ~= 4  & ...
            Dmeg.CRC.score{1,scorer} ~= 5)'*Dmeg.CRC.score{3,scorer};
        b = a - Dmeg.CRC.score{3,scorer};
        c = [b a];
        zerolist = [zerolist ; c];
    end
end

ideb = 1 ;
ifin = ideb + duration * fs ;

%% Deal with 1st window and create structure in object + .frq file
h = waitbar(0,'Please wait...');
for chan = 1:nchannels(Dmeg)
    Dtmp = getdata(Dmeg,chan,[ideb ifin],reference,ref1,ref2);
    Dtmp = filterlowhigh(Dtmp,Ah,Bh,Al,Bl);
    Dtmp = Dtmp-mean(Dtmp);
    [P,F] = pwelch(Dtmp,ifin-ideb,[],fs*4,fs) ;
    if chan==1
        chosen = find(F>=downlimit & F<=uplimit) ;
        up = find(F > F(max(chosen)));
        down = find(F < F(min(chosen)));
        data = zeros(nchannels(Dmeg),length(chosen)+2);
    end
    data(chan,:) = [sum(P(down)) P(chosen)' sum(P(up))];
end
% to save the data in a frq file -----------------
Dmeg = crc_freqsave_spm(Dmeg,data,reference,scorer);

Check_dat_sz = size(data,3)
Check_dmeg_frqsz = Dmeg.CRC.pwrspect.frqNsamples

%% Deal with the rest
maxmemload  = crc_get_defaults('mem.cps_maxmemload');
maxdouble   = maxmemload / 8 ;  % Maximum doubles in memory
time        = fs * duration ;           %#sample in window to process
unit        = time * nchannels(Dmeg) ;  %#samples to process over a set of channels
maxunit     = floor (maxdouble / unit); %#windows to process, accounting for memory, #channels, #samples/window
maxtime     = maxunit * time ;          %max #samples to produce spect pow, given memory
Nchunks     = ceil(nsamples(Dmeg) / maxtime);

for tt = 1 : Nchunks    % proceed in chunks for memory reasons
    tmp_mem = Dmeg(:, 1 + (tt-1) * maxtime : min (1 + tt * maxtime, nsamples(Dmeg) ) );
    tosub = (tt-1) * maxtime ; % the start of the present chunk
    dattoap = zeros(0,0,0);
    iitime = 1;
    
    while and(ifin + step * fs < nsamples(Dmeg) , ifin - tosub <= maxtime)
        ideb = ifin - step * fs ; % Starting after 1st window!
        ifin = ideb + duration * fs ;
        x = (ideb + ifin)/2 ;
        string = ['Please wait... ' num2str(round(100*x/nsamples(Dmeg))) ' %'];
        waitbar(round(x/nsamples(Dmeg)),h,string)
        concerned = find(or(or(or(... 
            and(ideb/fs > zerolist(:,1),ideb/fs < zerolist(:,2)),...
            and(ifin/fs > zerolist(:,1),ifin/fs < zerolist(:,2))),...
            and(ideb/fs < zerolist(:,1),ifin/fs > zerolist(:,1))),...
            and(ifin/fs > zerolist(:,2),ideb/fs < zerolist(:,2))));
        % List of artefacted window 'covered' by current data window, there
        % could be more than one overlapping...
        if ~isempty(concerned) % trouble for window
            concerned = concerned -1; % QUESTION: why -1 ???
            if concerned(end)>size(Dmeg.CRC.score{5,scorer},1) ...
                    || size(Dmeg.CRC.score{5,scorer},2)<3
                % either 2nd part of zerolist (after the manually defined 
                % artefacts -> all channels) or artefact on all channels
                % (as in old version manually defined artefacts applies on
                % all channels)
                for chan = 1:nchannels(Dmeg)
                    data(chan,:) = 0*[1 P(chosen)' 1];
                end
            else
                if any(Dmeg.CRC.score{5,scorer}(concerned,3)==0)
                    % All channels concerned for this window.
                    for chan = 1:nchannels(Dmeg)
                        data(chan,:) = 0*[1 P(chosen)' 1];
                    end
                else
                    % only single channels are concerned
                    for chan = 1:nchannels(Dmeg)
                        if any(chan == Dmeg.CRC.score{5,scorer}(concerned,3))
                            data(chan,:) = 0*[1 P(chosen)' 1];
                        else
                            if (ideb - tosub) < 0 || ifin-tosub > size(tmp_mem,2)
                                ifin_tmp = min(ifin,nsamples(Dmeg));
                                X = getdata(Dmeg,chan,[ideb ifin_tmp],reference,ref1,ref2);
                                Y = filterlowhigh(X,Ah,Bh,Al,Bl);
                                Y = Y-mean(Y);
                                try
                                    [P,F] = pwelch(Y,ifin-ideb,[],fs*4,fs) ;
                                catch
                                    [P,F] = pwelch(Y,ifin-ideb,[],[],fs) ;
                                end
                            else
                                X = getdata(tmp_mem,chan,[ideb ifin]-tosub,reference,ref1,ref2);
                                Y = filterlowhigh(X,Ah,Bh,Al,Bl);
                                Y = Y-mean(Y);
                                try
                                    [P,F] = pwelch(Y,ifin-ideb,[],fs*4,fs) ;
                                catch
                                    [P,F] = pwelch(Y,ifin-ideb,[],[],fs) ;
                                end
                            end
                            data(chan,:) = [sum(P(down)) P(chosen)' sum(P(up))];
                        end
                    end
                end
            end
        else % no trouble for window
            for chan = 1:nchannels(Dmeg)
                if (ideb - tosub) < 0 || ifin-tosub > size(tmp_mem,2)
                    ifin_tmp = min(ifin,nsamples(Dmeg));
                    X = getdata(Dmeg,chan,[ideb ifin_tmp],reference,ref1,ref2);
                    Y = filterlowhigh(X,Ah,Bh,Al,Bl);
                    Y = Y-mean(Y);
                    try
                        [P,F] = pwelch(Y,ifin-ideb,[],fs*4,fs) ;
                    catch
                        [P,F] = pwelch(Y,ifin-ideb,[],[],fs) ;
                    end
                else
                    X = getdata(tmp_mem,chan,[ideb ifin]-tosub,reference,ref1,ref2);
                    Y = filterlowhigh(X,Ah,Bh,Al,Bl);
                    Y = Y-mean(Y);
                    try
                        [P,F] = pwelch(Y,ifin-ideb,[],fs*4,fs) ;
                    catch
                        [P,F] = pwelch(Y,ifin-ideb,[],[],fs) ;
                    end
                end
                data(chan,:) = [sum(P(down)) P(chosen)' sum(P(up))];
            end
        end
        dattoap(:,:,iitime) = data ;
        iitime = iitime + 1;
    end
    Dmeg = crc_freqappnd_spm(Dmeg,dattoap) ;
    Check_dat_sz(end+1) = size(dattoap,3)
    Check_dmeg_frqsz(end+1) = Dmeg.CRC.pwrspect.frqNsamples
end

ideb = nsamples(Dmeg) - duration * fs ;
ifin = nsamples(Dmeg) ;
string = ['Please wait... ' num2str(100*1/nsamples(Dmeg)) ' %'];
waitbar(1,h,string)

% Last bit
for chan = 1:nchannels(Dmeg)
    Dtmp = getdata(Dmeg,chan,[ideb ifin],reference,ref1,ref2);
    Dtmp = filterlowhigh(Dtmp,Ah,Bh,Al,Bl);
    Dtmp = Dtmp-mean(Dtmp); 
    [P,F] = pwelch(Dtmp,ifin-ideb,[],fs*4,fs) ;
    data(chan,:,1) = [sum(P(down)) P(chosen)' sum(P(up))];
end
Dmeg = crc_freqappnd_spm(Dmeg,data) ;
close(h)

Check_dat_sz(end+1) = size(data,3)
Check_dmeg_frqsz(end+1) = Dmeg.CRC.pwrspect.frqNsamples



% if ~isfield(Dmeg.CRC,'pwrspect')
%     Dmeg.CRC.pwrspect = [];
% end
Dmeg.CRC.pwrspect.frqbins = [-Inf ; F(chosen) ; Inf];
Dmeg.CRC.pwrspect.step = step;
Dmeg.CRC.pwrspect.duration = duration;
save(Dmeg)

return

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filtering the data
function Y = filterlowhigh(X,Ah,Bh,Al,Bl)

% Apply Butterworth filter
Y = filtfilt(Bh,Ah,X);
Y = filtfilt(Bl,Al,Y);

return

% Extracting the data
function Dtmp = getdata(Dmeg,chan,ideb_ifin,reference,ref1,ref2)
% Get the data in time window and specified channel

ideb = ideb_ifin(1);
ifin = ideb_ifin(2);
switch reference
    case 1
        Dtmp = Dmeg(chan,ideb:ifin) ;
    case -1
        Dtmp = mean([Dmeg(chan,ideb:ifin) ; ...
            Dmeg(chan,ideb:ifin)-Dmeg(ref2,ideb:ifin)]);
    case -2
        Dtmp = mean([Dmeg(chan,ideb:ifin)-Dmeg(ref1,ideb:ifin); ...
            Dmeg(chan,ideb:ifin)-Dmeg(ref2,ideb:ifin)]);
    otherwise
        Dtmp = Dmeg(chan,ideb:ifin) - Dmeg(reference-1,ideb:ifin);
end
return