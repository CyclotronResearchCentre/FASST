function Dmeg = crc_spectcompute(args)

% FORMAT crc_spectcompute(args)
% Compute spectrogram for data using pwelsh and sliding window.
%
% INPUT
%   args : structure with the following fields
%       .file   - data file
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

% NOTES for programmers:
% - No filtering on 1st and last chunk!!!
% - BUGS
%   + when deciding the list of frequencies to keep
%       * chosen == list of frequencies between the up/down bounds, OK
%       * up/down == list of freq above/below bounds -> WAS NOT RIGHT!!!
%   + the *old* filter was wrongly setup: it needed a '/fsample(D)' to 
%     define the cut off frequencis in the filtering routine!



crcdef = crc_get_defaults('cps');
try 
    file = args.file ;
catch
    file = spm_select(1,'mat', 'Select imported EEG file','' ,pwd) ;
end

try
    uplimit = args.fmax ;
catch
    uplimit = crcdef.uplimit ;
end

try
    downlimit = args.fmin ;
catch
    downlimit = crcdef.downlimit ;
end

try
    duration = args.dur ;
catch
    duration = crcdef.duration ;
end

try
    step = args.step ;
catch
    step = crcdef.step ;
end

try
    scorer = args.scorer; 
catch
    scorer = crcdef.scorer ;
end

try
    reference = args.ref;
catch
    reference = crcdef.reference;
end

try
    Dmeg = args.D;
catch
    Dmeg = crc_eeg_load(file);
end

% Determine reference
% Coded reference :
% 1 Means original reference
% x (above 1 and below length(D.channels) + 2)
% -1 Mean of Ref and Ref 2
% -2 Mean of M1 & M2

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

fs = fsample(Dmeg);
zerolist=[-1 -1];

if isfield(Dmeg,'CRC')
    if isfield(Dmeg.CRC,'score')  
        % Remove the artefact & arousal from the computation
        if ~isempty(Dmeg.CRC.score{5,scorer})
            zerolist=[zerolist ; Dmeg.CRC.score{5,scorer}(:,1:2)];
        end
        zerolist=[zerolist ; Dmeg.CRC.score{6,scorer}];        
        % Remove the movement time and the 'unscorable' pages from computation
        a = find(Dmeg.CRC.score{1,scorer} ~= 0 & Dmeg.CRC.score{1,scorer} ~= 1 & Dmeg.CRC.score{1,scorer} ~= 2 &...
             Dmeg.CRC.score{1,scorer} ~= 3 & Dmeg.CRC.score{1,scorer} ~= 4  & Dmeg.CRC.score{1,scorer} ~= 5)'*Dmeg.CRC.score{3,scorer};
        b = a - Dmeg.CRC.score{3,scorer};
        c = [b a];
        zerolist = [zerolist ; c];
    end
end


ideb = 1 ;
ifin = ideb + duration * fs ;

h = waitbar(0,'Please wait...');
for chan = 1:nchannels(Dmeg)
    switch reference   
        case 1           
            [P,F] = pwelch(Dmeg(chan,ideb:ifin),ifin-ideb,[],fs*4,fs) ;           
        case -1
            meanofrefdat = mean([Dmeg(chan,ideb:ifin) ; ...
                               Dmeg(chan,ideb:ifin)-Dmeg(ref2,ideb:ifin)]);
            [P,F] = pwelch(meanofrefdat,ifin-ideb,[],fs*4,fs) ;
        case -2
            meanofrefdat = mean([Dmeg(chan,ideb:ifin)-Dmeg(ref1,ideb:ifin); ...
                               Dmeg(chan,ideb:ifin)-Dmeg(ref2,ideb:ifin)]);
            [P,F] = pwelch(meanofrefdat,ifin-ideb,[],fs*4,fs) ;        
        otherwise
            [P,F] = pwelch(Dmeg(chan,ideb:ifin) - ...
                        Dmeg(reference-1,ideb:ifin),ifin-ideb,[],fs*4,fs);            
    end
    chosen = find(F>=downlimit&F<=uplimit) ;
%     up = find(F>max(chosen)); % BUG!
%     down = find(F<min(chosen)); 
    up = find(F>F(max(chosen)));
    down = find(F<F(min(chosen))); 
    data(chan,:) = [sum(P(down)) P(chosen)' sum(P(up))];    
end
% to save the data in a frq file -----------------
Dmeg = crc_freqsave_spm(file,Dmeg,data,reference,scorer);
%%
maxmemload  = crc_get_defaults('mem.cps_maxmemload');
maxdouble   = maxmemload / 8 ;  % Maximum double number in memory
time        = fs * duration ;           %nbre d'echantillons contenu dans une fenetre à analyser
unit        = time * nchannels(Dmeg) ;  %nbre d'echantillons à analyser par tranche de fenetres (donc cumul des diff fenetres)
maxunit     = floor (maxdouble / unit); %combien de fenetres on peut analyser compte tenu de la memoire, du nombre de canaux et du nombre d'échantillons dans chacune de ces fenetres
maxtime     = maxunit * time ;          %le nombre max d'echantillons que l'on peut avoir pour calculer le spect pow par rapport à la memoire
Nchunks     = ceil(nsamples(Dmeg) / maxtime); 
% cut off frequencies
frqcut(1)= 0.1; % Low cutfrq to suppress DC components
frqcut(2) = uplimit+0.25; 

for tt = 1 : Nchunks    %for memory we have to compute the spectral power by chunks 
    tmp_mem = Dmeg(:, 1 + (tt-1) * maxtime : min (1 + tt * maxtime, nsamples(Dmeg) ) ); %tmp of memory
    tosub = (tt-1) * maxtime ; %the start of the present chunk
    dattoap = zeros(0,0,0);
    iitime = 1;
    
    while and (ifin + step * fs < nsamples(Dmeg) , ifin - tosub <= maxtime)   
        ideb = ifin - step * fs ;
        ifin = ideb + duration * fs ;
        x = (ideb + ifin)/2 ;
        string = ['Please wait... ' num2str(100*x/nsamples(Dmeg)) ' %'];
        waitbar(x/nsamples(Dmeg),h,string)  
        concerned = find(or(or(or(...
            and(ideb/fs > zerolist(:,1),ideb/fs < zerolist(:,2)),...
            and(ifin/fs > zerolist(:,1),ifin/fs < zerolist(:,2))),...
                    and(ideb/fs < zerolist(:,1),ifin/fs > zerolist(:,1))),...
                    and(ifin/fs > zerolist(:,2),ideb/fs < zerolist(:,2))));             
        if ~isempty(concerned)
            concerned = concerned -1;
            %old version of arefacts, there is no artefact on single
            %channel ||some kinds of artefact on all channels are concerned
            if concerned(length(concerned))>size(Dmeg.CRC.score{5,scorer},1)||size(Dmeg.CRC.score{5,scorer},2)<3
               for chan = 1:nchannels(Dmeg)
                    data(chan,:) = 0*[sum(P(down)) P(chosen)' sum(P(up))];
               end
            %There are undefined artefacts concerned
            else
                %It's undefined artefact on all channels
                if any(Dmeg.CRC.score{5,scorer}(concerned,3)==0)
                	for chan = 1:nchannels(Dmeg)
                    	data(chan,:) = 0*[sum(P(down)) P(chosen)' sum(P(up))];
                	end
                else 
                    %only single channels are concerned
                    for chan = 1:nchannels(Dmeg)
                    	if any(chan == Dmeg.CRC.score{5,scorer}(concerned,3))
                            data(chan,:) = 0*[sum(P(down)) P(chosen)' sum(P(up))]; 
                        else
                            if (ideb - tosub) < 0 || ifin-tosub > size(tmp_mem,2)
                                switch reference                        
                                    case 1                            
                                        X = Dmeg(chan,ideb:min(ifin,nsamples(Dmeg)));                        
                                    case -1
                                        X = mean( ...
                                            [Dmeg(chan,ideb:min(ifin,nsamples(Dmeg)));...
                                             Dmeg(chan,ideb:min(ifin,nsamples(Dmeg)))-...
                                             Dmeg(ref2,ideb:min(ifin,nsamples(Dmeg)))]);
                                    case -2
                                        X = mean( ...
                                            [Dmeg(chan,ideb:min(ifin,nsamples(Dmeg)))-...
                                             Dmeg(ref1,ideb:min(ifin,nsamples(Dmeg)));...
                                             Dmeg(chan,ideb:min(ifin,nsamples(Dmeg)))-...
                                             Dmeg(ref2,ideb:min(ifin,nsamples(Dmeg)))]);
                                    otherwise
                                        X = Dmeg(chan,ideb:min(ifin,nsamples(Dmeg)))-...
                                            Dmeg(reference-1,ideb:min(ifin,nsamples(Dmeg)));                            
                                end                                        
                                Y = filterlowhigh(Dmeg,X,frqcut);
                                try
                                    [P,F] = pwelch(Y,ifin-ideb,[],fs*4,fs) ;
                                catch
                                    [P,F] = pwelch(Y,ifin-ideb,[],[],fs) ;
                                end
                            else
                                switch reference                        
                                    case 1                        
                                        X = tmp_mem(chan,ideb-tosub:ifin-tosub);                                                        
                                    case -1
                                        X = mean( ...
                                            [tmp_mem(chan,ideb-tosub:ifin-tosub); ...
                                             tmp_mem(chan,ideb-tosub:ifin-tosub)-...
                                             tmp_mem(ref2,ideb-tosub:ifin-tosub)]);
                                    case -2
                                        X = mean( ...
                                            [tmp_mem(chan,ideb-tosub:ifin-tosub)-...
                                             tmp_mem(ref1,ideb-tosub:ifin-tosub); ...
                                             tmp_mem(chan,ideb-tosub:ifin-tosub)-...
                                             tmp_mem(ref2,ideb-tosub:ifin-tosub)]);
                                    otherwise
                                        X = tmp_mem(chan,ideb-tosub:ifin-tosub)-...
                                            tmp_mem(reference-1,ideb-tosub:ifin-tosub);                            
                                end
                                Y = filterlowhigh(Dmeg,X,frqcut); % Yo = OLDfilterlowhigh(Dmeg,X,frqcut);
                                try
                                    [P,F] = pwelch(Y,ifin-ideb,[],fs*4,fs) ; % [Po,Fo] = pwelch(Yo,ifin-ideb,[],fs*4,fs) ;
                                catch
                                    [P,F] = pwelch(Y,ifin-ideb,[],[],fs) ;
                                end
                            end
                            data(chan,:) = [sum(P(down)) P(chosen)' sum(P(up))]; % do(chan,:) = [sum(Po(down)) Po(chosen)' sum(Po(up))];
                        end
                    end
                end
            end
        else
             for chan = 1:nchannels(Dmeg)
                if (ideb - tosub) < 0 || ifin-tosub > size(tmp_mem,2) 
                    switch reference                        
                        case 1                            
                            X = Dmeg(chan,ideb:min(ifin,nsamples(Dmeg)));                        
                        case -1
                            X = mean( ...
                                [Dmeg(chan,ideb:min(ifin,nsamples(Dmeg)));...
                                 Dmeg(chan,ideb:min(ifin,nsamples(Dmeg)))-...
                                 Dmeg(ref2,ideb:min(ifin,nsamples(Dmeg)))]);
                        case -2
                            X = mean( ...
                                [Dmeg(chan,ideb:min(ifin,nsamples(Dmeg)))-...
                                 Dmeg(ref1,ideb:min(ifin,nsamples(Dmeg)));...
                                 Dmeg(chan,ideb:min(ifin,nsamples(Dmeg)))-...
                                 Dmeg(ref2,ideb:min(ifin,nsamples(Dmeg)))]);
                        otherwise
                            X = Dmeg(chan,ideb:min(ifin,nsamples(Dmeg)))-...
                                Dmeg(reference-1,ideb:min(ifin,nsamples(Dmeg)));                            
                    end                                        
                    Y = filterlowhigh(Dmeg,X,frqcut);
                    try
                        [P,F] = pwelch(Y,ifin-ideb,[],fs*4,fs) ;
                    catch
                        [P,F] = pwelch(Y,ifin-ideb,[],[],fs) ;
                    end
                else
                    switch reference                        
                        case 1                        
                            X = tmp_mem(chan,ideb-tosub:ifin-tosub);                                                        
                        case -1
                            X = mean( ...
                                [tmp_mem(chan,ideb-tosub:ifin-tosub); ...
                                 tmp_mem(chan,ideb-tosub:ifin-tosub)-...
                                 tmp_mem(ref2,ideb-tosub:ifin-tosub)]);
                        case -2
                            X = mean( ...
                                [tmp_mem(chan,ideb-tosub:ifin-tosub)-...
                                 tmp_mem(ref1,ideb-tosub:ifin-tosub); ...
                                 tmp_mem(chan,ideb-tosub:ifin-tosub)-...
                                 tmp_mem(ref2,ideb-tosub:ifin-tosub)]);
                        otherwise
                            X = tmp_mem(chan,ideb-tosub:ifin-tosub)-...
                                tmp_mem(reference-1,ideb-tosub:ifin-tosub);                            
                    end
                    Y = filterlowhigh(Dmeg,X,frqcut); % Yo = OLDfilterlowhigh(Dmeg,X,frqcut);
                    try
                        [P,F] = pwelch(Y,ifin-ideb,[],fs*4,fs) ; % [Po,Fo] = pwelch(Yo,ifin-ideb,[],fs*4,fs) ;
                    catch
                        [P,F] = pwelch(Y,ifin-ideb,[],[],fs) ;
                    end
                end
                data(chan,:) = [sum(P(down)) P(chosen)' sum(P(up))]; % do(chan,:) = [sum(Po(down)) Po(chosen)' sum(Po(up))];
            end
        end       
        dattoap(:,:,iitime) = data ;
        iitime = iitime + 1;        
    end
    Dmeg = crc_freqappnd_spm(Dmeg,dattoap) ;   
end

ideb = nsamples(Dmeg) - duration * fs ;
ifin = nsamples(Dmeg) ;
string = ['Please wait... ' num2str(100*1/nsamples(Dmeg)) ' %'];
waitbar(1,h,string)

for chan = 1:nchannels(Dmeg)
        switch reference
            case 1        
                [P,F] = pwelch(Dmeg(chan,ideb:ifin),ifin-ideb,[],fs*4,fs) ;            
            case -1                
                meanrefdat=mean([Dmeg(chan,ideb:ifin);Dmeg(chan,ideb:ifin)-Dmeg(ref2,ideb:ifin)]);
                [P,F] = pwelch(meanrefdat,ifin-ideb,[],fs*4,fs) ;
            case -2
                meanrefdat=mean([Dmeg(chan,ideb:ifin)-Dmeg(ref1,ideb:ifin);...
                    Dmeg(chan,ideb:ifin)-Dmeg(ref2,ideb:ifin)]);
                [P,F] = pwelch(meanrefdat,ifin-ideb,[],fs*4,fs) ;
            otherwise
                [P,F] = pwelch(Dmeg(chan,ideb:ifin)-Dmeg(reference-1,ideb:ifin),ifin-ideb,[],fs*4,fs) ;            
        end
    data(chan,:,1) = [sum(P(down)) P(chosen)' sum(P(up))];    
end
Dmeg = crc_freqappnd_spm(Dmeg,data) ;
close(h)
if ~isfield(Dmeg.CRC,'pwrspect')
    Dmeg.CRC.pwrspect = [];
end
Dmeg.CRC.pwrspect.frqbins = [-Inf ; F(chosen) ; Inf];
Dmeg.CRC.pwrspect.step = step; 
Dmeg.CRC.pwrspect.duration = duration; 
save(Dmeg)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Y = filterlowhigh(D,X,frqcut)

[B,A] = butter(3, ...
            [frqcut(1)/(fsample(D)/2), frqcut(2)/(fsample(D)/2)], ...
            'bandpass');

% Apply Butterworth filter
Y = filtfilt(B,A,X);

return

% % OLD version!!!
% function Y = OLDfilterlowhigh(D,X,frqcut)
% % function Y = filterlowhigh(X,frqcut)
% 
% flc = frqcut(1)/fsample(D);
% fhc = frqcut(2)/fsample(D);
% % flc = frqcut(1); % BUGGY version !!!
% % fhc = frqcut(2);
% 
% k = .7; % cut-off value
% 
% alphal = (1-k*cos(2*pi*flc)-sqrt(2*k*(1-cos(2*pi*flc))-k^2*sin(2*pi*flc)^2))/(1-k);
% alphah = (1-k*cos(2*pi*fhc)-sqrt(2*k*(1-cos(2*pi*fhc))-k^2*sin(2*pi*fhc)^2))/(1-k);
%    
% % Apply low pass filter
% Y = filtfilt(1-alphal,[1 -alphal],X);
% 	
% % Apply high pass filter
% Y = filtfilt(1-alphah,[1 -alphah],X-Y);
% 
% return