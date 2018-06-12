function D = DC_bademg(D)

% FORMAT DC_bademg(args)
% It build a new emg channels from emg channels available
%        new emg is saved in C.CRC.DC.Emg 
% INPUT
%       .file   - data file (.mat files)
% 
% NB: it has to be processed after DC_preprocessing, DC_badchannels, 
% DC_pwr and DC_popping
% _________________________________________________________________
% Copyright (C) 2014 Cyclotron Research Centre

% Written by D. Coppieters 't Wallant, 2014.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$
% ----------------------

% parameters depending on sleep recordings
winsize =   30; 
emg = find(strncmp(chanlabels(D),'Chin EMG',8));
fs      =   fsample(D);
nspl    =   nsamples(D);
Time    =   ceil(nspl/fs);
NofW    =   floor(Time / winsize);
nbr_sc_epoch    =   winsize*fs;
% Threshold fixed 
% tr = crc_get_defaults('qc');
tr_fm   = 0.1;%tr.bc.fm;      
unk = zeros(NofW,1);
h = waitbar(0,'Progression bad emgchannels detection...');
if ~isempty(emg)
    if  length(emg)== 2 % Ask to the user if he uses the bipolar display: if yes => mimic the display 
        check_flat = 0;
        for w = 1 : NofW
            sig_emg = D(emg,(w-1)*nbr_sc_epoch + 1: min(w*nbr_sc_epoch,nspl))- mean(D(emg,(w-1)*nbr_sc_epoch + 1: min(w*nbr_sc_epoch,nspl)),2)*ones(1,size( D(emg,(w-1)*nbr_sc_epoch + 1: min(w*nbr_sc_epoch,nspl)),2));
            ac = median(abs(sig_emg),2);
            if any(ac<tr_fm) && check_flat<4    %emg plat
                if check_flat<3
                    check_flat = check_flat+1; % tracé plat
                else 
                    check_flat = check_flat+1;
                    fprintf('one of both emg is flat. Only one available emg tracks can be used')
                end
                [gum g(w)]=max(ac);
                new_emg((w-1)*nbr_sc_epoch + 1: min(w*nbr_sc_epoch,nspl))  = sig_emg(g(w),:);
                unk(w) = 1;
            elseif ~any(ac<tr_fm) && (max(ac)/min(ac))>2 && check_flat<4   % pas de flat mais une incoherence entre les deux emg channels
                [gum g(w)]=min(ac);
                s = [];
                for itr = 1 : ceil(size(sig_emg,2)/fs) 
                    emg_tr = sig_emg(:,(itr-1)*fs+1 : min(itr*fs,length(sig_emg)));
                    s(:,itr) = mean(abs(emg_tr),2);
                end
                % y a t il des tranches avec des différences importantes?
                ck1 = find((max(s)./min(s))>2);
                gd  = setdiff(1:ceil(length(sig_emg)/fs),ck1);
                if ~isempty(gd) && length(gd)>ceil(size(sig_emg,2)/fs)/2
                  new_emg((w-1)*nbr_sc_epoch+1: min(w*nbr_sc_epoch,nspl))  = (sig_emg(1,:)-sig_emg(2,:));
                    if ((w-1)*winsize+gd(end))==Time
                        gd = gd(1:end-1);
                    end
                    Mnew_emg = zeros(numel(gd),fs);
                    for ick = 1 : numel(gd)                      
                        Mnew_emg(ick,:) = new_emg((w-1)*winsize*fs+(gd(ick)-1)*fs+1:(w-1)*winsize*fs+gd(ick)*fs);
                    end
                    Mtot = mean(Mnew_emg);
                    for ick = 1 : numel(ck1) 
                        new_emg((w-1)*winsize*fs+(ck1(ick)-1)*fs+1 : (w-1)*winsize*fs + ck1(ick)*fs) = Mtot;
                    end
                    unk(w) = 0;
                else 
                    new_emg((w-1)*nbr_sc_epoch + 1: min(w*nbr_sc_epoch,nspl))  = sig_emg(g(w),:);
                    unk(w) = 1;
                end
            elseif check_flat<4   
                 unk(w) = 0;
                 new_emg((w-1)*nbr_sc_epoch + 1: min(w*nbr_sc_epoch,nspl))  = (sig_emg(1,:)-sig_emg(2,:));
            else 
                [gum g(w)] = min(ac);
                new_emg((w-1)*nbr_sc_epoch + 1: min(w*nbr_sc_epoch,nspl))  = sig_emg(g(w),:);
                unk(w) = 1;
            end
            String  =  ['Progression bad emgchannels detection... : ' num2str(w/NofW*100) ' %'];
            waitbar((w/(NofW)),h,String);
        end
    elseif ~isempty(emg)
        new_emg = D(emg,:);
        for w = 1 : NofW
            raw = D(emg,(w-1)*nbr_sc_epoch + 1: min(w*nbr_sc_epoch,nspl));
            new_emg(:,(w-1)*nbr_sc_epoch + 1: min(w*nbr_sc_epoch,nspl))  = raw - mean(raw,2)*ones(1,size(raw,2));
            String  =  ['Progression bad emgchannels detection... : ' num2str(w/NofW*100) ' %'];
            waitbar((w/(NofW)),h,String);
        end
        if size(new_emg,1)>2
            new_emg = mean(new_emg);
        end
        unk = zeros(1,NofW);
    end    
    chg = zeros(1,Time);
    for tr = 1 : Time
        window_emg = new_emg((tr-1)*fs + 1: min(tr*fs,nspl));
        chg(tr) = max(abs(window_emg)); 
        String  =  ['Progress of threshold evaluation... : ' num2str(tr/floor(Time)*100) ' %'];
        waitbar((tr/(floor(Time))),h,String);
    end
    close(h); 
    known = find(~unk);
    if isempty(known) || known(4)>7
        wind = 1:4;
    else 
        wind = known(1:4);
    end
    slides = [];
    for iw = 1 : length(wind)
        slides = [slides (wind(iw)-1)*winsize+1 : wind(iw)*winsize];
    end
    inw = find(abs(zscore(chg(slides)))<3);
    wake_sup = mean(chg(inw))+2*std(chg(inw));
    %% Global changes
    high_emg = find(chg>wake_sup);
    % exclude whole scoring window
    Wart_idx = ceil(high_emg/winsize);
    dWart_idx1 = Wart_idx(1:end-floor(0.50*winsize)+1);
    dWart_idx3 = Wart_idx(floor(0.50*winsize):end);
    dWart_dur  = dWart_idx1 - dWart_idx3;
    MT = unique(Wart_idx(dWart_dur==0)); % Epochen mit > 50% Bewegung
    mt =[];
    for imt = 1 : length(MT)
        mt = [mt (MT(imt)-1)*winsize:min(MT(imt)*winsize,Time)];
    end
   
    %% Local changes
    art_emg = [];
    for w = 1 : NofW
        v_n = chg(max(1,(w-2)*winsize+1):min(Time,(w+1)*winsize))/(median(chg(max(1,(w-2)*winsize+1):min(Time,(w+1)*winsize))));
        high_v = find(v_n>1.5);
        if ~isempty(high_v)
            art_emg = [art_emg high_v+max(1,(w-2)*winsize+1)-1];
        end
    end
    Mvt_emg = union(high_emg,art_emg);
    %% Get additional parameter: intensity and duration
    ind_events = find(diff(Mvt_emg)~=1);
    fin = sort([Mvt_emg(ind_events) Mvt_emg(end)]);
    deb = sort([Mvt_emg(1) Mvt_emg(ind_events+1)]);
    mvtrem = [];
    for im = 1 : length(deb) % für jedes einzelne Bewegungs-event durchgehen
        ii = 3;
        c = 1;
        bef = [];
        while deb(im) - ii >0 && c < 11; % c = 10; nur 10 Sekunden vorher und nachher beurteilen
            if isempty(intersect(deb(im)-ii, Mvt_emg)) % ist Bewegungsbeginn mit schon hohem EMG assoziiert?
                bef = union(bef,ii);
                c = c+1;
            end 
            ii =ii +1;
        end
        ii = 1;
        c = 1;
        aft = [];
        while fin(im) + ii <Time && c < 11;
            if isempty(intersect(fin(im)+ii, Mvt_emg))
                aft = union(aft,ii);
                c = c+1;
            end 
            ii =ii +1;
        end
        gsl = [deb(im)-bef fin(im)+aft]; % 10 Sekunden VOR/NACH "event"
        if ~isempty(gsl)
            m_win = mean(chg(gsl));
            tr_emg =  find(abs(new_emg((deb(im)-1)*fs+1:min(nspl,fin(im)*fs)))>m_win);%ceil(find(abs(new_emg((w_in(im)-1)*winsize*fs+1:min(nspl,w_in(im)*winsize*fs)))>m_win)/fs);
            tr_un = unique(ceil(tr_emg/fs));
            if ~isempty(tr_un)
                dur_tr = [];
                ect_tr = [];
                for it = 1 : length(tr_un)
                    spl = tr_emg(ceil(tr_emg/fs)==tr_un(it));
                    dur_tr(it) = numel(spl);
                    ect_tr(it) = (spl(end)-spl(1));
                end
                art_emg = tr_un(and(dur_tr./ect_tr>0.05,ect_tr./fs>0.6));
                if ~isempty(art_emg)
                    mvtrem = [mvtrem art_emg+deb(im)-1];
                end
            end
        end
    end
    
    %% POWER 
    theta_s = D.CRC.DC.power.theta_s; 
    alpha_s = D.CRC.DC.power.alpha_s;
    spindle_s = D.CRC.DC.power.spindle_s;
    beta_s = D.CRC.DC.power.beta_s;

    %% Abnornal THETA Activity
    theta = theta_s;
    sp = spindle_s./(alpha_s+spindle_s+beta_s);
    msp = max(sp);
    rspindle = find(sp/msp>0.85);
    rtheta2 = [];
    mdt = median(theta);
    for w = 1 : NofW
        epo = (w-1)*winsize+1:min(Time,w*winsize); % 1 Scoring-window
        epoc = setdiff(epo,Mvt_emg);
        ii = 1;
        c = 1;
        bef = [];
        while epo(1)- ii >0 && c < 11;
            if isempty(intersect(epo(1)-ii, Mvt_emg))
                bef = union(bef,ii);
                c = c+1;
            end 
            ii =ii +1;
        end
        ii = 1;
        c = 1;
        aft = [];
        while epo(end) + ii <Time && c < 11 ;
            if isempty(intersect(epo(end)+ii, Mvt_emg))
                aft = union(aft,ii);
                c = c+1;
            end 
            ii =ii +1;
        end
        % 10 Sekunden VOR/NACH jedem scoring window
        sl = [epo(1)-bef epo(end)+aft epoc];
        v_n = theta(epo)/median(theta(sl));%z-score
        high_v = find(and(v_n>2,theta(epo)>mdt));
        if ~isempty(high_v)
            rtheta2 = [rtheta2 high_v+(w-1)*winsize];
        end
    end 
    %%%
    rtheta2 = setdiff(rtheta2,rspindle);
    ind_events = find(diff(rtheta2)~=1);
    artth = [];
    if ~isempty(ind_events)
        fint = sort([rtheta2(ind_events) rtheta2(end)]);
        debt = sort([rtheta2(1) rtheta2(ind_events+1)]);
        itt = find(fint-debt+1>3);
        for tt = 1 : length(itt)
            artth = [artth debt(itt(tt)) : fint(itt(tt))];
        end  
    end
    
    %% Abnormal BETA activity    
    rbeta2 = [];
    ind = isnan(beta_s);
    ind2 = find(ind == 0);
    big_beta = find(beta_s>prctile(beta_s(ind2),99));
    beta_s(big_beta) = median(beta_s(ind2));
    big_beta = union(big_beta,find(zscore(beta_s(ind2))>3));
    mdb = median(beta_s(ind2));
    for w = 1 : NofW
        epo = (w-1)*winsize+1:min(Time,w*winsize);
        epoc = setdiff(epo,Mvt_emg);
        ii = 1;
        c = 1;
        bef = [];
        while epo(1)- ii >0 && c < 11;
            if isempty(intersect(epo(1)-ii, Mvt_emg))
                bef = union(bef,ii);
                c = c+1;
            end 
            ii =ii +1;
        end
        ii = 1;
        c = 1;
        aft = [];
        while epo(end) + ii <Time && c < 11;
            if isempty(intersect(epo(end)+ii, Mvt_emg))
                aft = union(aft,ii);
                c = c+1;
            end 
            ii =ii +1;
        end
        sl = [epo(1)-bef epo(end)+aft epoc];
        v_n = beta_s(epo)/median(beta_s(sl));%zscore
        high_v = find(and(v_n>2,beta_s(epo)>mdb));
        if ~isempty(high_v)
            rbeta2 = [rbeta2 high_v+(w-1)*winsize];
        end
    end
    %%%
    rbeta2 = union(rbeta2,big_beta);
    rbeta2 = setdiff(rbeta2,rspindle);
    ind_events = find(diff(rbeta2)~=1);
    artbe = [];
    if ~isempty(ind_events)
        finb = sort([rbeta2(ind_events) rbeta2(end)]);
        debb = sort([rbeta2(1) rbeta2(ind_events+1)]);
        itb = find(finb-debb+1>3);
        for tt = 1 : length(itb)
            artbe = [artbe debb(itb(tt)) : finb(itb(tt))];
        end
    end
    
    %% Abnormal activity in ALPHA band
    ralpha2 = [];
    mda = median(alpha_s);
    for w = 1 : NofW
        epo = (w-1)*winsize+1:min(Time,w*winsize);
        epoc = setdiff(epo,Mvt_emg);
        ii = 1;
        c = 1;
        bef = [];
        while epo(1)- ii >0 && c < 11;
            if isempty(intersect(epo(1)-ii, Mvt_emg))
                bef = union(bef,ii);
                c = c+1;
            end 
            ii =ii +1;
        end
        ii = 1;
        c = 1;
        aft = [];
        while epo(end) + ii <Time && c < 11;
            if isempty(intersect(epo(end)+ii, Mvt_emg))
                aft = union(aft,ii);
                c = c+1;
            end 
            ii =ii +1;
        end
        sl = [epo(1)-bef epo(end)+aft epoc];  
        v_n = alpha_s(epo)/median(alpha_s(sl));%zscore
        high_v = find(and(v_n>2,alpha_s(epo)>mda));
        bdn = high_v;
        if ~isempty(high_v)
            ralpha2 = [ralpha2 bdn+(w-1)*winsize];
        end
    end
    %%%%
    ralpha2 = setdiff(ralpha2,rspindle); 
    ind_events = find(diff(ralpha2)~=1);
    artal = [];
    if ~isempty(ind_events)
        fina = sort([ralpha2(ind_events) ralpha2(end)]);
        deba = sort([ralpha2(1) ralpha2(ind_events+1)]);
        ita = find(fina-deba+1>3);
        for tt = 1 : length(ita)
            artal = [artal deba(ita(tt)) : fina(ita(tt))];
        end   
    end
    
    %% Combination
    arousal_f = unique([artth, artbe, artal]); 
    ind_events = find(diff(arousal_f)~=1);
    % EEG shift only
    arousal_nr =[];
    % EEG and EMG shift
    arou_chk = [];
    arou_chk_alt = [];
    if  ~isempty(ind_events)
        fin = sort([arousal_f(ind_events) arousal_f(end)]);
        deb = sort([arousal_f(1) arousal_f(ind_events+1)]);
        for ar = 1 : length(deb)
            arousal_nr = union(arousal_nr,deb(ar):fin(ar));
            pot_arou = deb(ar):fin(ar);
            ii = 3;
            c = 1;
            bef = [];
            while deb(ar)- ii >0 && c < 11;
                if isempty(intersect(deb(ar)-ii, Mvt_emg))
                    bef = union(bef,ii);
                    c = c+1;
                end 
                ii =ii +1;
            end
            ii = 3;
            c = 1;
            aft = [];
            while fin(ar) + ii <Time && c < 11;
                if isempty(intersect(fin(ar)+ii, Mvt_emg))
                    aft = union(aft,ii);
                    c = c+1;
                end 
                ii =ii +1;
            end
            sl = [deb(ar)-bef fin(ar)+aft];
            tonus_around = median(chg(sl));
            tonus_arou = median(chg(pot_arou));
            if tonus_arou>tonus_around*1.5;
                arou_chk = [arou_chk  pot_arou];
            end
            % ------------------------------------------------
            % Alternative: with only tonus BEFORE Arousal (FR)
            tonus_before = median(chg(deb(ar)-bef));
            if tonus_arou>tonus_before*1.5;
                arou_chk_alt = [arou_chk_alt pot_arou];
            end
            % ------------------------------------------------
        end
    end
    D.CRC.DC.artemg             = mvtrem;
    D.CRC.DC.MT                 = mt;
    D.CRC.DC.arousalr           = transpose(arou_chk); 
    D.CRC.DC.arousalrALT        = transpose(arou_chk_alt); 
    D.CRC.DC.arousal            = arousal_nr; 
    D.CRC.DC.Emg                = new_emg;
    D.CRC.DC.artalpha           = artal;
    save(D);
else
        fprintf(1,'-------------------------------------------------------  \n      error : Data no preprocessed !!!     \n------------------------------------------------------- \n');
end