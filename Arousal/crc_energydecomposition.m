function [burst teo] = crc_energydecomposition(sf,nrem_est,fs,cfg)

% parameters
tr_tinf = 0.4; % minimal time to consider a sigma burst
tr_tsup = 5;   % maximal time to consider a sigma burst
winsize = 30;
nspl = numel(sf);

% Filtering with the Gabor filter
% ss = DC_Gabor2(sf,fband,fc,fs);
[Bhs Ahs] = butter(3,cfg.hpfreq/(fs/2),'high');
[Bls Als] = butter(3,cfg.lpfreq /(fs/2),'low');
sh = filtfilt(Bls,Als,sf);
ss = filtfilt(Bhs,Ahs,sh);
% TEO
teo = DC_TEO(ss);

% Moving average avoid discontinuity due to jump in the signal (popping artifact) (Potamanios 1995)                        
 for ima = 1 : 6
    teo = MA(teo,1);
 end
 
high_sig = zeros(size(nrem_est));
if numel(nrem_est)>0
    for iw = 1 : numel(nrem_est)
      time_in = (nrem_est(iw)-1)*winsize*fs+1:min(nspl,nrem_est(iw)*winsize*fs+1);
      high_sig(iw) = mean(teo(time_in));
    end  
else 
    high_sig = mean(teo);
end   
Lim_act = mean(high_sig);    
Lim_inf = 0.1*Lim_act;

% Energy condition (Jump with a minimal amplitude)
denv = diff([0; teo]);
pot_min = (denv(1:end-1)<0 & denv(2:end)>0 & teo(1:end-1)<Lim_inf);
keypoint = [1; find(pot_min); nspl];
activity = teo>Lim_act;
deb_act = find(diff([0; activity])==1);
fin_act = find(diff([activity; 0])==-1);
chck1 = (fin_act-deb_act)./fs>tr_tinf;
deb_act = deb_act(chck1);
fin_act = fin_act(chck1);
deb_epo = zeros(size(deb_act));
fin_epo = zeros(size(deb_act));
for iact = 1 : numel(deb_act)
  deb_epo(iact) = keypoint(find(keypoint<=deb_act(iact),1,'last'));
  fin_epo(iact) = keypoint(find(keypoint>=fin_act(iact),1));
end      
burst(:,1) = unique(deb_epo(deb_epo>fs & fin_epo<nspl-fs & (fin_epo(:)-deb_epo(:))/fs>=tr_tinf & (fin_epo(:)-deb_epo(:))/fs<=tr_tsup));
burst(:,2) = unique(fin_epo(deb_epo>fs & fin_epo<nspl-fs & (fin_epo(:)-deb_epo(:))/fs>=tr_tinf & (fin_epo(:)-deb_epo(:))/fs<=tr_tsup)); 

% -------------------------------------------------------------------------
%               SUBFUNCTIONS
%           _____________________
% -------------------------------------------------------------------------
function mx = MA(x,win)   
    mx = x;
    if win == 1
        mx(1:end-1) = (mx(2:end) + x(1:end-1))/2;
    else
        v = -round(win/2)+1:round(win/2);
        for ix = round(win/2) : numel(mx)-round(win/2)
            mx(ix) = mean(x(ix+v));
        end 
end 