function fmode = crc_ARdecomposition(sigma_event,fs)%D,eeg_chan)


% Automatic Sleep spindle_ARMAs detection method:
% Main steps:
% * Modelization ARMA to estimate the oscillatory mode 

% -------------------------------------------------------------------------
% ----------------------------- Parameters --------------------------
%--------------------------------------------------------------------------
ar_order = 8;       % autoregressive model order
ma_order = 2;       % moving average model order
fdelta = 4;
fbeta = 18;

% -------------------------------------------------------------------------
% ----------------------------- Processing  --------------------------
%--------------------------------------------------------------------------
 % -------------- ANALYSE de la FREQUENCE ----------------
 % calcul de la puissance instantanée (ARMA model)  
 sigma_event(:) = sigma_event(:)-mean(sigma_event(:));
 idevent = iddata(sigma_event(:),[],1/fs);
 try 
     thm = armax(idevent,[ar_order ma_order]);
 catch
     ar_order = 5;
     thm = armax(idevent,[ar_order ma_order]);
 end
 ar = thm.A; 
 fqc = DC_ARfeatures(ar,fs);   
 inclu = fqc.oscillatory & fqc.frequency>fdelta & fqc.frequency<fbeta;
 [gum idmax] = max(fqc.magnitude(inclu));
 mode = fqc.frequency(inclu);
fmode = mode(idmax);