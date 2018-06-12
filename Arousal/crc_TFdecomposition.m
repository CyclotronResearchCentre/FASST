function [flags tempo fmode_mean fmode_arma intrafreq Pmax intrapower] = crc_TFdecomposition(event,fs,teo)
            
% -------------------------------------------------------------------------
% ----------------------------- Parameters --------------------------
%--------------------------------------------------------------------------
flags = 0;
tempo = [0 0];
intrafreq = {[]};
intrapower = {[]};
fmode_mean =0;
fmode_arma = 0;
Pmax = 0;
fmin = 8;
fmax = 30;
tr_inf = 0.4;
prct = prctile(teo(fs:end-fs),70);
MainPart = find(teo(fs:end-fs)>prct)+fs;
% -------------------------------------------------------------------------
% ----------------------------- ASDM --------------------------
%--------------------------------------------------------------------------  

[Tx , Fs] = sswt(event,fs,'fmin',fmin,'fmax',fmax,'display','off');  
ResFq = max(diff(Fs));
Ex = abs(Tx);
[TPower frq] = max(Ex);

% check the interruption
spi = 1;
it = 1;
intraspi = cell(1);
while it <= numel(frq)
    [intraspi{spi} it] = checknext(frq,it,ceil(0.8/ResFq),0.1*fs,Fs);
    spi = spi + 1;
end
debspi = 1;
for iwave = 1 : size(intraspi,2)
    endspi = debspi + numel(intraspi{iwave})-1;
    durin = endspi-debspi+1;
    if durin/fs>tr_inf && numel(intersect(debspi:endspi,MainPart))>0.70*numel(MainPart) && ~or(endspi<fs+0.5*fs,debspi>size(Ex,2)-fs-0.5*fs)
        fmode = crc_ARdecomposition(event(debspi:endspi),fs);  
        if abs((fmode-mean(intraspi{iwave}(intraspi{iwave}>0))))<1
            flags = flags+1;
            tempo(flags,:) = [debspi endspi]; 
            fmode_mean(flags) = mean(intraspi{iwave}(intraspi{iwave}>0));
            fmode_arma(flags) = fmode;
            intrafreq{flags} = intraspi{iwave};
            Pmax(flags) = max(TPower(debspi:endspi));
            intrapower{flags} = TPower(debspi:endspi);
        end
    end
    debspi = endspi +1;
end

% -------------------------------------------------------------------------
%               SUBFUNCTIONS
%           _____________________
% -------------------------------------------------------------------------
         
function [intraspi pointeur] = checknext(frq,it,Tf,Tt,Fs)
pointeur = it;
go = 1;
intraspi = [];
while and(pointeur+go<numel(frq), go <= Tt)
    if abs(frq(pointeur+go)-frq(pointeur))>Tf
        go = go+1;
    else
        intraspi = [intraspi zeros(1,go-1) Fs(frq(pointeur))];
        pointeur = pointeur+go;
        go = 1;
    end
end
pointeur = pointeur+1;