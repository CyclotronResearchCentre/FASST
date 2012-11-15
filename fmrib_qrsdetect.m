% fmrib_qrsdetect() - Detect QRS peaks from ECG channel using combined
%   adaptive thresholding [Christov04];
%
% This program detects QRS peaks from a ECG channel.  First a complex lead
%   is constructed by computing the Teager Energy Operator TEO [Kim04] and 
%   then computing an adaptive threshold for each sampling point using a 
%   slightly modified version of [Christov04]. Points passing the threshold 
%   are counted as peaks.  A correction algorithm is then run (qrscorrect.m) 
%   which corrects for false positives and negatives and align the  peaks 
%   using correlation of the original ECG data.
%
% Usage:
%   >> peaks=fmrib_qrsdetect(EEG,ecgchan)
%
% Inputs:
%   EEG: EEGLAB data structure
%   ecgchan: The number of the ECG channel in the EEG data structure
%
% Ouptut:
%   peaks: index of QRS peak locations
%
%
% [Christov04] Real time electrocardiogram QRS detection using combined 
%  adaptive threshold, Ivaylo I. Christov.  Biomedical Engineering Online, 
%  BioMed Central (2004). available at:
%  http://www.biomedical-engineering-online.com/content/3/1/28
%
% [Kim04] Improved ballistocardiac artifact removal from the 
%  electroencephalogram recored in FMRI, KH Kim, HW Yoon, HW Park. 
% J NeouroSience Methods 135 (2004) 193-203.
%
%
% Author: Rami Niazy, FMRIB Centre, University of Oxford.
%   
% Copyright (c) 2004 University of Oxford

%123456789012345678901234567890123456789012345678901234567890123456789012
%
% Copyright (C) 2004 University of Oxford
% Author:   Rami K. Niazy, FMRIB Centre
%           rami@fmrib.ox.ac.uk
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

%   JUNE 03, 2005
%   Released after testing

%   APR 11, 2005
%   Beta version with more adaptive F threshold

%   APR 6, 2005
%   Removed accidentaly pasted code (BIG BUG)

%  MAR 16, 2005
%   Fixed typos and Misc Bugs
%   error message for NaN Decimation results

%   FEB 6, 2005
%   optimised k selection

%   DEC 23, 2004
%   Update (c)

%   Dec 15, 2004
%   Handles original fs
%   multiple of 125 Hz

%   Nov 11, 2004
%   Fixed R init

%   Nov 3, 2004
%   Corrected decimation
%   problem of new fs < 100

%   Oct 26, 2004
%   Fixed bug in calculating rem
%   Fixed bug when non-int fs used

% July 2005
% incorporated to pm_sleepEEG, including check display

% November 2006
% Some little corrections by cp, 
% WARNING, use only 1 EKG channel !!!!

% December 2006
% Added the possibility to pass *2* ecg channels, in this case, use the
% difference between the 2 specified channels for QRS detection.
% Added by cp

% November 2008
% Made compatible with spm8 dataformat.
% Added by yl

% March 2008
% Replacing access to structure meeg by object.
% cp

function Peaks=fmrib_qrsdetect(EEG,ecgchan,display_check)

nargchk(2,2,nargin);
dLFlag=1;
if ~exist('decimate.m','file')
    error('QRS detection requires the DSP toolbox');
end

ofs = EEG.fsample;

datachans = EEG.nchannels;

if any(ecgchan>datachans) || any(ecgchan<1)
    error('ECG channel out of data range, fmrib_qrsdetect() error!');
end
if length(ecgchan)==1
    ECG=EEG(ecgchan,:);
    dif_ecg = 0;
elseif length(ecgchan)==2
    ECG=EEG(ecgchan(1),:)-EEG(ecgchan(2),:);
    dif_ecg = 1;
else
    error('ECG channel defined as single channel or difference or 2 channels!');
end

ECG=double(ECG);

if rem(ofs,128)==0
    dL=ofs/128;    
elseif rem(ofs,100)==0
    dL=ofs/100;
elseif rem(ofs,125)==0
    dL=ofs/125;  
else
    dL=round(ofs/100);
    if ofs/dL < 100
        dLFlag=0;
        dL=1;
    end
end


%Decimate signal
%---------------
if dLFlag
    if dL>4
        if rem(dL,2)==0
            Ecg=decimate(ECG,dL/2);
            Ecg=decimate(Ecg,2);
        elseif rem(dL,3)==0
            Ecg=decimate(ECG,dL/3);
            Ecg=decimate(Ecg,3);
        elseif rem(dL,5)==0
            Ecg=decimate(ECG,dL/5);
            Ecg=decimate(Ecg,5);
        elseif rem(dL,7)==0
            Ecg=decimate(ECG,dL/7);
            Ecg=decimate(Ecg,7);
        elseif rem(dL,9)==0
            Ecg=decimate(ECG,dL/9);
            Ecg=decimate(Ecg,9);
        else
            try 
                Ecg=decimate(ECG,dL);
            catch
                Ecg=ECG;
                dL=1;
%                 dLFlag=0;
            end
        end
    else
        Ecg=decimate(ECG,dL);
    end
else
    Ecg=ECG;
end
fs=ofs/dL;
if find(isnan(Ecg)==1)
    error('Decimation failed.  Downsample the data first and try again');
end



%MFR Settings and init
%----------------------
L=length(Ecg);
msWait=floor(0.55*fs);
ms1200=floor(1.2*fs);
ms350=floor(0.35*fs);
ms300=floor(0.3*fs);
ms50=floor(0.05*fs);
Mc=0.45;
s5=floor(5*fs);
DetectFlag=0;
timer1=0;
Peaks=[];
peakc=1;
firstdetect=1;
Ecg=Ecg(:);

%Allocate memory
%---------------
M=zeros(L,1);
R=zeros(L,1);
F=zeros(L,1);
MFR=zeros(L,1);
Y=zeros(L,1);
M5=ones(5,1);
R5=ones(5,1);
% F350=zeros(ms350,1);

%Pre-proc Filtering
%-----------------
fL  = round(fs/50);
b   = ones(1,fL)/fL;
Ecg = filtfilt(b,1,Ecg);
fL  = round(fs/35);
b   = ones(1,fL)/fL;
Ecg = filtfilt(b,1,Ecg);

%Estimate init R and k
%-----------------
FFTp=round(100*fs);
P2=ceil(log(FFTp)/log(2));
NFFT=2^P2;
Fecg=fft(detrend(Ecg(1:round(5*fs))).*hann(length(Ecg(1:round(5*fs)))),NFFT);
Pecg=Fecg.*conj(Fecg) / NFFT;
[MV,ML]=max(Pecg);
R5=R5*round(NFFT/ML);
k= round(fs*fs*pi/(2*2*pi*10*(R5(1))));


%Construct complex lead Y using TEO
%----------------------------------
f = [0 7/(fs/2) 9/(fs/2) 40/(fs/2) 42/(fs/2) 1];
a = [0 0 1 1 0 0];
wts = firls(100,f,a);
% ecgF=Ecg;
ecgF = filtfilt(wts,1,Ecg);
for n=(k+1):(L-k)
    Y(n)=ecgF(n)^2-ecgF(n-k)*ecgF(n+k);
end
Y(L)=0;
fL=round(fs/25);
b=ones(1,fL)/fL;
Y=filtfilt(b,1,Y);
Y(Y<0)=0;


%init M and F
%-------------
M5=Mc*max(Y(round(fs):round(fs+s5)))*M5;
M(1:s5)=mean(M5);
newM5=mean(M5);
F(1:ms350)=mean(Y(fs:fs+ms350));
F2(1:ms350)=F(1:ms350);

%Detect QRS
%----------
%wait bar
barth=5;
barth_step=barth;
Flag25=0; Flag50=0; Flag75=0;
fprintf('\nStage 1 of 5: Adaptive threshold peak detection.\n');

for n=1:L
       
    %-----------------calc------------------------------------------
    timer1=timer1+1;
        
    if length(Peaks)>=2
        if DetectFlag==1
            DetectFlag=0;
            M(n)=mean(M5);
            Mdec=(M(n)-M(n)*Mc)/(ms1200-msWait);
            Rdec=Mdec/1.4;
        elseif DetectFlag==0 && (timer1<=msWait || timer1 > ms1200)
            M(n)=M(n-1);
        elseif DetectFlag==0 && timer1 == msWait+1
            M(n)=M(n-1)-Mdec;
            newM5=Mc*max(Y(n-msWait:n));
            if newM5 > 1.5*M5(5)
                newM5=1.5*M5(5);
            end
            M5=[M5(2:end);newM5];
        elseif DetectFlag==0 && timer1 > msWait+1 && timer1 <= ms1200
            M(n)=M(n-1)-Mdec;
        end
    end
    
    if n>ms350
        F(n)=F(n-1)+(max(Y(n-ms50+1:n))-max(Y(n-ms350+1:n-ms300)))/150;
        F2(n)=F(n)-mean(Y(fs:fs+ms350))+newM5;
        %F(n)=mean(Y(n-ms350+1:n))+(max(Y(n-ms50+1:n))-...
         %   max(Y(n-ms350+1:n-ms300)))/150;
%          Me(n)=mean(Y(n-ms350+1:n));
%          Fl(n)=max(Y(n-ms50+1:n));
%          Fe(n)=max(Y(n-ms350+1:n-ms300));
    end            
    
    Rm=mean(R5);
    R0int=round(2*Rm/3);
    if timer1 <= R0int
        R(n)=0;
    elseif length(Peaks) >=2
        R(n)=R(n-1)-Rdec;
    end
        
    MFR(n)=M(n)+F2(n)+R(n);
    
        
    if (all(Y(n)>=MFR(n)) && timer1 > msWait) || (all(Y(n)>=MFR) && firstdetect==1)
        if firstdetect==1;
            firstdetect=0;
        end
        Peaks(peakc)=n;
        if peakc>1
            R5=[R5(2:end);Peaks(peakc)-Peaks(peakc-1)];
        end
        peakc=peakc+1;
        DetectFlag=1;
        timer1=-1;
    end
    
    %---------------------------------------------------------------
    %update wait bar
    percentdone=floor(n*100/L);
    if floor(percentdone)>=barth
        if percentdone>=25 && Flag25==0
            fprintf('25%% ')
            Flag25=1;
        elseif percentdone>=50 && Flag50==0
            fprintf('50%% ')
            Flag50=1;
        elseif percentdone>=75 && Flag75==0
            fprintf('75%% ')
            Flag75=1;
        elseif percentdone==100
            fprintf('100%%\n')
        else
            fprintf('.')
        end

        while barth<=percentdone
            barth=barth+barth_step;
        end
        if barth>100
            barth=100;
        end
    end 
end

%correct QRS Peaks
%-----------------
Peaks=qrscorrect(Peaks,Ecg,fs);
Peaks=Peaks*dL;

% Check display
if display_check == 1
    
    Fmenu = figure('IntegerHandle','on',...
        'Name',sprintf('%s%s','VISUAL CHECK of QRS DETECTION'),...
        'NumberTitle','off',...
        'Tag','Menu',...
        'Position',[275 476 1000 500],...
        'Resize','off',...
        'Color',[.5 .5 1],...
        'MenuBar','none',...
        'DefaultTextFontName','Helvetica',...
        'DefaultTextFontSize',10,...
        'DefaultUicontrolFontName','Helvetica',...
        'DefaultUicontrolFontSize',12,...
        'DefaultUicontrolInterruptible','on',...
        'Renderer','zbuffer',...
        'Visible','on');
    
    uicontrol(Fmenu,'Style','pushbutton',...
        'String','Done',...
        'Position',[920 150 50 200],...
        'FontSize',10,...
        'ToolTipString','quiting qrs_detect',...
        'UserData','Exit',...
        'CallBack','close(gcf);return;',...
        'Visible','on');                                                                              
    varargout = {Fmenu};
    
    set(Fmenu,'Visible','on')
    
    if dif_ecg
%         Ylim = [min(EEG.data.y(ecgchan(1),:)-EEG.data.y(ecgchan(2),:)) ...
%                 max(EEG.data.y(ecgchan(1),:)-EEG.data.y(ecgchan(2),:))] ;
        Ylim = [min(EEG(ecgchan(1),:)-EEG(ecgchan(2),:)) ...
                max(EEG(ecgchan(1),:)-EEG(ecgchan(2),:))] ;
    else
%         Ylim = [min(EEG.data.y(ecgchan,:)) max(EEG.data.y(ecgchan,:))];
        Ylim = [min(EEG(ecgchan,:)) max(EEG(ecgchan,:))];
    end
%     t = 0:1/EEG.Fsample:EEG.Nsamples/EEG.Fsample-1/EEG.Fsample;
%     N = zeros(1,EEG.Nsamples);
    t = 0:1/EEG.fsample:EEG.nsamples/EEG.fsample-1/EEG.fsample;
    N = zeros(1,EEG.nsamples);
    N(Peaks) = Ylim(2);
    hold
%     plot(t,-EEG.data.y(length(EEG.channels),:))
    plot(t,-EEG(EEG.nchannels,:))
    plot(t,N,'-r')
    grid on
    xlabel('time (s)')
    ylabel('Amplitude (microV)')
%     uiscroll([1 EEG.Nsamples],10,[],[],Ylim);
    uiscroll([1 EEG.nsamples],10,[],[],Ylim);
end

return;
