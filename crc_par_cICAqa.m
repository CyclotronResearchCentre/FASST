function [ch_Mi]=crc_par_cICAqa(Di)

% Picking the "best" correction matrix using a 'mutual information' (MI)
% criteria.
%__________________________________________________________________
% Copyright (C) 2011 Cyclotron Research Centre

% Written by Y. Leclercq & C. Phillips, 2010.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$


% check if MInfo directory is in the path, add it if necessary
if ~exist('mutualinfo.m','file')
    Pfasst = which('crc_main');
    P_MInfo = fullfile(spm_str_manip(Pfasst,'h'),'MInfo');
    addpath(P_MInfo);
end

crcdef = crc_get_defaults('par');

if nargin < 1
    fname = spm_select(1,'mat', 'Select imported EEG file','' ,pwd);
    Di = crc_eeg_load(fname);
end

fs = Di.fsample;

lgts=crcdef.length;
ideblist=[];
ifinlist=[];
Xlist=[];
Jlist=[];
if or(crcdef.useinitseg,and(~crcdef.useinitseg,~crcdef.additioseg))
    ideblist=[ideblist 1];
    ifinlist=[ifinlist lgts*fs];
end

if crcdef.additioseg
    newnbrs=rand(1,crcdef.additioseg);
    sznbrs=Di.nsamples-(2*lgts*fs);
    ideblist=[ideblist round(newnbrs*sznbrs)+lgts*fs];
    ifinlist=[ifinlist round(newnbrs*sznbrs)+2*lgts*fs];
end

mxlist=[];

Refnames = crcdef.Refnames ;

badindex = Di.CRC.cICA{4};
selected_chan = setdiff(1:Di.nchannels, badindex);
selnames = chanlabels(Di,selected_chan);
refind = [];
for kk = 1:length(Refnames)
    search = find(strcmp(upper(Refnames(kk)),upper(selnames)));
    if ~isempty(search)
        refind = [refind search];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for zz=1:length(ideblist)

    ideb=ideblist(zz);
    ifin=ifinlist(zz);

    data_o = Di(:,max(1,ideb):min(ifin,Di.nsamples));

    %%%%%%%%%Filter the very low frequencies

    zx=1.5*fs;
    zx=ones(1,zx)/zx;
    data2=conv2(data_o,zx,'same');
    data=data_o-data2;
    clear data2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Gives the peaks of ECG and BCG template for each channel
    allPeaks = Di.CRC.EKGPeaks;
    Peaks = allPeaks(and(ideb<allPeaks,allPeaks<ifin)==1)-ideb+1;
    data_oref = data_o(refind,:);

    %Check the peaks
    lim = mean(diff(Peaks))/3;
    list = find(diff(Peaks)>lim);

    Peaks = Peaks([list length(Peaks)]);

    fitted_art = mean_art(data_oref,Peaks);
    plotdat = data(refind,:);

    mxmui=zeros(size(Di.CRC.cICA{3},2)+1,length(refind));

    for ii=1:size(plotdat,1)
        mxmui(1,ii)=mutualinfo(fitted_art(ii,:),plotdat(ii,:));
    end

    for jj=1:length(Di.CRC.cICA{3})

        cormx = Di.CRC.cICA{3}{jj};
        plotdat = data;
        plotdat(selected_chan,:) = cormx*data(selected_chan,:);
        plotdat = plotdat(refind,:);

        for ii=1:size(plotdat,1)
            mxmui(jj+1,ii)=mutualinfo(fitted_art(ii,:),plotdat(ii,:));
        end

    end
    mxlist = [mxlist mxmui];
    [X, J] = min(mean(abs(mxlist),2));
    currentmui = J-1;
    Xlist = [Xlist X];
    Jlist = [Jlist currentmui];
end

[X, J] = min(mean(abs(mxlist),2));
[Y, I] = min(Xlist);

perc = mxlist;

for uu = 2:size(perc,1)
    perc(uu,:)=perc(uu,:)./perc(1,:);
end
perc(1,:)=perc(1,:)./perc(1,:);

% QUESTION:
% what did Yves with this choice of the optimal correction matrix ?!
ch_Mi = J-1;  % Matrix that reduces the mutual information in mean
ch_Mi = Jlist(I);    % Matrix that reaches the minimal mutual information on one of the segments
[minperc, minidx]=min(mean(perc,2));
ch_Mi = minidx-1;


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the mean artefact of the pulse artefact given ECG peaks.

function outdat=mean_art(d,pks)

Nit = min(15,length(pks)-1);

step = max(diff(pks));

outdat = zeros(size(d));
warhouse = zeros(size(d,1),(length(pks)-1)*step);

for ii = 1:length(pks) - 1

    warhouse(:,1+step*(ii-1):step*ii) = ...
        remap(d(:,pks(ii):pks(ii+1)-1),step);

end

Cube=reshape(warhouse(:,1:step*Nit),[size(warhouse,1) step Nit]);
iCube = 1;

for ii = 1:length(pks) - 1

    Cube(:,:,iCube) = warhouse(:,1+step*(ii-1):step*ii);
    warhouse(:,1+step*(ii-1):step*ii) = mean(Cube,3);

    avg = mean(warhouse(:,1+step*(ii-1):step*ii),2)*ones(1,step);

    warhouse(:,1+step*(ii-1):step*ii) = warhouse(:,1+step*(ii-1):step*ii) - avg;

    outdat(:,pks(ii):pks(ii+1)-1) = ...
        remap(warhouse(:,1+step*(ii-1):step*ii),pks(ii+1)-pks(ii));

    iCube = mod(iCube+1,20);
    if iCube == 0
        iCube = 20;
    end

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remap the data to the specified number of points.

function outdat=remap(data,step)

outdat = zeros(size(data,1),step);

sz = size(data,2)-1;
step=step-1;
for i=1:step

    w = mod((i-1)*sz/(step),1);
    idx = 1+floor((i-1)*sz/(step));
    outdat(:,i)= data(:,idx)*(1-w) + data(:,idx+1)*(w);
end

outdat(:,end)=data(:,end);

end