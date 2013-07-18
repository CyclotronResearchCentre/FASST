function Do = crc_par_cICAsav(Di,ch_Mi)
% FORMAT Do = crc_par_cICAsav(Di,ch_Mi)
%
% This function applies the selected correction matrix (index ch_Mi) on the
% data Di, and creates the cICA-corrected data set.
%__________________________________________________________________
% Copyright (C) 2011 Cyclotron Research Centre

% Written by C. Phillips
% Cyclotron Research Centre, University of Liege, Belgium
% $Id:$

% Load chosen correction matrix
if numel(ch_Mi)==1
    cormx = Di.CRC.cICA{3}{ch_Mi};
    one_corM = 1;
else
    N_corM = numel(ch_Mi);
    one_corM = 0;
end

% Prepare output data
crcdef = crc_get_defaults('par');
prefix = crcdef.cicapref;
fn_o = [prefix,Di.fname];
Do = clone(Di,fn_o);

% Size of the sample analyzed in one step in seconds (advice: 60 sec is the
% minimum.
dat_sz = (Di.nsamples-1)/Di.fsample; % in seconds
sz = min(Di.CRC.cICA{2},dat_sz);
% Step in seconds (advice: step = 2/3 size)
step = Di.CRC.cICA{1};
fs = Di.fsample;

badindex = sort(Di.CRC.cICA{4}); % gives the badindex
selected_chan = setdiff(1:Di.nchannels, badindex);

ideb = 1;
ind_it = 0;
ifin = min(floor(1+sz*fs),Di.nsamples);
while ideb<Di.nsamples
    if ~one_corM
        ind_it = ind_it+1;
        cormx = Di.CRC.cICA{3}{ch_Mi(min(ind_it,N_corM))};
    end
    cordata = Di(:,ideb:ifin); % load current data;
    cordata(selected_chan,:) = cormx*cordata(selected_chan,:); %Correct current data
    Do(:,ideb:ifin) = cordata;
    ideb = ideb+round(step)*fs;
    ifin = min(ifin+round(step)*fs,Do.nsamples);
end

return