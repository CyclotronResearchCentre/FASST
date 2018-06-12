function C = DC_extra(C)

% FORMAT DC_extra(args)
% It regroups all the artifacts detected in one vector:
%       D.CRC.DC.shortartf.total
% INPUT
%       .file   - data file (.mat files)
% 
% NB: it has to be processed after DC_preprocessing, DC_badchannels,
% DC_pwr, DC_popping, DC_emg and DC_reference
% 
% _________________________________________________________________
% Copyright (C) 2014 Cyclotron Research Centre

% Written by D. Coppieters 't Wallant, 2014.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$
% ----------------------

% if isempty(varargin)
%     C = spm_eeg_load; 
% else 
%     C = spm_eeg_load(varargin{1});
% end

% parametres généraux
arousal = C.CRC.DC.shortartf.arousal;
movement = C.CRC.DC.shortartf.movement;
popping = C.CRC.DC.shortartf.popping;
badchan = C.CRC.DC.shortartf.badchan;  
shortartef = unique([badchan(:); movement(:); popping(:); arousal(:)]);
while ~isempty(find(or(or(diff(shortartef)==2,diff(shortartef)==3),diff(shortartef)==4)))
    inter = diff(shortartef);
    inter_1 = inter == 2 ;
    inter_2 = inter == 3 ;
    inter_3 = inter == 4 ;
    to_be_add1 = shortartef(inter_1) + 1;
    to_be_add2 = shortartef(inter_2) + 2;
    to_be_add3 = shortartef(inter_3) + 3;
    shortartef = setdiff(unique([shortartef(:);to_be_add1(:); to_be_add2(:); to_be_add3(:)]),0);
    shortartef = unique(shortartef);
end
C.CRC.DC.shortartf.artefact = shortartef;

%%% to add artefact into scoring field
% artefact / arousal on all channel
tmp = diff(C.CRC.DC.shortartf.artefact);
ind=find(tmp > 1);
fin = [C.CRC.DC.shortartf.artefact(ind); C.CRC.DC.shortartf.artefact(end)];
debut = [C.CRC.DC.shortartf.artefact(1); C.CRC.DC.shortartf.artefact(ind+1)];
% fin = [C.CRC.DC.shortartf.artefact(ind)-1; C.CRC.DC.shortartf.artefact(end)-1];% to relign at the second
% debut = [C.CRC.DC.shortartf.artefact(1)-1; C.CRC.DC.shortartf.artefact(ind+1)-1];% to realign at th second
C.CRC.score{5} = [debut fin];% we place artefact and arousal on all channel in what is considered as artefact in FASST
% artefact on single channels
tmp2 = [C.CRC.DC.badchannels.chan_incoherent];%.chan_incoherent: noisy channel
for bcl = 1:size(C.CRC.DC.badchannels.chan_defaillant,1) %.chan_defaillant: flat channels
    if ~isempty(C.CRC.DC.badchannels.chan_defaillant{bcl,1}) == 1
        tmp2(bcl) = C.CRC.DC.badchannels.chan_defaillant(bcl,1); % tmp2 contains all artefacted channel (a channel cannot be noisy and flat at the same time)
    end
end
debut2 = [0:30:((size(tmp2,1)-1)*30)+1]';% start of each 30s window
fin2 = [30:30:((size(tmp2,1)-1)*30) nsamples(C)/fsample(C)]'; %end of each 30s window
tmp3 = [];
for ar = 1:size(tmp2,1)
    if ~isempty(tmp2{ar}) == 1
        for te = 1:size(tmp2{ar},2)
            tmp3 = [tmp3 ; debut2(ar) fin2(ar) tmp2{ar}(te)];% each bad channel has its onw start and end >> multiple line for windows for several bad channels
        end
    end
end

C.CRC.score{6} = tmp3;% we place artefact on single channels in what is considered as araousal in FASST
C.CRC.score = C.CRC.score';
save(C);
fprintf('all artefacts are in the vector D.CRC.DC.shortartf.artefact and arousals in D.CRC.DC.shortartf.arousal... %d \n', length(shortartef));