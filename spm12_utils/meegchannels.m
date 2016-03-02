function ind = meegchannels(this, modality)
% Return indices of M/EEG channels
% 
% This function is inspired from SPM8 and was added to FASST to make it
% compatible with SPM12!
% Therefore it relies on the SPM12 @meeg object methods 'indchantype'.
%__________________________________________________________________
% Copyright (C) 2016 Cyclotron Research Centre, ULg, Belgium

% Written by C. Phillips

type = chantype(this); % all chantypes in the data

if nargin == 1
    ind = find(ismember(upper(type), ...
        {'EEG', 'MEG', 'MEGMAG', 'MEGGRAD', 'MEGPLANAR', 'REF', 'REFMAG', 'REFGRAD', 'LFP'}));
else
    switch modality
        case 'EEG'
            ind = find(ismember(upper(type), {'EEG'}));
        case 'LFP'
            ind = find(ismember(upper(type), {'LFP'}));
        case 'MEG'
            ind = find(ismember(upper(type), {'MEG', 'MEGMAG', 'MEGGRAD'}));
        case 'MEGMAG'
            ind = find(ismember(upper(type), {'MEGMAG'}));
        case 'MEGPLANAR'
            ind = find(ismember(upper(type), {'MEGPLANAR'}));
        case 'MEGGRAD' 
            ind = find(ismember(upper(type), {'MEGGRAD'})); 
        case 'MEEG'
            ind = find(ismember(upper(type), {'EEG', 'MEG', 'MEGMAG', 'MEGGRAD', 'MEGPLANAR'}));
        otherwise
            error('Unsupported modality.');
    end
end

ind = ind(:)'; % must be row to allow to use it as loop indices

return

%% ORGINAL SPM8 METHOD FOR THE MEEG OBJECT

% % FORMAT ind = meegchannels(this, modality)
% %
% %  this      - MEEG object
% %  modality  - one of EEG, LFP, MEG (excluding planar), MEGPLANAR, MEEG [optional]
% %
% %  ind       - row vector of M/EEG channels
% %
% % See also eogchannels, ecgchannels, emgchannels
% %__________________________________________________________________________
% % Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
% 
% % Stefan Kiebel
% % $Id: meegchannels.m 5045 2012-11-09 15:53:23Z vladimir $
% 
% 
% type = chantype(this);
% 
% if nargin == 1
%     ind = find(ismember(upper(type), ...
%         {'EEG', 'MEG', 'MEGMAG', 'MEGGRAD', 'MEGPLANAR', 'REF', 'REFMAG', 'REFGRAD', 'LFP'}));
% else
%     switch modality
%         case 'EEG'
%             ind = find(ismember(upper(type), {'EEG'}));
%         case 'LFP'
%             ind = find(ismember(upper(type), {'LFP'}));
%         case 'MEG'
%             ind = find(ismember(upper(type), {'MEG', 'MEGMAG', 'MEGGRAD'}));
%         case 'MEGMAG'
%             ind = find(ismember(upper(type), {'MEGMAG'}));
%         case 'MEGPLANAR'
%             ind = find(ismember(upper(type), {'MEGPLANAR'}));
%         case 'MEGGRAD' 
%             ind = find(ismember(upper(type), {'MEGGRAD'})); 
%         case 'MEEG'
%             ind = find(ismember(upper(type), {'EEG', 'MEG', 'MEGMAG', 'MEGGRAD', 'MEGPLANAR'}));
%         otherwise
%             error('Unsupported modality.');
%     end
% end
% 
% ind = ind(:)'; % must be row to allow to use it as loop indices
