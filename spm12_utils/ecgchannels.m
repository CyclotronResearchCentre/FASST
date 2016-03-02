function ind = ecgchannels(this)
% Return indices of ECG channels
% FORMAT ind = ecgchannels(this)
%
% This function is inspired from SPM8 and was added to FASST to make it
% compatible with SPM12!
% Therefore it relies on the SPM12 @meeg object methods 'indchantype'.
%__________________________________________________________________
% Copyright (C) 2016 Cyclotron Research Centre, ULg, Belgium

% Written by C. Phillips

type = chantype(this);
ind = find(ismember(upper(type), {'ECG', 'EKG'}));
ind = ind(:)'; % must be row to allow to use it as loop indices

return

%% ORGINAL SPM8 METHOD FOR THE MEEG OBJECT

% function ind = ecgchannels(this)
% % Return indices of ECG channels
% % FORMAT ind = ecgchannels(this)
% %
% %  this      - MEEG object
% %  ind       - row vector of indices of ECG channels
% %
% % See also eogchannels, emgchannels, meegchannels
% %__________________________________________________________________________
% % Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
% 
% % Christophe Phillips & Stefan Kiebel
% % $Id: ecgchannels.m 2884 2009-03-16 18:27:25Z guillaume $
% 
% type = chantype(this);
% ind = find(ismember(upper(type), {'ECG', 'EKG'}));
% ind = ind(:)'; % must be row to allow to use it as loop indices
