function ind = emgchannels(this)
% Return indices of EMG channels
% FORMAT ind = emgchannels(this)
% 
% This function is inspired from SPM8 and was added to FASST to make it
% compatible with SPM12!
% Therefore it relies on the SPM12 @meeg object methods 'indchantype'.
%__________________________________________________________________
% Copyright (C) 2016 Cyclotron Research Centre, ULg, Belgium

% Written by C. Phillips

type = chantype(this);
ind = find(ismember(upper(type), {'EMG'}));
ind = ind(:)'; % must be row to allow to use it as loop indices

return

%% ORGINAL SPM8 METHOD FOR THE MEEG OBJECT

% function ind = emgchannels(this)
% % Return indices of EMG channels
% % FORMAT ind = emgchannels(this)
% %
% %  this      - MEEG object
% %  ind       - row vector of indices of EMG channels
% %
% % See also eogchannels, ecgchannels, meegchannels
% %__________________________________________________________________________
% % Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
% 
% % Christophe Phillips & Stefan Kiebel
% % $Id: emgchannels.m 2884 2009-03-16 18:27:25Z guillaume $
% 
% type = chantype(this);
% ind = find(ismember(upper(type), {'EMG'}));
% ind = ind(:)'; % must be row to allow to use it as loop indices