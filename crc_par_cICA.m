function Do = crc_par_cICA(Di)
% FORMAT Do = crc_par_cICA(Di)
%
% This function estimates the correction matrix by cICA, then picks the
% "best" one using an automatic mutual information criteria, and finally
% apllies the picked correction matrix onto the whole data set.
%__________________________________________________________________
% Copyright (C) 2011 Cyclotron Research Centre

% Written by C. Phillips
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

% Determine correction matrices, block by block
Di = crc_par_cICAmx(Di);
Di.cache = [];
save(Di);

% Find optimal correction matrix
ch_Mi = crc_par_cICAqa(Di);

% Creating the corrected data with selected correction matrix
Do = crc_par_cICAsav(Di,ch_Mi);

Di.CRC.cICA_ii = ch_Mi;
save(Di)

return