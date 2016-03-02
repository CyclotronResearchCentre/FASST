function sc_fasst = crc_convert_AseegaScore(sc_aseega)
% sc_fasst = crc_convert_AseegaScore(sc_aseega)
%
% Function that converts the score from Assega into the score from FASST
% because they use different conventions.
% Here only score 'hypno6' (following R&K) is considered.
%__________________________________________________________________________
% Copyright (C) 2014 Cyclotron Research Centre

% Written by C. Phillips, 2013.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

% Add Aseega scores as FASST scores
%  1 = stage 4      -> 4
%  3 = stage 3      -> 3
%  5 = stage 2      -> 2
%  7 = stage REM    -> 5
%  8 = stage 1      -> 1
%  10 = Wake        -> 0
%  11 = Artifact    -> 6 (MT)
% With 7 = unscorable windows (merged with MT for comparison)

sc_fasst = sc_aseega;
sc_fasst(sc_aseega==1) = 4;
sc_fasst(sc_aseega==5) = 2;
sc_fasst(sc_aseega==7) = 5;
sc_fasst(sc_aseega==8) = 1;
sc_fasst(sc_aseega==10) = 0;
sc_fasst(sc_aseega==11) = 6;

return
