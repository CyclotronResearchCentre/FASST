function [time str] = crc_time_converts(secondes)
% FORMAT [time str] = crc_time_converts(secondes)
% Turn duration expressed in seconds into hours, minutes and seconds
%
% Input
% seconds - duration to convert
%
% Output
% time    - 1x3 vector [hour, min, sec]
% str     - same but expressed as a string
%__________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by Y. Leclercq & C. Phillips, 2008.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

if size(secondes,1)<size(secondes,2)
    secondes=secondes';
end

h=round(floor(secondes/60^2));
reste=mod(secondes,60^2);

m=round(floor(reste/60));
s=floor(mod(reste,60));

sms=s+(secondes-h*3600-m*60-s);

time=[h m sms];

dim=length(secondes);
bla=blanks(dim)';

hou=char(ones(1,dim)'*'h');
min=char(ones(1,dim)'*'m');
sec=char(ones(1,dim)'*'s');

str=[num2str(h) hou bla num2str(m) min bla num2str(s) sec];