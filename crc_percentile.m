function v = crc_percentile(x,p)
% Estimate the percentile value for a vector of values
% No frill, no check!
%__________________________________________________________________
% Copyright (C) 2010 Cyclotron Research Centre

% Written by C. Phillips, 2010.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id:$


if size(x,1)==1
    x = x';
end

x = sort(x);
n = length(x);
q = [0 100*(0.5:(n-0.5))./n 100]';
xx = [x(1); x(1:n); x(n)];
v = interp1q(q,xx,p);

return



