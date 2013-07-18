% Copyright (C) 1999 Paul Kienzle
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; If not, see <http://www.gnu.org/licenses/>.

% Generate a butterworth filter.
% Default is a discrete space (Z) filter.
%
% [b,a] = butter(n, Wc)
%    low pass filter with cutoff pi*Wc radians
%
% [b,a] = butter(n, Wc, 'high')
%    high pass filter with cutoff pi*Wc radians
%
% [b,a] = butter(n, [Wl, Wh])
%    band pass filter with edges pi*Wl and pi*Wh radians
%
% [b,a] = butter(n, [Wl, Wh], 'stop')
%    band reject filter with edges pi*Wl and pi*Wh radians
%
% [z,p,g] = butter(...)
%    return filter as zero-pole-gain rather than coefficients of the
%    numerator and denominator polynomials.
%
% [...] = butter(...,'s')
%     return a Laplace space filter, W can be larger than 1.
%
% [a,b,c,d] = butter(...)
%  return  state-space matrices
%
% References:
%
% Proakis & Manolakis (1992). Digital Signal Processing. New York:
% Macmillan Publishing Company.

% Author: Paul Kienzle <pkienzle@user.sf.net>
% Modified by: Doug Stewart <dastew@sympatico.ca> Feb, 2003

function [a, b, c, d] = butter(n, W, varargin)

if (nargin>4 || nargin<2) || (nargout>4 || nargout<2)
    usage ('[b, a] or [z, p, g] or [a,b,c,d] = butter (n, W [, "ftype"][,"s"])');
end

% interpret the input parameters
if (~(length(n)==1 && n == round(n) && n > 0))
    error ('butter: filter order n must be a positive integer');
end

stop = 0;
digital = 1;
for i=1:length(varargin)
    switch varargin{i}
        case 's', digital = 0;
        case 'z', digital = 1;
        case { 'high', 'stop' }, stop = 1;
        case { 'low',  'pass', 'bandpass' }, stop = 0;
        otherwise,  error ('butter: expected [high|stop] or [s|z]');
    end
end


[r, c]=size(W);
if (~(length(W)<=2 && (r==1 || c==1)))
    error ('butter: frequency must be given as w0 or [w0, w1]');
elseif (~(length(W)==1 || length(W) == 2))
    error ('butter: only one filter band allowed');
elseif (length(W)==2 && ~(W(1) < W(2)))
    error ('butter: first band edge must be smaller than second');
end

if ( digital && ~all(W >= 0 & W <= 1))
    error ('butter: critical frequencies must be in (0 1)');
elseif ( ~digital && ~all(W >= 0 ))
    error ('butter: critical frequencies must be in (0 inf)');
end

% Prewarp to the band edges to s plane
if digital
    T = 2;       % sampling frequency of 2 Hz
    W = 2/T*tan(pi*W/T);
end

% Generate splane poles for the prototype butterworth filter
% source: Kuc
C = 1; % default cutoff frequency
pole = C*exp(1i*pi*(2*[1:n] + n - 1)/(2*n));
if mod(n,2) == 1, pole((n+1)/2) = -1; end  % pure real value at exp(i*pi)
zero = [];
gain = C^n;

% splane frequency transform
[zero, pole, gain] = sftrans(zero, pole, gain, W, stop);

% Use bilinear transform to convert poles to the z plane
if digital
    [zero, pole, gain] = bilinear(zero, pole, gain, T);
end

% convert to the correct output form
if nargout==2,
    a = real(gain*poly(zero));
    b = real(poly(pole));
elseif nargout==3,
    a = zero;
    b = pole;
    c = gain;
else
    % output ss results
    [a, b, c, d] = zp2ss (zero, pole, gain);
end

return

%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Copyright (C) 1999 Paul Kienzle
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; If not, see <http://www.gnu.org/licenses/>.

% usage: [Zz, Zp, Zg] = bilinear(Sz, Sp, Sg, T)
%        [Zb, Za] = bilinear(Sb, Sa, T)
%
% Transform a s-plane filter specification into a z-plane
% specification. Filters can be specified in either zero-pole-gain or
% transfer function form. The input form does not have to match the
% output form. 1/T is the sampling frequency represented in the z plane.
%
% Note: this differs from the bilinear function in the signal processing
% toolbox, which uses 1/T rather than T.
%
% Theory: Given a piecewise flat filter design, you can transform it
% from the s-plane to the z-plane while maintaining the band edges by
% means of the bilinear transform.  This maps the left hand side of the
% s-plane into the interior of the unit circle.  The mapping is highly
% non-linear, so you must design your filter with band edges in the
% s-plane positioned at 2/T tan(w*T/2) so that they will be positioned
% at w after the bilinear transform is complete.
%
% The following table summarizes the transformation:
%
% +---------------+-----------------------+----------------------+
% | Transform     | Zero at x             | Pole at x            |
% |    H(S)       |   H(S) = S-x          |    H(S)=1/(S-x)      |
% +---------------+-----------------------+----------------------+
% |       2 z-1   | zero: (2+xT)/(2-xT)   | zero: -1             |
% |  S -> - ---   | pole: -1              | pole: (2+xT)/(2-xT)  |
% |       T z+1   | gain: (2-xT)/T        | gain: (2-xT)/T       |
% +---------------+-----------------------+----------------------+
%
% With tedious algebra, you can derive the above formulae yourself by
% substituting the transform for S into H(S)=S-x for a zero at x or
% H(S)=1/(S-x) for a pole at x, and converting the result into the
% form:
%
%    H(Z)=g prod(Z-Xi)/prod(Z-Xj)
%
% Please note that a pole and a zero at the same place exactly cancel.
% This is significant since the bilinear transform creates numerous
% extra poles and zeros, most of which cancel. Those which do not
% cancel have a 'fill-in' effect, extending the shorter of the sets to
% have the same number of as the longer of the sets of poles and zeros
% (or at least split the difference in the case of the band pass
% filter). There may be other opportunistic cancellations but I will
% not check for them.
%
% Also note that any pole on the unit circle or beyond will result in
% an unstable filter.  Because of cancellation, this will only happen
% if the number of poles is smaller than the number of zeros.  The
% analytic design methods all yield more poles than zeros, so this will
% not be a problem.
%
% References:
%
% Proakis & Manolakis (1992). Digital Signal Processing. New York:
% Macmillan Publishing Company.

% Author: Paul Kienzle <pkienzle@users.sf.net>

function [Zz, Zp, Zg] = bilinear(Sz, Sp, Sg, T)

if nargin==3
    T = Sg;
    [Sz, Sp, Sg] = tf2zp(Sz, Sp);
elseif nargin~=4
    usage('[Zz, Zp, Zg]=bilinear(Sz,Sp,Sg,T) or [Zb, Za]=blinear(Sb,Sa,T)');
end;

p = length(Sp);
z = length(Sz);
if z > p || p==0
    error('bilinear: must have at least as many poles as zeros in s-plane');
end

% ----------------  -------------------------  ------------------------
% Bilinear          zero: (2+xT)/(2-xT)        pole: (2+xT)/(2-xT)
%      2 z-1        pole: -1                   zero: -1
% S -> - ---        gain: (2-xT)/T             gain: (2-xT)/T
%      T z+1
% ----------------  -------------------------  ------------------------
Zg = real(Sg * prod((2-Sz*T)/T) / prod((2-Sp*T)/T));
Zp = (2+Sp*T)./(2-Sp*T);
if isempty(Sz)
    Zz = -ones(size(Zp));
else
    Zz = [(2+Sz*T)./(2-Sz*T)];
    Zz = postpad(Zz, p, -1);
end

if nargout==2, [Zz, Zp] = zp2tf(Zz, Zp, Zg); end

return

%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Copyright (C) 1999 Paul Kienzle
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; If not, see <http://www.gnu.org/licenses/>.

% usage: [Sz, Sp, Sg] = sftrans(Sz, Sp, Sg, W, stop)
%
% Transform band edges of a generic lowpass filter (cutoff at W=1)
% represented in splane zero-pole-gain form.  W is the edge of the
% target filter (or edges if band pass or band stop). Stop is true for
% high pass and band stop filters or false for low pass and band pass
% filters. Filter edges are specified in radians, from 0 to pi (the
% nyquist frequency).
%
% Theory: Given a low pass filter represented by poles and zeros in the
% splane, you can convert it to a low pass, high pass, band pass or
% band stop by transforming each of the poles and zeros individually.
% The following table summarizes the transformation:
%
% Transform         Zero at x                  Pole at x
% ----------------  -------------------------  ------------------------
% Low Pass          zero: Fc x/C               pole: Fc x/C
% S -> C S/Fc       gain: C/Fc                 gain: Fc/C
% ----------------  -------------------------  ------------------------
% High Pass         zero: Fc C/x               pole: Fc C/x
% S -> C Fc/S       pole: 0                    zero: 0
%                   gain: -x                   gain: -1/x
% ----------------  -------------------------  ------------------------
% Band Pass         zero: b ? sqrt(b^2-FhFl)   pole: b ? sqrt(b^2-FhFl)
%        S^2+FhFl   pole: 0                    zero: 0
% S -> C --------   gain: C/(Fh-Fl)            gain: (Fh-Fl)/C
%        S(Fh-Fl)   b=x/C (Fh-Fl)/2            b=x/C (Fh-Fl)/2
% ----------------  -------------------------  ------------------------
% Band Stop         zero: b ? sqrt(b^2-FhFl)   pole: b ? sqrt(b^2-FhFl)
%        S(Fh-Fl)   pole: ?sqrt(-FhFl)         zero: ?sqrt(-FhFl)
% S -> C --------   gain: -x                   gain: -1/x
%        S^2+FhFl   b=C/x (Fh-Fl)/2            b=C/x (Fh-Fl)/2
% ----------------  -------------------------  ------------------------
% Bilinear          zero: (2+xT)/(2-xT)        pole: (2+xT)/(2-xT)
%      2 z-1        pole: -1                   zero: -1
% S -> - ---        gain: (2-xT)/T             gain: (2-xT)/T
%      T z+1
% ----------------  -------------------------  ------------------------
%
% where C is the cutoff frequency of the initial lowpass filter, Fc is
% the edge of the target low/high pass filter and [Fl,Fh] are the edges
% of the target band pass/stop filter.  With abundant tedious algebra,
% you can derive the above formulae yourself by substituting the
% transform for S into H(S)=S-x for a zero at x or H(S)=1/(S-x) for a
% pole at x, and converting the result into the form:
%
%    H(S)=g prod(S-Xi)/prod(S-Xj)
%
% The transforms are from the references.  The actual pole-zero-gain
% changes I derived myself.
%
% Please note that a pole and a zero at the same place exactly cancel.
% This is significant for High Pass, Band Pass and Band Stop filters
% which create numerous extra poles and zeros, most of which cancel.
% Those which do not cancel have a 'fill-in' effect, extending the
% shorter of the sets to have the same number of as the longer of the
% sets of poles and zeros (or at least split the difference in the case
% of the band pass filter).  There may be other opportunistic
% cancellations but I will not check for them.
%
% Also note that any pole on the unit circle or beyond will result in
% an unstable filter.  Because of cancellation, this will only happen
% if the number of poles is smaller than the number of zeros and the
% filter is high pass or band pass.  The analytic design methods all
% yield more poles than zeros, so this will not be a problem.
%
% References:
%
% Proakis & Manolakis (1992). Digital Signal Processing. New York:
% Macmillan Publishing Company.

% Author: Paul Kienzle <pkienzle@users.sf.net>

% 2000-03-01 pkienzle@kienzle.powernet.co.uk
%       leave transformed Sg as a complex value since cheby2 blows up
%       otherwise (but only for odd-order low-pass filters).  bilinear
%       will return Zg as real, so there is no visible change to the
%       user of the IIR filter design functions.
% 2001-03-09 pkienzle@kienzle.powernet.co.uk
%       return real Sg; don't know what to do for imaginary filters
function [Sz, Sp, Sg] = sftrans(Sz, Sp, Sg, W, stop)

if (nargin ~= 5)
    usage('[Sz, Sp, Sg] = sftrans(Sz, Sp, Sg, W, stop)');
end;

C = 1;
p = length(Sp);
z = length(Sz);
if z > p || p == 0
    error('sftrans: must have at least as many poles as zeros in s-plane');
end

if length(W)==2
    Fl = W(1);
    Fh = W(2);
    if stop
        % ----------------  -------------------------  ------------------------
        % Band Stop         zero: b ? sqrt(b^2-FhFl)   pole: b ? sqrt(b^2-FhFl)
        %        S(Fh-Fl)   pole: ?sqrt(-FhFl)         zero: ?sqrt(-FhFl)
        % S -> C --------   gain: -x                   gain: -1/x
        %        S^2+FhFl   b=C/x (Fh-Fl)/2            b=C/x (Fh-Fl)/2
        % ----------------  -------------------------  ------------------------
        if (isempty(Sz))
            Sg = Sg * real (1./ prod(-Sp));
        elseif (isempty(Sp))
            Sg = Sg * real(prod(-Sz));
        else
            Sg = Sg * real(prod(-Sz)/prod(-Sp));
        end
        b = (C*(Fh-Fl)/2)./Sp;
        Sp = [b+sqrt(b.^2-Fh*Fl), b-sqrt(b.^2-Fh*Fl)];
        extend = [sqrt(-Fh*Fl), -sqrt(-Fh*Fl)];
        if isempty(Sz)
            Sz = [extend(1+rem([1:2*p],2))];
        else
            b = (C*(Fh-Fl)/2)./Sz;
            Sz = [b+sqrt(b.^2-Fh*Fl), b-sqrt(b.^2-Fh*Fl)];
            if (p > z)
                Sz = [Sz, extend(1+rem([1:2*(p-z)],2))];
            end
        end
    else
        
        % ----------------  -------------------------  ------------------------
        % Band Pass         zero: b ? sqrt(b^2-FhFl)   pole: b ? sqrt(b^2-FhFl)
        %        S^2+FhFl   pole: 0                    zero: 0
        % S -> C --------   gain: C/(Fh-Fl)            gain: (Fh-Fl)/C
        %        S(Fh-Fl)   b=x/C (Fh-Fl)/2            b=x/C (Fh-Fl)/2
        % ----------------  -------------------------  ------------------------
        Sg = Sg * (C/(Fh-Fl))^(z-p);
        b = Sp*((Fh-Fl)/(2*C));
        Sp = [b+sqrt(b.^2-Fh*Fl), b-sqrt(b.^2-Fh*Fl)];
        if isempty(Sz)
            Sz = zeros(1,p);
        else
            b = Sz*((Fh-Fl)/(2*C));
            Sz = [b+sqrt(b.^2-Fh*Fl), b-sqrt(b.^2-Fh*Fl)];
            if (p>z)
                Sz = [Sz, zeros(1, (p-z))];
            end
        end
    end
else
    Fc = W;
    if stop
        % ----------------  -------------------------  ------------------------
        % High Pass         zero: Fc C/x               pole: Fc C/x
        % S -> C Fc/S       pole: 0                    zero: 0
        %                   gain: -x                   gain: -1/x
        % ----------------  -------------------------  ------------------------
        if (isempty(Sz))
            Sg = Sg * real (1./ prod(-Sp));
        elseif (isempty(Sp))
            Sg = Sg * real(prod(-Sz));
        else
            Sg = Sg * real(prod(-Sz)/prod(-Sp));
        end
        Sp = C * Fc ./ Sp;
        if isempty(Sz)
            Sz = zeros(1,p);
        else
            Sz = [C * Fc ./ Sz];
            if (p > z)
                Sz = [Sz, zeros(1,p-z)];
            end
        end
    else
        % ----------------  -------------------------  ------------------------
        % Low Pass          zero: Fc x/C               pole: Fc x/C
        % S -> C S/Fc       gain: C/Fc                 gain: Fc/C
        % ----------------  -------------------------  ------------------------
        Sg = Sg * (C/Fc)^(z-p);
        Sp = Fc * Sp / C;
        Sz = Fc * Sz / C;
    end
end

return

%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Copyright (C) 1994, 1995, 1996, 1997, 1998, 2000, 2002, 2004, 2005,
%               2006, 2007, 2008, 2009 John W. Eaton
%
% This file is part of Octave.
%
% Octave is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or (at
% your option) any later version.
%
% Octave is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Octave; see the file COPYING.  If not, see
% <http://www.gnu.org/licenses/>.

% -*- texinfo -*-
% @deftypefn {Function File} {} postpad (@var{x}, @var{l}, @var{c})
% @deftypefnx {Function File} {} postpad (@var{x}, @var{l}, @var{c}, @var{dim})
% @seealso{prepad, resize}
% @end deftypefn

% Author: Tony Richardson <arichard@stark.cc.oh.us>
% Created: June 1994

function y = postpad (x, l, c, dim)

if nargin < 2 || nargin > 4
    %print_usage ();
    error('wrong number of input arguments, should be between 2 and 4');
end

if nargin < 3 || isempty(c)
    c = 0;
else
    if ~isscalar(c)
        error ('postpad: third argument must be empty or a scalar');
    end
end

nd = ndims(x);
sz = size(x);
if nargin < 4
    % Find the first non-singleton dimension
    dim  = 1;
    while dim < nd+1 && sz(dim)==1
        dim = dim + 1;
    end
    if dim > nd
        dim = 1;
    elseif ~(isscalar(dim) && dim == round(dim)) && dim > 0 && dim< nd+1
        error('postpad: dim must be an integer and valid dimension');
    end
end

if ~isscalar(l) || l<0
    error ('second argument must be a positive scalar');
end

if dim > nd
    sz(nd+1:dim) = 1;
end

d = sz(dim);

if d >= l
    idx = cell(1,nd);
    for i = 1:nd
        idx{i} = 1:sz(i);
    end
    idx{dim} = 1:l;
    y = x(idx{:});
else
    sz(dim) = l-d;
    y = cat(dim, x, c * ones(sz));
end

return
