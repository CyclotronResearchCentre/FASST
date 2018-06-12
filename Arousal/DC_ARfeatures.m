function feature = DC_ARfeatures(varargin)

if nargin<1
    'error --- no data or ar parameters'
elseif nargin == 1
    ar = varargin{1};
    fs = 1;
elseif nargin == 2
    ar = varargin{1};
    fs = varargin{2};
end 


% From Franaszczuk and Blinowska 1985
% Searching for parameters C and alpha to easily
% interpret the transfert function of the system. The
% system found as describing the neuronal system: a
% superposition of the natural frequency modes of diff
% amplitude and damping coefficients.
% Fonction de transfert du système linéaire AR: A^-1 = H
% A = feature.ar*z.^[0:-1:-order]';
% les poles de la fonction de transfert sont les
% solutions du dénominateur.
sol = 1./roots(fliplr(ar));
feature.oscillatory = real(sol)>0 & imag(sol)~=0;
feature.sol = sol;
%%----------------------------------------------
% Demonstration avec explication - version plus rapide
% dessous
% ----------------------------------------------
% On va écrire la fonction de transfert avec son
% dénominateur sous forme factorielle grâce aux
% solutions trouvés.
% den = 1;
% for j = 1:order
%     den = den*(z-sol(j));
% end
% H = z^order/den;
% % Supposant que les pôles soient tous simples, on veut
% % réécrire H sous forme d'une somme => nécessité de
% % trouver la Constante C d'amortissement propre à
% % chaque pôle. 
% for j = 1:order
%     C(j) = double(limit((z-sol(j))*H/z,z,sol(j)));
% end
% alpha = log(sol);
% % system stable when Beta > 0 => real(A)<0
% Beta = -real(alpha);
% Phi = double(abs(atan(imag(C)./real(C))));
% Omega = abs(imag(alpha))/(2*pi)*fs;
% B = double(2*abs(C));
% Omega_o2 = Beta.^2 + Omega.^2;
% Omega_o = (Omega_o2).^0.5;
% S = varnoise./(abs(ar*exp(-1i*[0:order])'*exp(w)).^2);
% plan s: usefull for the detection and follow-up of
% the EEG changes under the influence of different
% experimental condition
% plot(Beta,Omega_o,'*')
% The system can be represented as being a set of parallel filters % corresponding to the generators of the different rhythms. 

feature.magnitude = abs(sol);
[Maxmag imaxmag] = max(abs(sol));
ws = abs(imag(log(sol./abs(sol))));

%  w = 2*pi*f;
%  w = ws*fs
%  f = ws*fs/2pi
feature.frequency  = fs*ws./(2*pi);
damping = -log(feature.magnitude);
feature.naturalfrequency = (damping.^2 + feature.frequency.^2).^0.5;
feature.damping = fs*damping;    
