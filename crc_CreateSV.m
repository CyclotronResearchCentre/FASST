function SV = crc_CreateSV(data,t_wind,avg,el_set,type,maxphi,ndisp,el_label)

%_________________________________________________________________________
%
% FORMAT SV = crc_CreateSV(data,t_wind,avg,el_set,type,maxphi)
%
% Function that creates the "Spline View" (SV) structure
% The values at the electrodes are interpolated over the surface of a sphere.
% Interpolation can also be obtained for the surface Laplacian.
% IN:
%   - data   : data matrix (Nel x Ntb)
%   - t_wind : time window for data are 
%   - avg    : data can be averaged over t_wind and only one image produced (avg=1)
%              or one image is created per time bin (avg=0), SV(1...Ntb).
%   - el_set : electrodes set, an index refering to the sets defined 
%               in spm_eegfp_electrset, or the electrodes coordinates (on a shere),
%               or the name of the electrodes as in the data.
%   - types  : 'simple' EEG interpolation (1) or Laplacian (2)
%   - maxphi : Maximum phi angle used for the interpolation and display,
%               by defaut 7*pi/12 = 105°.
%   - ndisp  : defines the resoltion of the image (2*ndisp+1, 2*ndisp+1)
%   - el_label: labels of the electrodes when posistion is given in el_set
%
% Sub-functions
%--------------
% function [ x, y, z ] = FlatProjection( maxphi, n )
%   Generate a grid of location on a *sphere* of radius 1.
%   It is put into (2n+1)x(2n+1) matrices
% function SS = CreateSS(f,el_xyz,order,Niter)
%   Prepare the spline interpolation
% function y = Interpol(a, x, opt)
%   Proceed to the interpolation of the values at the electrodes,
%   Calculation is based on a spherical spline interpolation, 
%   i.e. electrodes are supposed to be on a sphere of radius 1
%
% All this has been adapted from routines written by Morgan Willis
% when he was working with Mick Rugg at the ICN, London, UK.
%_________________________________________________________________________
% Adapted by c.phillips@ulg.ac.be, 2002.12.11
% Last modified,
% - 2004/09/06, use of electr_sets function
% - 2004/11/09, adapted for spm_eeg
% - 2005/03/16, input the name of the channels, no worry about the 
%               order defined in spm_eeg_electr.
% - 2008/01/19, adapted for SPM8 (jschrouff@doct.ulg.ac.be)
%__________________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by J.Schrouff & C. Phillips, 2009.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$


gui=0;
if nargin<8, el_label=[]; end
if nargin<7, ndisp = 50;   end
if nargin<6, maxphi = 10*pi/12;   end
if nargin<5, type = 1;          end
if nargin<4, el_set = {};       end
if nargin<3, avg = 0;           end
if nargin<2, t_wind = [];       end
if nargin<1
    data = [];
    gui=1;
end

% Interpolation parameters
order = 3;
Niter = 20;

if gui
    Pdata = spm_select(1, 'mat', 'Select cleaned EEG file','' ,pwd,'i_*.*');
    load(Pdata);
end
[DNchan,Dtb] = size(data);

t_series = [];
if gui
    % Use GUI
	flag=0; pos = 1;
	while flag==0
        [wind_be,pos] = spm_input(['Time window to use'],pos,'e',[1 Dtb]);
        if length(wind_be)==1
            if wind_be(1)>0 & wind_be(1)<=Dtb
                wind_be(2) = wind_be(1); 
                flag = 1;
            end
        elseif length(wind_be)==2
            if wind_be(2)>=wind_be(1) & all(wind_be>0) & all(wind_be<=Dtb)
                flag = 1;
            end
        elseif length(wind_be)>2
            if wind_be(2)>=wind_be(1) & all(wind_be>0) & all(wind_be<=Dtb)
                flag = 1;
                t_series = wind_be;
            end
        end
	end
    t_wind = [wind_be(1) wind_be(end)];
elseif isempty(t_wind)
    % Use complete dataset
    t_wind = [1 Dtb];
elseif length(t_wind)==1
    % Use a single time point
    t_wind = [t_wind t_wind];
elseif length(t_wind)>2
    % Use a subset of time points
    t_series = t_wind;
    t_wind = [t_wind(1) t_wind(end)];
end
if isempty(t_series)
    t_series = t_wind(1):t_wind(2);
end
NdispD = length(t_series);

if isempty(el_set)
    % No electrode info provided, use spm_eegfp_electr
    [set_Nel,set_name] = spm_eegfp_electrset;
    el_set = find(set_Nel==DNchan);
    if isempty(el_set)
        error('Sorry can''t find right electrodes set');
    elseif length(el_set)==1
        [el_xyz,el_name] = spm_eegfp_electrset(el_set) ;
    else
        if gui
            [sel_set,pos] = spm_input('Electrode sets','+1','m',set_name(el_set,:))
            el_set = el_set(sel_set);
        else
            warning('2 sets at least have the same number of electrodes.')
        end
        [el_xyz,el_name] = spm_eegfp_electrset(el_set(1)) ;
    end
elseif length(el_set)==1
    % Select electrode set in spm_eegfp_electr
    [el_xyz,el_name] = spm_eegfp_electrset(el_set) ;
elseif size(el_set,1)==3 && size(el_label,1)==DNchan
    % Coordinates are directly provided

    el_xyz = el_set;
    el_name=el_label;
elseif ischar(el_set)
    % Case where user provides electrode names in order of the data
    % => make sure order coresponds to setup in spm_eegfp_electr
    ENchan = size(el_set,2);
    if ENchan~=DNchan,
        error('Nr of electrode names does correspond to nr of channels in data!')
    end
    [set_Nel,set_name] = spm_eegfp_electrset;
    el_set_i = find(set_Nel==DNchan);
    if isempty(el_set_i)
        error('Sorry can''t find right electrodes set');
    elseif length(el_set_i)==1
        [el_xyz,el_name] = spm_eegfp_electrset(el_set_i) ;
    end
%     Re-order channels.
    ord = zeros(ENchan,1);
    for ii=1:ENchan
        for j=1:size(el_name,2)
             if ~iscell(el_name{j})
                 if strcmpi(deblank(el_label(ii,:)),(el_name{j}));
                     ord(ii)=j;
                 end
             else
                for k=1:size(el_name{j},2)
                    if strcmpi(deblank(el_label(ii,:)),(el_name{j}{k}));
                        ord(ii)=j;
                    end
                end
            end
        end
    end
    el_xyz = el_xyz(:,ord);
    el_name = el_set;
else
    error('Wrong electrode set');
end

if iscell(el_name)
    el_name = char(el_name);
end

if gui
    [avg,pos]  = spm_input('Average over time window :','+1','yes|no',[1 0],1);
    [type,pos] = spm_input('Type of interpolation :','+1','EEG|CurrentDensity',[1 2],1);
end

SV = struct('maxphi',maxphi,'ndisp',ndisp,'latency',t_wind, ...
            't_series',t_series,'range',[0 0],'auto',1,'type',type,'normalised',[], ...
            'avg',avg,'hFigure',[],'hAxis',[],'hColourBar',[],'el_name',el_name);
[x, y, z ] = FlatProjection( SV.maxphi, SV.ndisp ); %grid of locations on a sphere of radius 1.
SV.normalised = [ x(:), y(:), z(:) ]';
s = [2*SV.ndisp+1, 2*SV.ndisp+1];
if avg
    f = mean(data(:,t_series),2);	
    SV.spline = CreateSS(f,el_xyz,order,Niter);
    interpolated = zeros( s );
    interpolated(:) = NaN;
    SV.interpolated = interpolated;
    SV.interpolated = reshape(Interpol(SV.spline,SV.normalised,SV.type),s(1),s(2));
    if SV.auto
        vM = max(abs(SV.interpolated(:)));
        SV.range = [-vM vM];
    else
        SV.range = [-1 1];
    end
else
    SV.spline = CreateSS(data(:,t_series),el_xyz,order,Niter);
    
    interpolated = zeros( [s,NdispD] );
    interpolated(:) = NaN;
    SV.interpolated = interpolated;
    SV.interpolated = reshape(Interpol(SV.spline,SV.normalised,SV.type),s(1),s(2),NdispD);
    if SV.auto
        vM = max(max(max(abs(SV.interpolated))));
        SV.range = [-vM vM];
    else
        SV.range = [-1 1];
    end
end

if gui % Draw the result if I use the gui.
    crc_DrawSV(SV);
end

return
%________________________________________________________________________
%________________________________________________________________________
%
% SUBFUNCTIONS
%________________________________________________________________________
%________________________________________________________________________
%________________________________________________________________

%________________________________________________________________
%
% Coordinates on a sphere
%------------------------

function [ x, y, z ] = FlatProjection( maxphi1, n )
% Generate a grid of location on a *sphere* of radius 1.
% It is put into (2n+1)x(2n+1) matrices

s = [ 2*n+1, 2*n+1 ];
x = zeros( s );
x(:,:) = NaN;
y = zeros( s );
y(:,:) = NaN;
z = zeros( s );
z(:,:) = NaN;

for iy = -n:n
    nx = round( sqrt( n*n - iy*iy ) );
    for ix = -nx:nx

        ir = sqrt( ix*ix + iy*iy );
        phi = maxphi1 * ir / n;
        costheta = 0;
        sintheta = 0;
        if ir > 0
            costheta = iy/ir;
            sintheta = -ix/ir;
        end
        
        x(ix+n+1,iy+n+1) = costheta*sin(phi);
        y(ix+n+1,iy+n+1) = sintheta*sin(phi);
        z(ix+n+1,iy+n+1) = cos(phi);
 
    end
end
return
%________________________________________________________________
%
% Prepare the spline interpolation
%---------------------------------

function SS = CreateSS(f,el_xyz,order,Niter)

% if ( size(f,2) ~= 1 )
%     error( 'Values to interpolate must be a column vector' );
% elseif ( size(el_xyz,1) ~= 3 )
%     error( 'Electrodes coordinates as a 3xNel matrix' );
% elseif ( size(f,1) ~= size(el_xyz,2) )
%     error( 'There must be one value to interpolate per electrode' );
% end

SS.X = el_xyz;
SS.m = order;
SS.n = Niter;
SS.nPoints = size(f,1);
SS.c0 = 0;
SS.c = [];

% Calculate c0 and c
nSize = SS.nPoints + 1;
A = zeros( nSize, nSize );
A( 1, 2:nSize ) = 1;
A( 2:nSize, 1 ) = 1;
A( 2:nSize, 2:nSize ) = g( SS.X' * SS.X, SS.m, SS.n );
Nf = size(f,2);
if Nf==1
    B = zeros( nSize, 1 );
    B( 2:nSize ) = f;
    x = A\B;
    SS.c0 = x(1);
    SS.c  = x(2:nSize);
else
    B = zeros( nSize, Nf );
    B( 2:nSize,: ) = f;
    x = A\B;
    SS.c0 = x(1,:);
    SS.c  = x(2:nSize,:)';
end

%________________________________________________________________
%
% Interpolation function
%-----------------------
function y = Interpol(a, x, opt)

% Proceed to the interpolation of the values at the electrodes,
% Calculation is based on a spherical spline interpolation, 
% i.e. electrodes are supposed to be on a sphere of radius 1
% The values at the electrodes can be interpolated
% or the spatial Laplacian.
% IN:
%   - a     = Spherical spline structure
%   - x     = coordinates of points where function is interpolated
%   - opt   = chose between normal (1) or Laplacian (2) interpolation

if nargin<3, opt=1; end

Ndat = length(a.c0);
y = zeros(size(x,2),Ndat);
if opt==1        % Simple interpolation
    GG = g( x' * a.X, a.m, a.n );
    if Ndat==1
        y = a.c0 + GG * a.c;
    else
        for ii=1:Ndat
            y(:,ii) = a.c0(ii) + GG * a.c(ii,:)';
        end
    end
elseif opt==2    % Interpolation of the Laplacian
    GG = g( x' * a.X, a.m-1, a.n );
    if Ndat==1
        y = GG * a.c;
    else
        for ii=1:Ndat
            y(:,ii) =  GG * a.c(ii,:)';
        end
    end
else
    error('Wrong interpolation type');
end
return

%________________________________________________________________
%
% Core of the interpolation.
%---------------------------

function G = g( X, m, n )

G = zeros( size(X) );
P = Legendre( X, n );
for i = 1:n
   a =  (2*i + 1) / ( i^m * (i+1)^m );
   G = G +  P{i+1}*a;
end
G = G / (4*pi);
    
return;

%________________________________________________________________
%
% Estimate the Legendre Polynomes
%--------------------------------
function P = Legendre( X, n )

P = cell(1,n+1);
P{1} = ones( size(X) );
P{2} = X;

for i = 2:n
    P{i+1} = ( (2*i-1)* X .* P{i} - (i-1)* P{i-1})/ i;
end


