function D = crc_par_cICAmx(Di)

% FORMAT function D = crc_par_cICAmx(Di)
%
% This function estimates the correction matrices, no corrected data are
% saved. To clean the data, up to the user to pick the "best" correction 
% matrix, either manually or using some automatic MI criteria.
%__________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by Y. Leclercq & C. Phillips
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

% NOTE:
% There may be a problem with the latest block of data if it's too short:
% because the loop goes block by block it could happen that the last block
% is only a few seconds long!
% To be addressed later.

crcdef = crc_get_defaults('par');

% Recover Peaks and badchan
badindex = Di.cache.badch;

mxes = cell(0);

% Size of the sample analyzed in one step in seconds (advice: 60 sec is the
% minimum.
sz = crcdef.bcgrem.size;
% Step in seconds (advice: step = 2/3 size)
step = crcdef.bcgrem.step;

D = Di;
fs = Di.fsample;

% Initialize struff
ideb=1;
ifin=1+round(sz)*fs;
Nit = 0;
totaltime = 0;
% Loop
while ifin<=Di.nsamples
    [W, time ] = crc_bcgrem(Di,ideb,ifin,badindex);
    mxes = [mxes W];
    totaltime = totaltime + time;
    ideb = ideb+round(step)*fs;
    ifin = ifin+round(step)*fs;
    if ifin<Di.nsamples & (Di.nsamples-ifin)<round(sz)*fs
        % ensuring last block is long enough by sticking it to previous one
        % if necessary.
        ifin = Di.nsamples;
    end
    Nit = Nit+1;
end
if ~Nit, error('Data too short to be corrected with cICA'); end

if ~isfield(D,'CRC')
    D.CRC = [];
end

if ~isfield(D.CRC,'cICA')
    cICA=cell(4,1);
    cICA{1,1}=step;
    cICA{2,1}=sz;
    cICA{3,1}=mxes;
    cICA{4,1}=badindex;
    D.CRC.cICA = cICA;
else
    Z=cell(4,1);
    Z{1,1}=step;
    Z{2,1}=sz;
    Z{3,1}=mxes;
    Z{4,1}=badindex;
    D.CRC.cICA=[D.CRC.cICA Z];
end

crc_time_converts(totaltime)

end

%% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BCGrem
function [Wcorr, time, z, mu, lambda] = ...
                   crc_bcgrem(Dmeeg,ideb,ifin,badindex,z1,mu1,lambda1,Wold)

% FORMAT function [Wcorr, time, z, mu, lambda] = ...
%              crc_bcgremoval(Dmeeg,ideb,ifin,badindex,z1,mu1,lambda1,Wold)
% Remove the BCG artefact from a stretch of data
%
% INPUT:
% - Dmeeg       : data structure
% - ideb        : index of beginning time bin \_ of time window to use
% - ifin        : index of ending time bin    /
% - badindex    : index of channels to "skip"
% - [z1, mu1, lambda1, Wold] : previous estimates of those variables
%__________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

crcdef = crc_get_defaults('par');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Channels used to build the constrain vector.
% They're taken from around the head
Refnames = crcdef.Refnames ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Loading of the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clock1 = clock;
fs     = Dmeeg.fsample;
data_o = Dmeeg(:,max(1,ideb):min(ifin,Dmeeg.nsamples));

%%%%%%%%%Filter the very low frequencies

zx = 1.5*fs;
zx = ones(1,zx)/zx;
data2 = conv2(data_o,zx,'same');
data = data_o-data2;
clear data2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Removes the badchannels
%

% %Find EOG channels
% EOGchan = eogchannels(Dmeeg);
selected_chan = setdiff(1:Dmeeg.nchannels, badindex);

% Gives the peaks of ECG and BCG template for each channel
allPeaks = Dmeeg.CRC.EKGPeaks;

Peaks = allPeaks(and(ideb<allPeaks,allPeaks<ifin)==1)-ideb+1;

data_osel = data_o(selected_chan,:);

%Check the peaks
lim = mean(diff(Peaks))/3;
list = find(diff(Peaks)>lim);
Peaks = Peaks([list length(Peaks)]);

fitted_art = mean_art(data_osel,Peaks);

data_sel = data(selected_chan,:);
selnames = chanlabels(Dmeeg,selected_chan);
refind = [];
for kk = 1:length(Refnames)
    search = find(strcmp(Refnames(kk),upper(selnames)));
    if ~isempty(search)
        refind = [refind search];
    end
end

% The multi-cICA algorithm is used to determine the Independent Component
% closer to BCG template for each channel
try
    %Detects NAN in fitted_art and replaces them with 0

    lNAN=find(isnan(fitted_art));
    fitted_art(lNAN)=zeros(size(fitted_art(lNAN),1),size(fitted_art(lNAN),2));


    % find movements artefact and reject
    %===================================
    rejdat = data_osel;
    sizerejwin = 0.5;
    tmp = find(abs(detrend(rejdat','constant'))>max(std(rejdat,0,2))*crcdef.bcgrem.scSNR);

    countbad = zeros(1,size(rejdat,1));
    Ntmp = size(rejdat,2);
    for ii=1:size(rejdat,1)
        countbad(ii) = sum( (tmp<ii*Ntmp+1) & (tmp >(ii-1)*Ntmp) );
    end

    tmp2 = mod(tmp,size(rejdat,2));
    tmp2(tmp2==0) = size(rejdat,2);
    tmp3 = unique(tmp2);

    kil=1;
    if ~isempty(tmp3)
        tmp4 = max(1,tmp3(kil)-sizerejwin*fs):tmp3(kil);
    else
        tmp4 =[];
    end

    while kil<length(tmp3)
        inbubble=find(and(tmp3<tmp3(kil)+sizerejwin*fs,tmp3>tmp3(kil)));
        if ~isempty(inbubble)
            tmp4 = [tmp4 (tmp4(end))+1:tmp3(inbubble(end))];
            kil = inbubble(end);
        else
            tmp4 = [tmp4  max(1,tmp3(kil)-sizerejwin*fs):min(tmp3(kil)+sizerejwin*fs,size(data,2))];
            kil=kil+1;
        end
    end
    toreject = unique(tmp4);
    tokeep=setdiff(1:size(data_sel,2),toreject);

    if numel(toreject)/numel(tokeep)>.5
        disp('Rejecting more than 50% of time bins from current time window!');
        error('Probably too much noise in some channel(s).');
    end
    
    icadata=data_sel(:,tokeep);
    icaref=fitted_art(refind,tokeep);

    if nargin>4
        [z mu lambda]=multicica(icadata,icaref,z1,mu1,lambda1);
    else
        [z mu lambda]=multicica(icadata,icaref);
    end

    disp('cICA done')

    % NB:   z are the weight vectors corresponding to BCG independent source
    %
    %       mu and lambda don't serve any purpose in the current execution.
    %       But they can be used to initialize the multicica algorithm for an
    %       input signal close to the one which has been analysed. It results
    %       in a faster convergence towards the solution.
    %
    %       For example, if we have a signal of total length 2 minutes 30
    %       seconds and analyse the signal from time 0 to 1 min 30 seconds, it
    %       will results in a particular z, mu and lambda. If we initialize the
    %       cICA algorithm with these to analyse signal from 1 min to 2 min 30
    %       seconds, we will have a faster converge than without specifying any
    %       initialization.
    %

    % Clustering using multiple execution of k-means clustering

    if size(data,1)>16
        Nclus=5;
    else
        Nclus=3;
    end

    [vmax,cmax,codemax]=multikmeans(z,Nclus,crcdef.NitKmeans);

    disp('Clustering done.')

    if vmax<0.55
        fprintf('Warning : no BCG artifact clearly identified')
    end

    nbreart=size(cmax,1);

    % Construction of an orthonormal base in the whithened data space
    w=cmax';
    W=icaorth(data_sel,w);
    %W=icaorth(data,w);
    disp('Gramm-Schmidt base reconstruction done.')
    disp(' ')
    % Computation of the independent component

    Winv=pinv(W);

    % Computation of the cleaned signal
    Winvcorr=Winv;
    Winvcorr(:,1:nbreart)=zeros(size(Winv,1),nbreart);
    Wcorr=Winvcorr*W;
catch
    if exist('Wold','var')
        Wcorr=Wold;
        z=z1;
        mu=mu1;
        lambda=lambda1;
    else
        disp(' ')
        disp('Warning: Algorithm fails to correct the EEG recordings, probably due to')
        disp('too much artefacts')
        Wcorr=ones(numel(selected_chan));
        z=[];
        mu=[];
        lambda=[];
    end
end

% tmp = Wcorr*data_sel;
% newdata = zeros(size(data_o));
% newdata(selected_chan,:) = tmp;
% other_chan = setdiff(1:size(data_o,1),selected_chan);
% newdata(other_chan,:) = data_o(other_chan,:);

time=clock-clock1;
time=time(length(time)-2)*3600+time(length(time)-1)*60+time(length(time));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MultiICA
function [w mu lambda]=multicica(data,ref,init,mu,lambda)
% FORMAT [wout muout lambdaout] = multicica(data,ref,init,mu,lambda)
%
% Calculate the "constrained ICA" (cICA) decomposition of the dat, given
% some reference signal.
%
% Input argument :
% data   - the mixed signals matrix (N x time)
% ref    - the M references matrix (M x time)
%
% Optional input argument :
% init   - the initialization of the weight matrix w
% mu     - \ are initialization of parameters
% lambda - /
% (cf Paper of Lu & Rajapakse "Approach and applications of constrained
% ICA" for more details)
%
% Output argument :
% wout      - the value of the weighting matrix after convergenece
% muout     - \ are the value of mu and lambda after convergence
% lambdaout - /
%
% Reference:
% Lu & Rajapakse, Approach and applications of constrained ICA
% IEEE Transactions on neural networks, 16(1), pp 203-212, 2005
%__________________________________________________________________
% close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Whitening of the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[E, D]=pcamat(data);
[newvect,wmx,dwmx]=whitenv(data,E,D);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

siz = size(newvect);
vectorsize = siz(1); % Number of channels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Core of cICA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Building reference
r2=ref;
nbreref=size(r2,1); %Number of references

% Correlation between w(i-1) et w
corrw=[];

% Number of Iteration
it=0;
normgold=0;

%%%Initialization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% What we want
if nargin<3
    w=randn(nbreref,siz(1));                     % weight vectors
else
    w=dwmx'*init';
    w=w';
end
if nargin<4
    mu=zeros(nbreref,nbreref)';                   % 2nd Lagrangian multiplier
end
if nargin<5
    lambda=zeros(nbreref,nbreref)';               % 1st Lagrangian multiplier
end
% Constant

rho=eye(nbreref);                  % positive number => function ???
gamma=eye(nbreref);                % learning rate of the lagrangian multiplier.

% dynamic 'constant'
ksi=0.1;                    % Threshold of the difference between the reference and the output
neta=0.5;                   % Positive learning rate to avoid uncertainty of convergence

notok=true;
diffrgold=-inf;

while notok
    y=w*newvect;           %Calculation of the output based on the new w
    Rxx=newvect*newvect';   %Autocorrelation mx

    g=[];
    dg=[];
    ddg=[];

    i=1;
    while i<nbreref+1

        g=[g ; norm(y(i,:)-r2(i,:)) - ksi];
        dg=[dg ; ((y(i,:)-r2(i,:))/(g(i)+ksi))];
        ddg=[ddg ; mean(((g(i)+ksi)^2-((y(i,:)-r2(i,:))))/(g(i)+ksi)^3)];
        i=i+1;
    end

    %g=norm(nor-r2)-ksi;       % Constrain function
    %dg=(nor-r2)/(g+ksi);      % 1st derivative of the constrain function
    %ddg=mean(((g+ksi)^2-(nor-r2))/(g+ksi)^3); %Expectation of the 2nd derivative of the cfun

    h=mean((y.^2)')'-1;     %Bound the negentropy and the weight vector (y=w*x)
    dG=tanh(y);             %1der of Nonlinear function used in the approx of negentropy
    ddG=1 - dG.^2;          %2der of Nonlinear function used in the approx of negentropy
    rhob=rho;               %*sign(sum(log(cosh(y))-log(cosh(randn(size(y))))));

    Lw=rhob*dG*newvect'-0.5*mu*dg*newvect'-lambda*y*newvect';

    s=rhob*mean(ddG')'-0.5*mu*ddg-lambda*(ones(1,nbreref))';
    s=1./s;
    s=s*ones(1,length(s)).*eye(length(s));

    %update
    wold=w;
    muold=mu;
    lambdaold=lambda;

    w=w-neta*(s*Lw*Rxx^(-1));
    mu=max(0,mu+gamma.*(g*ones(1,length(g)).*eye(length(g))));
    lambda=lambda+gamma*(h*ones(1,length(g)).*eye(length(g)));

    %find NaN in the w matrix and replace with random value
    z=find(not(w>0|w<0|w==0));
    w(z)=randn(size(w(z),1),size(w(z),2));
    corrwi=norm(wold*w')/(norm(wold)*norm(w'));
    diffrg=abs((normgold-norm(g))/normgold);
    normgold=norm(g);

    % Convergence test
    if or((1-corrwi)<10e-6&(abs(diffrg)<25e-6), diffrg<0.1 & (1-corrwi)<1e-6 & it>3 & diffrg>diffrgold)
        notok=false;
    end
    diffrgold=diffrg;
    it=it+1;
    if it>5 & not(diffrg>0|diffrg<0|diffrg==0)
        w=[];
        break
    end
end

%normalization of weight vectors
for f=1:nbreref
    w(f,:)=w(f,:)/norm(w(f,:));
end

%Project w in the dewithen space.
if length(w)>1
    w=w*wmx;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ICAorthogonalization
function W=icaorth(datanotw,guess)

% Orthogonalization of the ICA components.
%__________________________________________________________________
%The guess vector should be column-vectors

[E, D] = pcamat(datanotw);
[data, wm, dwm] = whitenv(datanotw, E, D);


TrueNumberofIC=size(datanotw,1);

NumberofIC=size(data,1);
numsamples=size(data,2);

W=zeros(NumberofIC,NumberofIC);

if nargin<2
    guess=[];
else
    W(1:NumberofIC,1:size(guess,2))=dwm'*guess;
    for i = 1:size(guess,2)
        W(:,i)=W(:,i)/norm(W(:,i));
    end

    initguess=W(1:NumberofIC,1:size(guess,2));
    for i = 2:size(guess,2)
        W(:,i)=W(:,i)-W(:,1:i-1)*W(:,1:i-1)'*W(:,i);
        W(:,i)=W(:,i)/norm(W(:,i));
    end
end


i=size(guess,2)+1;
failures=0;
%greatfailure=0;

while i<NumberofIC+1
    diffw=inf;
    w=ones(NumberofIC,1);
    w=w-W*W'*w;
    randn('state',sum(100*clock));
    w=w.*randn(NumberofIC,1);

    w=w/norm(w);
    winit=w;
    wold=ones(NumberofIC,1);
    %wold=wold/norm(wold);

    it=1;
    %For one ic
    while diffw>1e-12

        %if greatfailure==0
        %        hypTan = tanh(data' * w);
        %       w = (data * hypTan - sum(1 - hypTan .^ 2)' * w)/numsamples;
        %end

        w=w/norm(w,2);
        w=w-W*W'*w;
        w=w/norm(w,2);

%         wold;
        diffw=1-abs((w'*wold)/(norm(wold,2)*norm(w,2)));
        %pause
        %fprintf('Correlation between w and wold (%d).%d.%d\n', diffw,it,failures)
        wold=w;
        if it>1000
            if failures>3
                break
            end
            w=randn(NumberofIC,1);
            w=w/norm(w,2);
            wold=ones(NumberofIC,1);
            it=1;
            failures=failures+1;
        end
        it=it+1;
    end

    if failures>3
        break
        %greatfailure=1;
    end
    if diffw==0
        winit';
    end
    W(:,i)=w;
    i=i+1;
    failures=0;
end

initguess=initguess*diag(1./(sqrt(sum(initguess.^2))));
W(1:NumberofIC,1:size(guess,2))=initguess;

W = W' * wm;
W=W(1:i-1,:);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Multi-kmeans
function [vmax,cmax,codemax]=multikmeans(set,k,nbreex)

% FORMAT [vmax,cmax,codemax]=multikmeans(set,k,nbreex)
% Apply k-mean clustering multiple times to obtain the best clustering.
%
% INPUT
% set     - data set in matrix format, each column being a data vector
% k       - number of cluster
% nbreex  - number of execution of the clustering
%
% OUTPUT
% vmax    - variance of clusters, maximum of all repetitions
% cmax    - cluster centroids, optimum of all repetitions
% codemax - vector clustering, optimum of all repetitions
%__________________________________________________________________________

if nargin < 2
    k=2;
end

if nargin < 3
    nbreex=100;
end

critere=true;
diff=[];
vectmax=[];
cmax=[];
vmax=0;
codemax=[];

while critere
    lastcmax=cmax;
    lastvmax=vmax;
    lastcodemax=codemax;
    i=1;
    while i<nbreex+1
        [var,centroids,code]=kmeans(set',k);
        if var > vmax
            cmax=centroids;
            vmax=var;
            codemax=code;
        end
        i=i+1;
        var;
    end
    vectmax=[vectmax vmax];
    if length(vectmax)==1
    else
        diff=[diff (vectmax(length(vectmax))-vectmax(length(vectmax)-1))];
        if diff(length(diff))<0.06
            critere=false;
        end
    end
    k=k+1;
end

cmax=lastcmax';
vmax=lastvmax;
codemax = lastcodemax;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% K-mean clustering
function [variance,centroids,code]=kmeans(set,k)

% K-mean clustering
% FORMAT [variance,centroids,code]=kmeans(set,k)
%
% Input:
% set   - (m,n) set of vectors to cluster merged into a big matrix
%           v1(1) v2(1) ... vn(1)
%           v1(2) v2(2) ... vn(2)
%           ...   ...   ... ...
%           v1(m) v2(m) ... vn(m)
% k     - the number of cluster you want
%
% Output:
% variance  - variance of the clusters
% centroids - centroids of the k clusters
% code      - cluster for the n vectors
%
% NB: This algorithm is quite rapid and the results
% depends on its initialization. So it is recommended
% to run this function multiple times.
%__________________________________________________________________

[dim,number] = size(set);
%   dimension of the space we work in
%   number of vector in the set
vari=var(set,0,2);

% Initialization : random centroids
newcentr=[];
i=1;

while i<k+1
    rand('state',sum(100*clock));
    ind=round((number-1)*rand(1))+1;
    newcentr=[newcentr set(:,ind)];
    i=i+1;
end

centroids=0*newcentr;
arraylast=[0 0 0 0];
while or(centroids~=newcentr,min(arraylast)==0)
    centroids=newcentr;
    i=1;
    j=1;
    array=zeros(k,number);
    arraylast=zeros(k,1);
    while i < number+1
        vect=centroids(:,1)-set(:,i);
        diffmin=norm(vect,dim);
        kmin=1;
        j=2;
        while j < k+1
            vect=centroids(:,j)-set(:,i);
            diff=norm(vect,dim);
            if diff<diffmin
                diffmin=diff;
                kmin=j;
            end
            j=j+1;
        end
        array(kmin,arraylast(kmin)+1)=i;
        arraylast(kmin)=arraylast(kmin)+1;
        i=i+1;
        j=1;
    end

    %arraylast
    newcentr=[];
    i=1;
    j=1;

    while i<k+1
        temp=zeros(dim,1);
        while j<arraylast(i)+1
            temp=temp+set(:,array(i,j));
            j=j+1;
        end

        if temp==zeros(dim,1)
            temp=set(:,1+round(abs((number-1)*rand(1))));
            while centroids(:,i)==temp
                temp=set(:,1+round((number-1)*rand(1)));
            end
            newcentr=[newcentr temp];
        else
            newcentr=[newcentr temp/arraylast(i)];
        end
        i=i+1;
        j=1;
    end
end

% End of clusterisation - Estimation of the cluster's variance
i=1;
j=1;
% varmx=[];
variance=zeros(dim,1);
code=zeros(1,number);

while i<k+1
    x=[];
    while j<arraylast(i)+1
        x=[x set(:,array(i,j))];
        code(array(i,j))=i;
        j=j+1;
    end
    if arraylast(i)>2
        variance=variance+arraylast(i)*var(x')';
    end
    i=i+1;
    j=1;
end

variance=mean((ones(dim,1)-(variance/number)./vari));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the mean artefact of the pulse artefact given ECG peaks.
function outdat = mean_art(d,pks)

Nit = 20;

step = max(diff(pks));

outdat = zeros(size(d));

warhouse = zeros(size(d,1),(length(pks)-1)*step);

for ii = 1:length(pks) - 1

    warhouse(:,1+step*(ii-1):step*ii) = ...
        remap(d(:,pks(ii):pks(ii+1)-1),step);

end

Cube = reshape(warhouse(:,1:step*Nit),[size(warhouse,1) step Nit]);
iCube = 1;

for ii = 1:length(pks) - 1

    Cube(:,:,iCube) = warhouse(:,1+step*(ii-1):step*ii);
    warhouse(:,1+step*(ii-1):step*ii) = mean(Cube,3);

    avg = mean(warhouse(:,1+step*(ii-1):step*ii),2)*ones(1,step);

    warhouse(:,1+step*(ii-1):step*ii) = warhouse(:,1+step*(ii-1):step*ii) - avg;

    outdat(:,pks(ii):pks(ii+1)-1) = ...
        remap(warhouse(:,1+step*(ii-1):step*ii),pks(ii+1)-pks(ii));

    iCube = mod(iCube+1,20);
    if iCube == 0
        iCube = 20;
    end

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remap the data to the specified number of points.
function outdat=remap(data,step)

outdat = zeros(size(data,1),step);
sz = size(data,2)-1;
step=step-1;
for i=1:step

    w = mod((i-1)*sz/(step),1);
    idx = 1+floor((i-1)*sz/(step));

    outdat(:,i)= data(:,idx)*(1-w) + data(:,idx+1)*(w);

end
outdat(:,end)=data(:,end);

end