% fmrib_pas() - Remove pulse
%   artifacts from EEG data collected inside an MRI machine.
%
%   This program removes pulse artifacts from EEG data collected
%   collected inside the MRI scanner. A choice of methods is available.
%
%   The first choice is Optimal Basis Set.
%   This method aligns all the pulse artifacts, in each EEG channel
%   separately, in a matrix and performs a Principal Component Analysis
%   (PCA) on the data.  The first N PCs (the Optimal Basis Set) are then
%   fitted to each artifact in that channel. The process is repeated for
%   each channel.  The other three methods are based on [Allen98] with
%   improvements to better capture and subtract the artifacts.
%   Basically, in these methods a statistical measurement is used to find a
%   template for the artifact  at each heart beat.  A window of 30 artifacts
%   centred around the artifact being processed is used.  The 30 artifacts
%   are processed by taking the mean ('mean' method), the median
%   ('median' method) or a Gaussian-weighted
%   ('gmean' method to emphasise shape of current artifact) mean.  It is
%   recommended to use the first ('obs') method as it generally better
%   fits the true artifact.
%
% Usage:
%    >> EEGOUT= fmrib_pas(EEG,qrsevents,method)
% or >> EEGOUT= fmrib_pas(EEG,qrsevents,method,npc)
%
% Inputs:
%   EEG: EEGLAB data structure
%   qrsevents: vector of QRS event locations specified in samples.
%   method: 'obs' for artifact principal components.  You need to specify
%               'npc', which is the number of PC to use in fitting the
%               artifacts.  If unsure, use 4.
%           'mean' for simple mean averaging.
%           'gmean' for Gaussian weighted mean.
%           'median' for median filter.
%
%
%
%   [Allen98] Allen et.al., 1998, Identification of EEG events in the MR
%   scanner: the problem of pulse artifact and a method for its
%   subtraction. NeuroImage8,229-239.
%
%
%
% Also See pop_fmrib_pas
%
%   Author:  Rami K. Niazy
%
%   Copyright (c) 2004 University of Oxford

%123456789012345678901234567890123456789012345678901234567890123456789012
%
% Copyright (C) 2004 University of Oxford
% Author:   Rami K. Niazy, FMRIB Centre
%           rami@fmrib.ox.ac.uk
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
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% JUN 3, 2005
% Released

% APR 20, 2005
% in OBS mode artifact matrix now 'constant'  demeaned instead
%  of linear trend removal

% DEC 23, 2004
% updated (c)

% DEC 22, 2004
% fixed end artifact bug
% fixed copyright

% Dec16, 2004
% Added RAPCO
% re-written help

% Nov 4, 2004
% Updated brog bar
% for extra steps

% Oct 26, 2004
% Changed input to accept
% vector of event location
% instead of event type to
% allow for scripting.

% November 2008
% Made compatible with SPM8 structure
% added by yl

% March 2011
% Switching from meeg structure to object.
% cp


function [fitted_art, clean] = fmrib_pas(EEG,QRSevents,method,npc)

nargchk(3,4,nargin);
fs          = EEG.fsample;
channels    = EEG.cache.l2cor; % get list of EEG channels to use
nchannels   = numel(channels);
t_chunk     = EEG.cache.chunk_ii; % get chunk limits
samples     = diff(t_chunk)+1;
dat_o       = EEG(channels,t_chunk(1):t_chunk(2)); % data under processing

switch method

    %% Using the PCA style correction
    case 'obs'

        if nargin < 4
            error('Incorrect number of input arguments');
        end
        %init
        %-----
        pcs=npc-1;

        %standard delay between QRS peak and artifact (allen,1998)
        delay=round(0.21*fs);

        Gwindow=20;
        GHW=floor(Gwindow/2);
%         rcount=0;
%         firstplot=1;

        % memory allocation
        %------------------
%         avgdata=zeros(Gwindow,round(fs*2));
%         drift=zeros(1,round(fs*2)*Gwindow);
        fitted_art=zeros(nchannels,samples);
%         mPA=zeros(1,round(2*fs));
        peakplot=zeros(1,samples);
%         mEye=eye(nchannels);

        %Extract QRS events
        %------------------

        peakplot(QRSevents)=1;
        sh=zeros(1,delay);
        np=length(peakplot);
        peakplot=[sh peakplot(1:np-delay)];
        peakI=find(peakplot>0);
        peakcount=length(peakI);
        pa=peakcount;

        %make filter
        %-----------
        a=[0 0 1 1];
        f=[0 0.4/(fs/2) 0.9/(fs/2) 1];
        ord=round(3*fs/0.5);
        fwts=firls(ord,f,a);

        %Artifact Subtraction
        %---------------------

        % for the first window/2 points use arthemitic mean for averageing.
        % findg mean QRS peak-to-peak (R-to-R) interval
        %     for ch=1:channels
        for ch=1:nchannels % use the EEG channel list
            if ch==1
                %init waitbar
                %------------
                barth=5;
                barth_step=barth;
                Flag25=0;
                Flag50=0;
                Flag75=0;
                fprintf('\nPulse artifact subtraction in progress...Please wait!\n');
            end

            RR = diff(peakI);
            mRR=mean(RR);
            sRR=std(RR);
            PArange=round(1.25*(mRR+sRR)/2);
            %         PArange=round((mRR+sRR)/2); % half RR
            midP=PArange+1;

            if ch==1
                while ((QRSevents(pa-1)+2*PArange-1) > samples)
                    pa=pa-1;
                end
                steps=nchannels*pa;
            end

            eegchan = filtfilt(fwts,1,dat_o(ch,:));
            pcamat = zeros(pa-2,2*PArange+1);
%             dpcamat=pcamat;
            %         for p=2:pa
            for p=3:(pa-1) % sometimes 2 beats are not longer than PArange (pm 130606)
                pcamat(p-2,:)=eegchan(peakI(p)-PArange:peakI(p)+PArange);
            end
            pcamat=detrend(pcamat','constant')';
            meaneffect=mean(pcamat);
            dpcamat=detrend(pcamat,'constant');
            [apc,ascore,asvar]=pca_calc(dpcamat');

            papc=[meaneffect' ascore(:,1:pcs)];


            try
                pad_fit=papc*(papc\...
                    detrend(dat_o(ch,peakI(1)-PArange:...
                    peakI(1)+PArange)','constant'));

                fitted_art(ch,peakI(1)-PArange:peakI(1)+PArange)=...
                    pad_fit';
            catch
            end


            %-----------------------------------------------------
            for p=2:GHW+1
                pad_fit=papc*(papc\...
                    detrend(dat_o(ch,peakI(p)-PArange:...
                    peakI(p)+PArange)','constant'));

                fitted_art(ch,peakI(p)-PArange:peakI(p)+PArange)=...
                    pad_fit';

                %update bar
                %----------
                percentdone=floor(((ch-1)*pa+p)*100/steps);
                if floor(percentdone)>=barth
                    if percentdone>=25 & Flag25==0
                        fprintf('25%% ')
                        Flag25=1;
                    elseif percentdone>=50 & Flag50==0
                        fprintf('50%% ')
                        Flag50=1;
                    elseif percentdone>=75 & Flag75==0
                        fprintf('75%% ')
                        Flag75=1;
                    elseif percentdone==100
                        fprintf('100%%\n')
                    else
                        fprintf('.')
                    end

                    while barth<=percentdone
                        barth=barth+barth_step;
                    end
                    if barth>100
                        barth=100;
                    end
                end
            end

            %---------------- Processing of central data ---------------------
            %cycle through peak artifacts identified by peakplot
%             rcount=GHW;
            for p=GHW+2:peakcount-GHW-1
                PreP=ceil((peakI(p)-peakI(p-1))/2);
                PostP=ceil((peakI(p+1)-peakI(p))/2);
                if PreP > PArange
                    PreP=PArange;
                end
                if PostP > PArange
                    PostP=PArange;
                end

                pad_fit=papc(midP-PArange:midP+PArange,:)*(papc(midP-PArange:midP+PArange,:)\...
                    detrend(dat_o(ch,peakI(p)-PArange:...
                    peakI(p)+PArange)','constant'));

                fitted_art(ch,peakI(p)-PreP:peakI(p)+PostP) = ...
                    pad_fit(midP-PreP:midP+PostP)';



                %update bar
                %----------
                percentdone=floor(((ch-1)*pa+p)*100/steps);
                if floor(percentdone)>=barth
                    if percentdone>=25 & Flag25==0
                        fprintf('25%% ')
                        Flag25=1;
                    elseif percentdone>=50 & Flag50==0
                        fprintf('50%% ')
                        Flag50=1;
                    elseif percentdone>=75 & Flag75==0
                        fprintf('75%% ')
                        Flag75=1;
                    elseif percentdone==100
                        fprintf('100%%\n')
                    else
                        fprintf('.')
                    end

                    while barth<=percentdone
                        barth=barth+barth_step;
                    end
                    if barth>100
                        barth=100;
                    end
                end
            end

            %---------------- Processing of last GHW peaks------------------
            sectionPoints=samples-(peakI(peakcount-GHW)-PreP)+1;


            for p=peakcount-GHW:peakcount
                try
                    pad_fit=papc(midP-PArange:midP+PArange,:)*(papc(midP-PArange:midP+PArange,:)\...
                        detrend(EEG(ch,peakI(p)-PArange:...
                        peakI(p)+PArange)','constant'));

                    fitted_art(ch,peakI(p)-PreP:peakI(p)+PostP)=...
                        pad_fit(midP-PreP:midP+PostP)';


                catch
                end


                %update bar
                %----------
                percentdone=floor(((ch-1)*pa+p)*100/steps);
                if floor(percentdone)>=barth
                    if percentdone>=25 & Flag25==0
                        fprintf('25%% ')
                        Flag25=1;
                    elseif percentdone>=50 & Flag50==0
                        fprintf('50%% ')
                        Flag50=1;
                    elseif percentdone>=75 & Flag75==0
                        fprintf('75%% ')
                        Flag75=1;
                    elseif percentdone==100
                        fprintf('100%%\n')
                    else
                        fprintf('.')
                    end

                    while barth<=percentdone
                        barth=barth+barth_step;
                    end
                    if barth>100
                        barth=100;
                    end
                end
            end
        end
        if size(dat_o,2)<size(fitted_art,2)
            fitted_art(:,(size(dat_o,2)+1):end) = [];
        end
        clean = dat_o - fitted_art;

        %% Using the AAS style correction
    otherwise

        %init
        %-----
        %standard delay between QRS peak and artifact (allen,1998)
        delay=round(0.21*fs);

        Gwindow=20;
        GHW=floor(Gwindow/2);
        rcount=0;
%         firstplot=1;

        % memory allocation
        %------------------
        avgdata=zeros(nchannels,round(fs*2),Gwindow);
        drift=zeros(nchannels,round(fs*2)*Gwindow);
        clean=zeros(nchannels,samples);
        fitted_art=zeros(nchannels,samples);
        mPA=zeros(nchannels,round(2*fs));
        mMag=zeros(nchannels,Gwindow);
        peakplot=zeros(1,samples);
%         mEye=eye(nchannels);

        %Extract QRS events
        %------------------

        peakplot(QRSevents)=1;
        sh=zeros(1,delay);
        np=length(peakplot);
        peakplot=[sh peakplot(1:np-delay)];
        peakI=find(peakplot>0);
        peakcount=length(peakI);

        %Artifact Subtraction
        %---------------------

        % for the first window/2 points use arthemitic mean for averageing.
        % findg mean QRS peak-to-peak (R-to-R) interval
        for p=2:GHW+1
            RR(p-1)=peakI(p)-peakI(p-1);
        end
        mRR=mean(RR);
        sRR=std(RR);
        PArange=round((mRR+sRR)/2); % half RR

        t=(0:peakI(GHW+1)+PArange-1)/fs;

        % subtraction of low freq trend; could use detrend.m function instead
        for n=1:nchannels
            pdrift = polyfit(t,dat_o(n,1:peakI(GHW+1)+PArange),1);
            drift(n,1:peakI(GHW+1)+PArange)=polyval(pdrift,t);
        end

        dat_o(:,1:peakI(GHW+1)+PArange,1) = ...
                            dat_o(:,1:peakI(GHW+1)+PArange,1) - ...
                            drift(:,1:peakI(GHW+1)+PArange,1);

        for p=2:GHW+1
            avgdata(:,1:2*PArange+1,p-1)=...
                dat_o(:,peakI(p)-PArange:peakI(p)+PArange);
        end

        for chan =1:nchannels
            avgdata(chan,:,:)=detrend(squeeze(avgdata(chan,:,:)),'constant');
        end

        mPA(:,1:2*PArange+1)=mean(avgdata(:,1:1+2*PArange,1:GHW),3);

        if peakI(1) > PArange

            alphaV=...
                sum(detrend(dat_o(:,peakI(1)-PArange:peakI(1)+PArange)','constant')'...
                .*mPA(:,1:2*PArange+1),2)./...
                sum(mPA(:,1:2*PArange+1).*mPA(:,1:2*PArange+1),2);

            alphaD=diag(alphaV);

        else

            alphaV=sum(detrend(dat_o(:,1:peakI(1)+PArange)','constant')'...
                .*mPA(:,PArange-peakI(1)+2:2*PArange+1),2)./...
                sum(mPA(:,PArange-peakI(1)+2:2*PArange+1).*...
                mPA(:,PArange-peakI(1)+2:2*PArange+1),2);

            alphaD=diag(alphaV);

        end


        dat_o(:,1:peakI(GHW+1)+PArange) = ...
                                        dat_o(:,1:peakI(GHW+1)+PArange)+ ...
                                        drift(:,1:peakI(GHW+1)+PArange);

        %init waitbar
        %------------
        barth=5;
        barth_step=barth;
        Flag25=0;
        Flag50=0;
        Flag75=0;
        fprintf('\nPulse artifact subtraction in progress...Please wait!\n');


        if peakI(1) > PArange


            clean(:,peakI(1)-PArange:peakI(1)+PArange)=...
                dat_o(:,peakI(1)-PArange:peakI(1)+PArange)-alphaD*mPA(:,1:2*PArange+1);
            startpoints=peakI(1)-PArange-1;
            clean(:,1:startpoints)=dat_o(:,1:startpoints);
            fitted_art(:,peakI(1)-PArange:peakI(1)+PArange) = alphaD*mPA(:,1:2*PArange+1);
        else

            clean(:,1:peakI(1)+PArange)=...
                dat_o(:,1:peakI(1)+PArange)-alphaD*mPA(:,PArange-peakI(1)+2:2*PArange+1);
            fitted_art(:,1:peakI(1)+PArange)= alphaD*mPA(:,PArange-peakI(1)+2:2*PArange+1);

        end

        %update bar
        %----------
        percentdone=floor(1*100/peakcount);
        if floor(percentdone)>=barth
            if percentdone>=25 & Flag25==0
                fprintf('25%% ')
                Flag25=1;
            elseif percentdone>=50 & Flag50==0
                fprintf('50%% ')
                Flag50=1;
            elseif percentdone>=75 & Flag75==0
                fprintf('75%% ')
                Flag75=1;
            elseif percentdone==100
                fprintf('100%%\n')
            else
                fprintf('.')
            end

            while barth<=percentdone
                barth=barth+barth_step;
            end
            if barth>100
                barth=100;
            end
        end


        %-----------------------------------------------------
        for p=2:GHW+1

            alphaV=sum(dat_o(:,peakI(p)-PArange:peakI(p)+PArange).*...
                         mPA(:,1:2*PArange+1),2)./...
                         sum(mPA(:,1:2*PArange+1).*mPA(:,1:2*PArange+1),2);

            alphaD=diag(alphaV);

            clean(:,peakI(p)-PArange:peakI(p)+PArange)=...
                            dat_o(:,peakI(p)-PArange:peakI(p)+PArange)-...
                            alphaD*mPA(:,1:2*PArange+1);
            fitted_art(:,peakI(p)-PArange:peakI(p)+PArange) = alphaD*mPA(:,1:2*PArange+1);
            %update bar
            %----------
            percentdone=floor(p*100/(peakcount-1));
            if floor(percentdone)>=barth
                if percentdone>=25 & Flag25==0
                    fprintf('25%% ')
                    Flag25=1;
                elseif percentdone>=50 & Flag50==0
                    fprintf('50%% ')
                    Flag50=1;
                elseif percentdone>=75 & Flag75==0
                    fprintf('75%% ')
                    Flag75=1;
                elseif percentdone==100
                    fprintf('100%%\n')
                else
                    fprintf('.')
                end

                while barth<=percentdone
                    barth=barth+barth_step;
                end
                if barth>100
                    barth=100;
                end
            end
        end

        %---------------- Processing of central data ---------------------
        %cycle through peak artifacts identified by peakplot
        for p=GHW+2:peakcount-GHW-1
            PreP=ceil((peakI(p)-peakI(p-1))/2); % interval before peak artifact
            PostP=ceil((peakI(p+1)-peakI(p))/2);% interval after peak artifact
            sectionPoints=peakI(p+GHW)+PostP-(peakI(p-GHW)-PreP-1);

            % subtract drift
            t=(0:sectionPoints-1)/fs;
            for n=1:nchannels
                pdrift=...
                    polyfit(t,dat_o(n,peakI(p-GHW)-PreP:peakI(p+GHW)+PostP),1);
                drift(n,1:sectionPoints)=polyval(pdrift,t);
            end


            dat_o(:,peakI(p-GHW)-PreP:peakI(p+GHW)+PostP) = ...
                dat_o(:,peakI(p-GHW)-PreP:peakI(p+GHW)+PostP) - ...
                drift(:,1:sectionPoints);

            rr=1;

            for r=p-GHW:p+GHW
                mMag(:,rr)=mean(abs(dat_o(:,peakI(r)-PreP:peakI(r)+PostP)),2);
                rr=rr+1;
            end

            minMag=min(mMag,[],2);

            %only count those EEG data with no motion or line artifact by rejecting
            %excessively large variations.  Use channel 1 as ref (Cz).
            oldavgdata=avgdata;
            oldrcount=rcount;
            rcount=0;
            for r=p-GHW:p+GHW
                if (mean(abs(dat_o(1,peakI(r)-PreP:peakI(r)+PostP))) ...
                        < 4*minMag(1) ) ...
                        & (max(abs(dat_o(1,peakI(r)-PreP:peakI(r)+PostP)))...
                        < 2* 175)
                    avgdata(:,1:PreP+PostP+1,rcount+1)=...
                        dat_o(:,peakI(r)-PreP:peakI(r)+PostP);
                    rcount=rcount+1;
                end
            end

            for chan = 1:nchannels
                avgdata(chan,:,:) = detrend(squeeze(avgdata(chan,:,:)),'constant');
            end

            %     rcount
            %if more than 8 good data sections found then do averaging, else use
            %previous.
            if rcount > 5
                switch method
                    case 'gmean'
                        if p==peakcount-GHW-1
                            mPA=mean(avgdata(:,:,1:rcount),3);
                        elseif rcount==(2*GHW+1)
                            G=gausswin(rcount);
                            G=G/sum(G);
                            for r=1:rcount
                                avgdata(:,:,r)=avgdata(:,:,r)*G(r);
                            end
                            mPA=sum(avgdata(:,:,1:rcount),3);
                        else
                            mPA=median(avgdata(:,:,1:rcount),3);
                        end
                    case 'mean'
                        mPA=mean(avgdata(:,:,1:rcount),3);
                    case 'median'
                        mPA=median(avgdata(:,:,1:rcount),3);
                end
            else
                switch method
                    case 'mean'
                        mPA=mean(oldavgdata(:,:,1:oldrcount),3);
                    case 'median'
                        mPA=median(oldavgdata(:,:,1:oldrcount),3);
                    case 'gmean'
                        mPA=mean(oldavgdata(:,:,1:oldrcount),3);
                end
            end

            alphaV=sum(detrend(dat_o(:,peakI(p)-PreP:peakI(p)+PostP)','constant')'.*...
                mPA(:,1:PreP+PostP+1),2)./...
                sum(mPA(:,1:PreP+PostP+1).*mPA(:,1:PreP+PostP+1),2);

            alphaD=diag(alphaV);

            dat_o(:,peakI(p-GHW)-PreP:peakI(p+GHW)+PostP) = ...
                                        dat_o(:,peakI(p-GHW)-PreP:peakI(p+GHW)+PostP)+ ...
                                        drift(:,1:sectionPoints);

            clean(:,peakI(p)-PreP:peakI(p)+PostP)=...
                                dat_o(:,peakI(p)-PreP:peakI(p)+PostP)-...
                                alphaD*mPA(:,1:PreP+PostP+1);
            fitted_art(:,peakI(p)-PreP:peakI(p)+PostP)= alphaD*mPA(:,1:PreP+PostP+1);
            %update bar
            %----------
            percentdone=floor(p*100/(peakcount-1));
            if floor(percentdone)>=barth
                if percentdone>=25 & Flag25==0
                    fprintf('25%% ')
                    Flag25=1;
                elseif percentdone>=50 & Flag50==0
                    fprintf('50%% ')
                    Flag50=1;
                elseif percentdone>=75 & Flag75==0
                    fprintf('75%% ')
                    Flag75=1;
                elseif percentdone==100
                    fprintf('100%%\n')
                else
                    fprintf('.')
                end

                while barth<=percentdone
                    barth=barth+barth_step;
                end
                if barth>100
                    barth=100;
                end
            end
        end

        %---------------- Processing of last GHW peaks------------------
        sectionPoints=samples-(peakI(peakcount-GHW)-PreP)+1;

        if samples-peakI(peakcount) >= PostP

            alphaV=sum(detrend(dat_o(:,peakI(peakcount)-...
                        PreP:peakI(peakcount)+PostP)','constant')'.*...
                    mPA(:,1:PreP+PostP+1),2)./...
                    sum(mPA(:,1:PreP+PostP+1).*mPA(:,1:PreP+PostP+1),2);

            alphaD=diag(alphaV);

            clean(:,peakI(peakcount)-PreP:peakI(peakcount)+PostP)=...
                dat_o(:,peakI(peakcount)-PreP:peakI(peakcount)+PostP)- ...
                alphaD*mPA(:,1:PreP+PostP+1);
            fitted_art(:,peakI(peakcount)-PreP:peakI(peakcount)+PostP)= alphaD*mPA(:,1:PreP+PostP+1);

            endpoints=samples-peakI(peakcount)-PostP;

            clean(:,samples-endpoints+1:samples)=...
                dat_o(:,samples-endpoints+1:samples);

        else

            alphaV=sum(detrend(dat_o(:,peakI(peakcount)-PreP:samples)','constant')'.*...
                mPA(:,1:samples-(peakI(peakcount)-PreP)+1),2)./...
                sum(mPA(:,1:samples-(peakI(peakcount)-PreP)+1).*...
                mPA(:,1:samples-(peakI(peakcount)-PreP)+1),2);

            alphaD=diag(alphaV);

            clean(:,peakI(peakcount)-PreP:samples)=...
                dat_o(:,peakI(peakcount)-PreP:samples)-...
                alphaD*mPA(:,1:samples-(peakI(peakcount)-PreP)+1);
            fitted_art(:,peakI(peakcount)-PreP:samples)= alphaD*mPA(:,1:samples-(peakI(peakcount)-PreP)+1);
        end



        for p=peakcount-GHW:peakcount-1

            alphaV=sum(detrend(dat_o(:,peakI(p)-PreP:peakI(p)+PostP)','constant')'.*...
                mPA(:,1:PreP+PostP+1),2)./...
                sum(mPA(:,1:PreP+PostP+1).*mPA(:,1:PreP+PostP+1),2);

            alphaD=diag(alphaV);

            clean(:,peakI(p)-PreP:peakI(p)+PostP)=...
                dat_o(:,peakI(p)-PreP:peakI(p)+PostP)-...
                alphaD*mPA(:,1:PreP+PostP+1);
            fitted_art(:,peakI(p)-PreP:peakI(p)+PostP)= alphaD*mPA(:,1:PreP+PostP+1);

            %update bar
            %----------
            percentdone=floor(p*100/(peakcount-1));
            if floor(percentdone)>=barth
                if percentdone>=25 & Flag25==0
                    fprintf('25%% ')
                    Flag25=1;
                elseif percentdone>=50 & Flag50==0
                    fprintf('50%% ')
                    Flag50=1;
                elseif percentdone>=75 & Flag75==0
                    fprintf('75%% ')
                    Flag75=1;
                elseif percentdone==100
                    fprintf('100%%\n')
                else
                    fprintf('.')
                end

                while barth<=percentdone
                    barth=barth+barth_step;
                end
                if barth>100
                    barth=100;
                end
            end
        end

end
return;


function [EVec, Eload, EVal] = pca_calc(vecs)
%
%  PCA_CALC   Principal Component Analysis
%
%   [EVEC, ELOAD, EVAL] = PCA_CALC(X) takes a column-wise de-meaned
%       data matrix X of size n-by-m and calculates the n
%       Eigenvectors (EVEC of size n-by-n) of the data covariance
%       matrix, the factor loadings (ELOAD of size n-by-m) and the
%       corresponding Eigenvalues (EVAL of size n-by-1).
%
%
% Author: Rami K. Niazy, FMRIB Centre, University of Oxford
%
% Copyright (c) 2004 University of Oxford.
%

[m,n]=size(vecs);
[Us,S,EVec] = svd(vecs,0);

if m == 1
    S = S(1);
else
    S = diag(S);
end
Eload = Us .* repmat(S',m,1);
S = S ./ sqrt(m-1);
if m <= n
    S(m:n,1) = 0;
    S(:,m:n) = 0;
end
EVal = S.^2;