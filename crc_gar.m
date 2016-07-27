function [file2] = crc_gar(prefile)

% FORMAT crc_gar(prefile)
% Gradient artefact rejection routine
% Input:
% prefile - EEG data file to correct, it can either be an imported SPM
%           format data set, or the raw BrainProducts data file.
%__________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by Y. Leclercq & C. Phillips, 2008.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

crcdef= crc_get_defaults('gar');

% Definition of the file to treat
if nargin<1 || (size(prefile,1)==1 && ~exist(prefile,'file'))
    prefile = spm_select(Inf, 'any', 'Select imported EEG file','' ,pwd, ...
        '\.[MVEemv][HADdha][DTFfdt]');
end

for fileii=1:size(prefile,1)
%     try
        D = crc_eeg_load(deblank(prefile(fileii,:)));
        file = fullfile(D.path,D.fname);
        if ~exist(file,'file'),
            warning('No valid file selected')
            return
        end
        
        % Set parameters
        %---------------
        % Sampling frequency of the output file
        fs = crcdef.output_fs;
        % Number of scan used to compute the average
        Nscav = crcdef.Nsc_aver;
        % Number of scan ignored at the beginning of the file.
        Nscig = crcdef.Nsc_skipped_be;
        %Define the channels to treat
        totreat = 1:D.nchannels;
        % To fix by something more appropriate? E.g. setdiff(1:D.nchannels,ecgchannels(D));
        
        % Definition of the subsampling factor
        n = round(D.fsample/fs);
        newfsample = round(D.fsample/n);
        %-Test for the presence of Signal Processing Matlab toolbox
        if n>1
            flag_TBX = license('checkout','signal_toolbox');
            if ~flag_TBX
                disp(['warning: using homemade resampling routine ' ...
                    'as signal toolbox is not available.']);
            end
        end
        
        ev = events(D);
        if crcdef.UseScanMrk
            
            % Create event code
            evts = round([ev(abs([ev.value])==crcdef.ScanMrk1).time]*D.fsample); % CP, why not use true tb index???
            if isempty(evts)
                evts = round([ev(abs([ev.value])==crcdef.ScanMrk2).time]*D.fsample); % CP, why not use true tb index???
            end
            
            if length(evts)<crcdef.Nsc_aver
                disp(sprintf(['ERROR: Cannot find sufficient scan '...
                    'trigger (%d or %d) in the events field, correction'...
                    'cannot be performed'],crcdef.ScanMrk1,crcdef.ScanMrk2));
                return
            end
            
            %Check that all evts are in databounds
            checkoutofbounds = find(evts(:) < D.nsamples);
            %             checkoutofbounds = find(evts(:) < D.Nsamples);
            evts = evts(1:checkoutofbounds(end));
            
            evts = evts(1:crcdef.MrkStep:end);
            %stephr = round(mean(evts(2:end)-evts(1:end-1)));
            stephr = round(median(evts(2:end)-evts(1:end-1)));
            steplr = round(stephr/n);
            
            beg = evts(1+Nscig); % skip first few Nscig scans
            nd = evts(end)+stephr; % catch the last scan too
            
        else
            if crcdef.AutoChk
                Threshold = crcdef.Threshold;
                Step    = crcdef.DetStep;
                Chan = crcdef.DetChan;
                
                % Determination of beg
                criterion = 1;
                dbchk = 1;
                while criterion
                    if mean(abs(D(Chan,dbchk:dbchk+Step*D.fsample)))> Threshold;
                        criterion = 0;
                        beg = dbchk + Step*D.fsample/2 + Nscig*crcdef.TR*D.fsample;
                    end
                    
                    dbchk=dbchk+Step*D.fsample;
                    
                end
                
                % Determination of nd
                criterion = 1;
                ndchk = D.nsamples;
                while criterion
                    if mean(abs(D(Chan,ndchk-Step*D.fsample:ndchk)))> Threshold;
                        criterion = 0;
                        nd = ndchk - Step*D.fsample;
                    end
                    
                    ndchk = ndchk - Step*D.fsample;
                    
                end
                
            else
                beg = crcdef.beg*newfsample*n;
                nd = crcdef.nd*newfsample*n;
            end
            stephr = round(crcdef.TR*newfsample*n);
            steplr = round(stephr/n);
            
        end
        
        nTR = floor((nd-beg)/stephr);
        % WARNING could lead to the loss of the last TR, if nd bound falls 
        % a bit short...
        if nTR<=Nscav
            error('Not enough TR''s for correction!')
        end
        
        % Create new object,
        % with new sample rate and data type set to float
        file2 = fullfile(D.path ,[crcdef.prefix D.fname]);
        Dtmp = D;
        if isempty(strfind(lower(dtype(D)),'float32')) % not already in floats
            try % using modified dtype function from meeg object
                Dtmp = dtype(Dtmp,'float32-LE');
            catch % or use ugly structure-fix
                Dts = struct(Dtmp);
                Dts.data.y.dtype = 'float32-LE';
                Dts.data.datatype = 'float32-LE';
                Dtmp = meeg(Dts);
            end
        end
        Do = clone(Dtmp,file2,[Dtmp.nchannels nTR*steplr 1]);
        clear Dtmp
        Do = fsample(Do,newfsample);
        
        % Update structure and save file
        try
            hour = D.info.hour;
            date = D.info.date;
            
            newbeg = datevec(datenum([date hour]) + ...
                datenum([0 0 0 crc_time_converts(evts(1+Nscig)/D.fsample/n)]));
            
            Do.info.hour = newbeg(4:6);
            Do.info.date = newbeg(1:3);
        catch
        end
        
        % New events
        new_ev = ev;
        if ~isempty(new_ev)
            for mm=1:length(new_ev)
                new_ev(mm).time = new_ev(mm).time - beg/(n*newfsample) ;
            end
            
            % Look for events beginning before the beginning of the new file or
            % after the end of the new file and suppressing them.
            tosuppress = find([new_ev.time]<1);
            new_ev = new_ev(tosuppress(end):end);
            tosuppress = find([new_ev.time]>nd/(n*newfsample));
            if tosuppress, new_ev = new_ev(1:(tosuppress(1)-1)); end
            % Put "new" events back in place
            Do = events(Do,1,new_ev);
        end
        
        h = waitbar(0,'Please wait...');
        set(h,'Name',D.fnamedat)
        
        nbeaf = floor(Nscav/2); % half #TR to do average
        ch1 = true; %ch2 = true; % flags
        
        for ii=1:nTR
            prog = round(100*ii/nTR);
            string = ['Please wait... ' num2str(prog) ' %'];
            waitbar(prog/100,h,string)
            c_ind1 = mod(ii-1,2*nbeaf+1) + 1 ;
            
            if ch1
                % first nbeaf+1 TR's same estimate of artefact
                data = D(totreat,beg+(1:(stephr*(2*nbeaf+1))));
                cube = reshape(data,length(totreat),stephr,nbeaf*2+1);
                avg_artef_hr = mean(cube,3);
                if n>1
                    if flag_TBX % Signal Proc. Toolbox
                        avg_artef = resample(avg_artef_hr', 1, n)';
                    else
                        avg_artef = sthsub(avg_artef_hr,n);
                    end
                else
                    avg_artef = avg_artef_hr;
                end
                ch1 = false;
                
            elseif ii>nbeaf+1 %&& ii<=nTR-nbeaf
                % next TR's, update the estimate of artefact
                % load only the new TR for averaging
                if ii<=nTR-nbeaf
                    c_ind2 = mod(ii-1+nbeaf,2*nbeaf+1) + 1 ;
                    cube(:,:,c_ind2) = D(totreat,beg+(ii+nbeaf-1)*stephr+(1:stephr));
                    avg_artef_hr = mean(cube,3);
                    if n>1
                        if flag_TBX % Signal Proc. Toolbox
                            avg_artef = resample(avg_artef_hr', 1, n)';
                        else
                            avg_artef = sthsub(avg_artef_hr,n);
                        end
                    else
                        avg_artef = avg_artef_hr;
                    end
                end
            end
            if n>1
                if flag_TBX % Signal Proc. Toolbox
                    dat_ii = resample(squeeze(cube(:,:,c_ind1))', 1, n)';
                else
                    dat_ii = sthsub(squeeze(cube(:,:,c_ind1)),n);
                end
                
            else
                dat_ii = squeeze(cube(:,:,c_ind1));
            end
            
            Do(totreat,((ii-1)*steplr+1):min(ii*steplr,Do.nsamples)) = ...
                dat_ii - avg_artef;
            
        end
        % save version number of routine
        [v,r] = crc_fasst_utils('Ver',mfilename);
        Do.info.ver_gar = struct('v_nr',v,'rel',r);        

        save(Do)
        close(h)
%     catch
%         disp([' !!! ERROR while processing file:' deblank(prefile(fileii,:)),' !!!']);
%     end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DECLARATION OF SUBFUNCTION
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function newdata = sthsub(data,n)

% Function to smooth and subsample data by a factor n
% Smoothing is done using an hamming window
%
% Usage newdata = sthsub(data,n);
%
% Where data are the data to subsample/smooth and have the dimension Nch x
% length. Where Nch is the number of channels and length the number of
% samples for each channel.
%
% n is the subsampling factor.
%

window = hamming(n)/sum(hamming(n));
lgdata = size(data,2);
Nch = size(data,1);
nbrewin = floor(lgdata/n);
bigwin = reshape(window*ones(1,nbrewin),1,n*nbrewin);

biggywin = ones(Nch,1)*bigwin;

tmp0 = (biggywin.*data(:,1:size(biggywin,2)));
tmp1 = reshape(tmp0',n,nbrewin*Nch);
tmp2 = sum(tmp1,1);

newdata = reshape(tmp2,nbrewin,Nch)';

return
