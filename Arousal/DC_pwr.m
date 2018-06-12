
function D = DC_pwr(D, config) 

% FORMAT DC_pwr(args)
% It computes the power spectrum (after bad channels
% detection) in every 1-s epoch, in the following
% frequency bands:
% * delta : 0.5-3Hz
% * theta : 4-7Hz
% * alpha : 8-11Hz
% * sigma : 12-14Hz
% * beta : 16-30Hz
% * pi : 0.5-30Hz
% It computes the mean power spectral density over all channels 
%
% INPUT
%       .file   - data file (.mat files)
%__________________________________________________________________
% Copyright (C) 2014 Cyclotron Research Centre

% Written by D. Coppieters 't Wallant, 2014.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$
% ----------------------

% load parameters
channels    =   config.channels;
winsize     =   config.winsize;
fs      =   fsample(D);
nspl    =   nsamples(D); 
NofW    =   ceil(nspl/(fs*winsize));
Time    =   ceil(nspl/fs);

delta_s =   zeros(1,Time);
theta_s =   zeros(1,Time);
alpha_s =   zeros(1,Time);
spindle_s = zeros(1,Time);
beta_s  =   zeros(1,Time);
    % to take into account the bad channels previously detected
    badchannels = cell(NofW,1);
%     for ibc = 1 : NofW
% 	badchannels{ibc} = union(D.CRC.DC.badchannels.chan_defaillant{ibc},D.CRC.DC.badchannels.chan_incoherent{ibc});
%     end
    h = waitbar(0,'Wait during the computation of power spectrum values ...');
    %%% new version
    for win = 1 : NofW
         if ~isempty(badchannels{win})
            goodchan =  setdiff(channels,badchannels{win});
        else 
            goodchan =  channels;
        end
        MeanEEG = mean(D(goodchan,(win-1)*winsize*fs+1:min(win*winsize*fs,nspl)));
        for kk = 1 : ceil(length(MeanEEG)/fs)
            small  =   MeanEEG(max((kk-1.2)*fs,1):min((kk+0.2)*fs,length(MeanEEG)));
            if length(small)>fs/2
                [P,F]   =   pmtm(small,[],0.1:0.5:30,fs);   
                delta_s(kk + (win-1)*winsize) =  mean(pmtm(small,[],0.5:0.1:4,fs));
                theta_s(kk + (win-1)*winsize) =  mean(pmtm(small,[],4:0.1:8,fs));
                alpha_s(kk + (win-1)*winsize) =  mean(pmtm(small,[],8:0.1:12,fs));
                spindle_s(kk + (win-1)*winsize) = mean(pmtm(small,[],12:0.1:16,fs));
                beta_s(kk + (win-1)*winsize) = mean(pmtm(small,[],16:0.1:30,fs));
            else
                delta_s(kk + (win-1)*winsize) = delta_s(end);
                theta_s(kk + (win-1)*winsize) = theta_s(end);
                alpha_s(kk + (win-1)*winsize) = alpha_s(end);
                spindle_s(kk + (win-1)*winsize) = spindle_s(end);
                beta_s(kk + (win-1)*winsize) = beta_s(end);
             end
        end
    end 
    D.CRC.DC.power.delta_s = delta_s;
    D.CRC.DC.power.theta_s = theta_s;
    D.CRC.DC.power.alpha_s = alpha_s;
    D.CRC.DC.power.spindle_s = spindle_s;
    D.CRC.DC.power.beta_s = beta_s;
    save(D); 
    close(h);
    fprintf('* power spectrum evaluated in six frequency bands: \n \b - slow wave activity (0.5-2Hz) - delta (3-5Hz) \n \b - theta (5-7Hz) \n \b - alpha (8-13Hz) \n \b - sigma (11-16Hz) \n \b - beta (16-30Hz) \n \b total power (0.5 - 30Hz)\n \b')
end

