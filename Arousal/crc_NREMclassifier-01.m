function D = crc_NREMclassifier(D,channels,out_spect)

% from D
filename = fname(D);
fs = fsample(D);
nspl = nsamples(D);
winsize = 30;
NofW = floor(nspl/(fs*winsize));
ne = numel(channels);

% parameters
band_delta = [0.5 3.5];
band_alpha = [8 10];
band_sigma = [12 15];

% load artefacts detected earlier in order to take them into account in the
% power spectrum computation
artspl = [];
if isfield(D,'CRC') && isfield(D.CRC,'DC') && isfield(D.CRC.DC,'shortartf')
    art = union_several(D.CRC.DC.shortartf.badchan,D.CRC.DC.shortartf.popping);
%     art = union_several(D.CRC.DC.shortartf.badchan,D.CRC.DC.shortartf.movement,D.CRC.DC.shortartf.arousal,D.CRC.DC.shortartf.popping,D.CRC.DC.shortartf.artefact);
    % ERROR with one person: for whatever reason, there is a 0 included in D.CRC.DC.shortartf.movement
    check = find(art == 0);
    if ~isempty(check)
        art(check) = [];
    end
    artspl = zeros(nspl,1);
    for iart = 1: numel(art)
        artspl((art(iart)-1)*fs+1 : art(iart)*fs) = (art(iart)-1)*fs+1 : art(iart)*fs;
    end
end

% load power spectrum values
try 
    load(fullfile(out_spect, ['spindle_sp_', filename]))
    load(fullfile(out_spect, ['spindle_fspect_', filename]))    
    % compute the power spectrum for all 'channels' derivation in D
    % remove artifacts if they were detected from crc_artifact method 
    % the power spectrum computation used a multi-taper approach (fct: pmtm)
    % with 50% of overlap between 4s-time windows.
    if size(Spect,3) == ne
        fprintf(1,'Power spectrum loaded \n')
    else
        fprintf(1,'Power spectrum computation... \n')
        Spect = zeros(NofW,numel(0.5:0.25:30),ne);
        for ie = 1:ne 
            fprintf(1,'Channel %d processed.. \n', ie)
            sf = D(channels(ie),:);
              for iw = 1:NofW
                subwin = (iw-1)*winsize*fs+1:iw*winsize*fs;
                clean_epoch = setdiff(subwin,artspl);
                if ~isempty(clean_epoch) && numel(clean_epoch)>winsize*fs/2
                    [Spect(iw,:,ie) fspect] = pmtm(sf(clean_epoch),4,0.5:0.25:30,fs);
                end
            end
        end
        save(fullfile(out_spect,['spindle_sp_' fname(D)]),'Spect','-v7.3')
        save(fullfile(out_spect,['spindle_fspect_' fname(D)]),'fspect','-v7.3')  
    end
catch
    % power spectrum computation
    Spect = zeros(NofW,numel(0.5:0.25:30),ne);
    fprintf(1,'Power spectrum computation... \n')
    for ie = 1:ne 
        fprintf(1,'Channel %d processed.. \n', ie)
        sf = D(channels(ie),:);
        for iw = 1:NofW
            subwin = (iw-1)*winsize*fs+1:iw*winsize*fs;
            clean_epoch = setdiff(subwin,artspl);
            if ~isempty(clean_epoch) && numel(clean_epoch)>winsize*fs/2
                [Spect(iw,:,ie) fspect] = pmtm(sf(clean_epoch),4,0.5:0.25:30,fs);
            end
        end
    end
    save(fullfile(out_spect,['spindle_sp_' fname(D)]),'Spect','-v7.3')
    save(fullfile(out_spect,['spindle_fspect_' fname(D)]),'fspect','-v7.3')
end

% initialization
weights = zeros(size(Spect,1),1);
weightd = zeros(size(Spect,1),1);
weighttot= zeros(size(Spect,1),1);
nrem_est = cell(ne,1);

for ie = 1 : ne
    for iw = 1 : size(Spect,1)
        weighttot(iw) = sum(Spect(iw,:,ie));
        weights(iw) = sum(Spect(iw,fspect>band_sigma(1) & fspect<band_sigma(2),ie))/sum(Spect(iw,(fspect>band_alpha(1) & fspect<band_alpha(2))|(fspect>20 & fspect<30),ie));
        weightd(iw) = sum(Spect(iw,fspect>band_delta(1) & fspect<band_delta(2),ie))/sum(Spect(iw,(fspect>band_alpha(1) & fspect<band_alpha(2))|(fspect>20 & fspect<30),ie)); 
    end
    lim_art = median(weighttot)+3*std(weighttot); 
    winart = find(weighttot>lim_art);
    winart = union(winart,unique(ceil(art/winsize)));

    sigma_act = setdiff(weights(weights>0&~isnan(weights)),winart);
    delta_act = setdiff(weightd(weightd>0&~isnan(weightd)),winart);
    objs = gmdistribution.fit(sigma_act(:),2,'Replicates',10); 
    objd = gmdistribution.fit(delta_act(:),2,'Replicates',10);

    [gum act_sig] = max(objs.mu);
    [gum act_sws] = max(objd.mu);
    sigm = cluster(objs,weights)==act_sig;
    sws = cluster(objd,weightd)==act_sws;
    nrem_est{ie} = setdiff(find(sigm|sws),winart);
end
  D.CRC.DC.classREM.Fz = nrem_est{1};
  D.CRC.DC.classREM.Cz = nrem_est{2};
  D.CRC.DC.classREM.Pz = nrem_est{3};
  save(D); 

