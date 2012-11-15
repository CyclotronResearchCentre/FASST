function crc_run_feedback(D, opt)
%CRC_RUN_FEEDBACK gives simple feedback after event detections
% Use as:
%   crc_run_feedback(D, opt)
% where opt has
%   .pad = padding (s)
%   .scl = scaling (in uV)
%   .met = methods ('interactive', GUI format, or 'print', print to file)
%   .wavtyp = 'slowwaves' or 'spindles'
%_______________________________________________________________________
% Copyright (C) 2010 Cyclotron Research Centre

% 10/11/16 avoid out of range errors
% 10/11/16 append comment
% 10/11/15 created

% Giovanni Piantoni
% $Id$


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check opt(ions)
if nargin == 1 || ~isfield(opt, 'wavtyp')
  disp('please, specify a wave type. Use opt.wavtyp = ''spindles'' or ''slowwaves''')
  return % don't throw an error, too cruel. 
end

if ~isfield(opt, 'pad'); opt.pad = crc_get_defaults('fbck.pad'); end % padding in s
if ~isfield(opt, 'scl'); opt.scl = crc_get_defaults('fbck.scl'); end % scaling
if ~isfield(opt, 'col'); opt.col = crc_get_defaults('fbck.col'); end % n of columns for 'print'

if ~isfield(opt, 'met'); opt.met = 'interactive'; end % interactive plot

if ~isfield(opt, 'comm'); opt.comm = ''; end % append comment

%%%%%%%%%%%%%%
% load data
if ischar(D)
  crc_eeg_load(D)
end

if isstruct(D)
  D = meeg(D);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare detection pnts
switch opt.wavtyp
  
  %%%%%%%%%%%%%%
  case 'spindles'
    
    if isfield(D.CRC, 'spindles')
      allbeg = D.CRC.spindles.bounds(:,1);
      allend = D.CRC.spindles.bounds(:,2);
      
      for k = 1 : size(D.CRC.spindles.maxelectrode, 1)
        [dum, allchn{k}] = intersect(chanlabels(D), D.CRC.spindles.maxelectrode(k,:));
      end
      
    else
      disp('could not find spindles in the file')
      return
    end
    
    %%%%%%%%%%%%%%
  case 'slowwaves'
    if isfield(D.CRC, 'SW')
      eegchn = meegchannels(D);
      for k = 1:numel(D.CRC.SW.SW)
        allbeg(k,1) = D.CRC.SW.SW(k).negmax_tp;
        allend(k,1) = D.CRC.SW.SW(k).posmax_tp;
        allchn{k} = eegchn(D.CRC.SW.SW(k).channels);
      end
    else
      disp('could not find slow waves in the file')
      return
    end
    
  otherwise
    disp('Unknown method')
    return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nevt  = size(allbeg,1); % n events

switch opt.met
  
  %%%%%%%%%%%%%%
  case 'interactive'
    
    for k = 1 : nevt
      detbeg = allbeg(k,1);
      detend = allend(k,1);
      chn    = allchn{k};
      h = figure;
      subfeedback(D, detbeg, detend, chn, opt)
      title(sprintf('%s n. %1.f/%1.f', opt.wavtyp, k, nevt))
      waitfor(h)
    end
    
  %%%%%%%%%%%%%%
  case 'print'
    
    nrow  = ceil( nevt / opt.col);
    
    h = figure;
    for k = 1 : nevt
      subplot( nrow, opt.col, k)
      detbeg = allbeg(k,1);
      detend = allend(k,1);
      chn    = allchn{k};
      subfeedback(D, detbeg, detend, chn, opt)
      set(gca, 'xtick', [], 'ytick', [])
%      title(sprintf('%1.f/%1.f', k, nevt))
    end
    
    %%%%%%%%%%%%%%
    % fix scaling
    pos = get(h, 'pos');
    pos(4) = pos(3) / opt.col * nrow; % every subplot should be square
    set(h, 'pos', pos)
    
    %%%%%%%%%%%%%%
    % print to png
    Dname    = fname(D); % to remove .mat at the end
    filename = sprintf('%s%sfbck_%s_%s_%s.png', path(D), filesep, Dname(1:end-4), opt.wavtyp, opt.comm);
    print(h, filename, '-dpng', '-r600')
    close(h)
  otherwise
    disp('Unknown method')
    return
    
end

function subfeedback(D, detbeg, detend, chn, opt)
% do the actual plotting for one event

padbeg = detbeg - opt.pad * fsample(D);
padend = detend + opt.pad * fsample(D);

% avoid out of range errors
if padbeg <= 0;           padbeg = 1; end
if padend >= nsamples(D); padend = nsamples(D); end

Tdet   = (detbeg:detend) / fsample(D);
Tpad   = (padbeg:padend) / fsample(D);

paddat = D(chn, padbeg:padend, 1);
detdat = D(chn, detbeg:detend, 1);

plot(Tpad, paddat, 'k');
hold on
plot(Tdet, detdat, 'r')

plot(Tpad, zeros(size(Tpad)), '--') % zero-line

%%%%%%%%%%%%%%
xlim([Tpad(1) Tpad(end)])

% y-axis is defined as starting from 10uV below lowest point,
% height is equal to scaling (opt.scl)
% optimized for slow waves
offset = min(min(detdat)) - 10;
ylim([0 opt.scl] + offset)

