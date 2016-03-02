function CM = crc_compareScore(sc1,sc2,flags)
%
% Compare 2 sleep scores, i.e. hypnograms.
% sc1 is assumed as the reference and sc2 is compared to it.
%
%_______________________________________________________________________
% Copyright (C) 2014 Cyclotron Research Centre

% Written by C. Phillips, 2014
% Cyclotron Research Centre, University of Liege, Belgium
% check flags

%% Checks
if numel(sc1)~=numel(sc2)
    error('CRC:compSc','Scores provided have different length');
end

%% Calculations
[CM, cOrder] = confusionmat(sc1,sc2);

figure;
suptitle(['', char(flags.name)])
subplot(211)
plot([sc1' sc2']),ylim([-0.5 6.6]), legend('FASST','Aseega'),
subplot(212),
plot([sc1 - sc2]),ylim([-6.6 6.6]), legend('sc1 - sc2'),
end


