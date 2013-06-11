function [mtch,perc_match] = crc_hypnocompare(h_ax, v12, score, handles)

% FORMAT [mtch,perc_match] = crc_hypnocompare(h_ax, v12, score, handles)
% Estimating the merged hypnogram and plotting it.
% 
% OUTPUT:
% - mtch        : vector of difference between the 2 hypnograms
% - perc_mtch   : percentage of match of the 2 histogram based on (1) the
%                 total length and (2) what was scored by the 1st scorer
%________________________________________________________________________
% Copyright (C) 2013 Cyclotron Research Centre

% Written by C. Phillips, 2013.
% Cyclotron Research Centre, University of Liège, Belgium
% $Id$

if length(score{1,v12(1)})==length(score{1,v12(2)})

    mtch = score{1,v12(1)}-score{1,v12(2)};
    matchvect = mtch;
    matchvect(mtch ~= 0) = 7;
    perc_match(1) = sum(mtch == 0)/length(score{1,v12(2)})*100;
    perc_match(2) = sum(mtch == 0)/sum(~isnan(score{1,v12(1)}))*100;

    cla(h_ax(1))
    set(h_ax(2),'CurrentAxes',h_ax(1))
    crc_hypnoplot(h_ax(1),...
        handles,...
        score{3,1},...
        matchvect,...
        'Merge')

    labl{1}=' ';
    labl{2}='No match';
    labl{3}=' ';
    labl{4}=' ';
    labl{5}=' ';
    labl{6}=' ';
    labl{7}=' ';
    labl{8}='Match';

    % Display Wake state (st0 => 7)
    set(h_ax(1),'YTickLabel',labl);
else
    cla(h_ax(1))
    errordlg({'The window size of the score are not the same. They cannot be compared'},'Error')
    mtch = [];
    perc_match = 0;
end

fprintf('\n\nMatching between scores: ')
fprintf('\n\t%d%% of total\t , \n\t%d%% of scored by %s\n', ...
    round(perc_match(1)),round(perc_match(2)),score{2,v12(1)})

end
