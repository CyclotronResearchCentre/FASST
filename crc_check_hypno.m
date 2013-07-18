function [n_empty,l_empty] = crc_check_hypno(score, scorer, i_win, n_disp)

% FORMAT [n_empty,l_empty] = crc_check_hypno(score, scorer, i_win, n_disp)
% Checking the hypnogram of a scorer for windows without score.
% 
% INPUT:
% - score       : score cell array
% - scorer      : index of scorer (1 by default)
% - i_win       : index of current window (1 by default)
% - n_disp      : number of indexes of empty windows printed at the prompt (10 by default)
%
% OUTPUT:
% - n_empty     : index of the next window without score
% - l_empty     : list of windows without scores
%________________________________________________________________________
% Copyright (C) 2013 Cyclotron Research Centre

% Written by C. Phillips, 2013.
% Cyclotron Research Centre, University of Liège, Belgium
% $Id$

if nargin<4
    n_disp = 10 ;
end
if nargin<3 || isempty(i_win)
    i_win = 1 ;
end
if nargin<2 || isempty(scorer)
    scorer = 1 ;
end

sc = score{1,scorer};
l_empty = find(isnan(sc));
tmp = find(l_empty>=i_win);
if ~isempty(tmp)
    n_empty = l_empty(tmp(1));
else
    fprintf('\n\t NO empty window after widow #%d.\n',i_win);
end

if n_disp>0
    l_next = find(l_empty>=i_win);
    en_disp = min(n_disp,length(l_next));
    fprintf('List of %d (out of %d) window(s) without score following #%d:\n',...
        en_disp,length(l_empty),i_win)
    for ii=1:en_disp
        fprintf('\t window #%d\n',l_empty(l_next(ii)))
    end
    fprintf('\n')
end

end 