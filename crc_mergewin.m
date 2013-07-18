function mw_excl = crc_mergewin(w_excl)

%%
% FORMAT mw_excl = crc_mergewin(w_excl)
%
% Function to merge a table of time windows. Any overlapping windows are
% merged together, producing a shorter list of the union of all the
% windows.
%
% INPUT
% - w_excl  : original time windows to merge, Nx2 matrix of [beg end]
%
% OUTPUT
% - mw_excl : merged time windows, Mx2 matrix of [beg end] with M<=N
%_______________________________________________________________________
% Copyright (C) 2013 Cyclotron Research Centre

% Written by C. Phillips, 2013.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

mw_excl = zeros(1,2)-Inf;
while ~isempty(w_excl)
    % treat windows from the beginning, i.e. look at start point
    [be,be_i] =  min(w_excl(:,1));
    be(2) = w_excl(be_i,2);
    % remove current window in all cases
    w_excl(be_i,:) = [];
    
    if size(mw_excl,1)==1 || be(2)>max(mw_excl(:,2))
        % if it is the 1st window or it goes beyond the end of existing 
        % new windows, then add it
        l_intersect = find( be(2)>w_excl(:,1) & be(2)<w_excl(:,2) );
        if numel(l_intersect)
            % check if there is any overlap with another window
            % and possibly use the end of that one
            [be(2),rem_i] = max(w_excl(l_intersect,2));
            % remove 2nd overlapping window
            w_excl(rem_i,:) = [];
        end
        % add this new window to list.
        mw_excl = [mw_excl ; be]; %#ok<AGROW>
    end
end
mw_excl(1,:) = []; % remove 1st initial line

end

%%
% % Testing examples
% w_excl = [ 1 5 ; 2 6 ; 3 8 ]
% w_excl = [ 1 8 ; 5 6 ; 10 12 ; 2 3 ]
% mw_excl

