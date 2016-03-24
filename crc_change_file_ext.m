function out = crc_change_file_ext(fn_in,ext_new)
%
% Function to change the extension of a file.
% If possible it avoids using the 'movefile' command from Matlab as this is
% super slow but rather relies on built-in java function.
%
% FORMAT
%   out = crc_change_file_ext(fn_in,new_ext)
%
% INPUT:
% - fn_in :   input filename. If no path is provided, the function assumes
%             it is located inthe current working directory
% - ext_new : new extension for the input filename.
%
% OUTPUT:
% - out : flag indicating if it worked (>0) or not (0). The exact value
%         depends on how the extension was changed
%           . 1, if java was sucessfully used
%           . 2, if Matlab's movefile was used
%__________________________________________________________________________
% Copyright (C) 2016 Cyclotron Research Centre

% Written by C. Phillips, 2016.
% Cyclotron Research Centre, University of Liege, Belgium

% Check if java is available
persistent flag_java

if isempty(flag_java)
    flag_java = usejava('jvm');
end

% Create new filename
[pth_orig,fn_orig,~] = fileparts(fn_in);
fn_out = fullfile(pth_orig,fn_orig,ext_new);

try
    if flag_java
        out = 1;
        java.io.File(fn_in).renameTo(java.io.File(fn_out));
    else
        out = 2;
        movefile(fn_in,fn_out)
    end
catch err
    out = 0;
    display(err);
end

end
