
% Little script that lets me export the updated files, using the list of
% updated file from the text log-file.
%__________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by Y. Leclercq & C. Phillips, 2008.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

Dlog = spm_select(1,'.*\.txt','Select update log-file');
dir_target = spm_str_manip(Dlog,'h');
dir_orig = spm_str_manip(which('crc_main'),'h');

file = textread(Dlog,'%s','delimiter','\n','whitespace','');
Nl = length(file);

txt_Target = {'Modified : ','Added : '};
Ntxt_Tar = length(txt_Target);
for jj=1:Ntxt_Tar
    txt_target = txt_Target{jj};
    Ntxt_tar = length(txt_target);
    for ii=Nl:-1:1
        if length(file{ii})>=Ntxt_tar & strcmp(txt_target,file{ii}(1:Ntxt_tar))
            f2move = file{ii}(Ntxt_tar+2:end);
            f2move_full = fullfile(dir_orig,file{ii}(Ntxt_tar+2:end));
            f_dest = fullfile(dir_target,f2move);
            if exist(f2move_full,'file')==2
                try
                    copyfile(f2move,f_dest)
                catch
                    mkdir(spm_str_manip(f_dest,'h'))
                    copyfile(f2move,f_dest)
                end
            end
        end
    end
end
clear Dlog dir_target dir_orig file Nl txt_target Ntxt_target f2move ...
        f2move_full f_dest ii jj txt_Target Ntxt_Tar Ntxt_tar
