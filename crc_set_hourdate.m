function crc_set_hourdate(files)
% (Re-)set the "real world" recording time and date of BP data, if the
% information is still available in the .vhdr header  file
% FORMAT D = crc_set_hourdate(files)
%
% files       - filenames of SPM8-meeg format data file
%_______________________________________________________________________
%
% The routine simply reloads the header and marker files (.vhdr/.vmrk) with
% the same filename as the .mat file, then stores again the information in
% the .info field.
% This can be useful, for example, when data were imported wih SPM8 and one
% would like to bring the hour/date info in the object to use it with
% FASST.
%_______________________________________________________________________
% Copyright (C) 2011 Cyclotron Research Centre

% Written by C. Phillips, 2011.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

if nargin==0
    % select some data file
    files = spm_select([1 Inf],'mat','Select M/EEG data files');
end

Nf = size(files,1);
for ii=1:Nf
    Di = crc_eeg_load(deblank(files(ii,:)));
    if ~isfield(Di,'info') || isempty(Di.info) || ...
            ~isfield(Di.info,'hour') || isempty(Di.info.hour) || ...
            ~isfield(Di.info,'date') || isempty(Di.info.date) || ...
            (all(Di.info.hour==[0 0 0]) && all(Di.info.date==[1 1 1]))
        [pth,fn,ext] = fileparts(Di.fname);
        lD = dir(fullfile(Di.path,[fn,'.vhdr']));
        if ~isempty(lD)
            % reset time/date
            [header,marker] = crc_eeg_readHeaderBrainP(fullfile(Di.path,lD(1).name));
            Di.info = header.info; % just in case, save the ASCII header
            try
                text_dt = marker.Date{1};
                year = str2double(text_dt(1:4));
                month = str2double(text_dt(5:6));
                day = str2double(text_dt(7:8));
                hour = str2double(text_dt(9:10));
                minutes = str2double(text_dt(11:12));
                secondes = str2double(text_dt(13:end))/1e6;

                Di.info.date = [year month day];
                Di.info.hour = [hour minutes secondes];
                save(Di)
                fprintf('\n Time and date added to file: \n\t%s\n',...
                    Di.fname);

            catch
                fprintf('\n WARNING: Could not extract the date/time of file: \n\t%s\n',...
                    Di.fname);
            end
        else
            % no header file in directory
            fprintf('\n WARNING: can''t find header file for file: \n\t%s\n',...
                Di.fname);
        end
    else
        fprintf('\n Nothing to do for file: \n\t%s\n',...
            Di.fname);
    end

end

return


