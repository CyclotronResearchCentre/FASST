function [scale]=crc_scales(D,indchan)
%get the scale(j)e of one channel depending on its units or its type. Set the
%channels units or types before hand using SPM8 to use it.
%--------------------------------------------------------------------------
%written by J.Schrouff, CRC, 05/31/2011

if nargin<1
    D = crc_eeg_load;
end
if nargin<2
    %get the scale(j)es of all channels
    indchan=1:nchannels(D);
end
scale=ones(numel(indchan),1);
def=crc_get_defaults('scales');

for j=1:numel(indchan)
    if  any(indchan(j)==emgchannels(D))  %strfind(testcond{:},'EMG')
        unemg=units(D,indchan(j));
        try
            scale(j)=eval(unemg{1});
        catch
            if strcmpi(unemg{1},'V')
                scale(j)=def.emg;
            else
                scale(j)=1;
            end
        end
    elseif any(indchan(j)==eogchannels(D)) %strfind(testcond{:},'EOG')
        unemg=units(D,indchan(j));
        try
            scale(j)=eval(unemg{1});
        catch
            if strcmpi(unemg{1},'V')
                scale(j)=def.eog;
            else
                scale(j)=1;
            end
        end
    elseif any(indchan(j)==ecgchannels(D)) %strfind(testcond{:},'ECG')
        unemg=units(D,indchan(j));
        try
            scale(j)=eval(unemg{1});
        catch
            if strcmpi(unemg{1},'V')
                scale(j)=def.ecg;
            else
                scale(j)=1;
            end
        end
    elseif any(indchan(j)==meegchannels(D)) %'EEG', 'MEGMAG', MEGPLANAR','LFP'
        chtyp=chantype(D,indchan(j));
        if strcmpi(chtyp,'EEG')
            unemg=units(D,indchan(j));
            try
                scale(j)=eval(unemg{1});
            catch
                if strcmpi(unemg{1},'V')
                    scale(j)=def.eeg;
                else
                    scale(j)=1;
                end
            end
        elseif strcmpi(chtyp,'MEGMAG')
            unemg=units(D,indchan(j));
            try
                scale(j)=eval(unemg{1});
            catch
                if strcmpi(unemg{1},'T')
                    scale(j)=def.megmag;
                else
                    scale(j)=1;
                end
            end
        elseif strcmpi(chtyp,'MEGPLANAR')
            unemg=units(D,indchan(j));
            try
                scale(j)=eval(unemg{1});
            catch
                if strcmpi(unemg{1},'T/m')
                    scale(j)=def.magplan;
                else
                    scale(j)=1;
                end
            end
         elseif strcmpi(chtyp,'LFP')
            unemg=units(D,indchan(j));
            try
                scale(j)=eval(unemg{1});
            catch
                if strcmpi(unemg{1},'V')
                    scale(j)=def.lfp;
                else
                    scale(j)=10 ;
                end
            end
        end
    else
        chtyp=chantype(D,indchan(j));
        if strcmpi(chtyp,'other')
            unemg=units(D,indchan(j));
            try
                scale(j)=eval(unemg{1});
            catch
                if strcmpi(unemg{1},'V')
                    scale(j)=def.other;
                else
                    scale(j)=1;
                end
            end
        else
            disp('Unknown type of channel, please edit')
            disp('Assuming scale is one')
            scale(j)=1;
        end
    end
end