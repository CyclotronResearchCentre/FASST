function comp_res = crc_comparison(flags)

% Comparison method for ASEEGA and FASST.
%
% ASEEGA becomes a new scorer in the FASST file called. ASEEGA scores are
% visible in FASST.
%
% To compare scores of ASEEGA and FASST, some statistics made at three levels:
%   matches on 2 conditions : wake and sleep
%   matches on 2 conditions : REM  - NREM
%   matches on 5 conditions : wake - S1 - S2 - S3&S4 - REM
% NB: MT and unscorable windows are not included in statistic
% 
% Flags is a structure composed of options:
% ----------------------------------------
%  name        default   description
%
% 'addspin'     true    includes spindles events from ASEEGA as events in FASST
% 'addab'       true    includes bursts events from ASEEGA as events in FASST
% 'pathname',   false   the path name where FASST data are
% 'addxl',      true    to save statict parameters in an excel file 
%                       the excel file is saved in the path name directory
%                       with name: 'result_comparison'.
% 'addpie',     true    to display pie charts 
% 'display',    false   to display results on command window 
%
% Statistic parameters:
% ---------------------
%
% Po = Observed agreement : number of matches over total number of windows scored
% kappa of Cohen = (P0 - Pe) / (1 - Pe) with Pe an expected probablity
%
% With FASST chosen as gold standard, and for each state:
%   TPR = True Positive Ration (sensitivity), TPR = matches in state S over
%   the total number of windows scored as S with FASST
%
%   PPV = Positive Predictive Ratio, PPV = matches in state S over the
%   total number of windows scored as S with ASEEGA


% check the flags

def_flags = struct(...
    'addspin', true,...
    'addab', true,...
    'save', true,...
    'pathname',false,...
    'addxl', true,...
    'addpie', true,...
    'display', false...
    );

if nargin
    flags = crc_check_flag(def_flags,flags);
else
    flags = def_flags;
end

% choose folder 
fn = spm_select(Inf,'mat','Select ASEEGA score file(s)', ...
    '',pwd,'.*AseegaAnalysis.A4R.mat');
[dat,Ndat] = subfct_loaddata(fn);


if ~flags.pathname
    flags.pathname = spm_select(1,'dir','Select folder with files with Fasst');
end
glob = 0;
if and(flags.addpie,size(fn,1) > 1)
    flags.addpie = 0;
    glob = 1;
end
    
for ifi = 1 : size(fn,1) 
    % import assega data 
    Do = crc_import_aseega(deblank(fn(ifi,:)), flags);
    Df = Do{1};
    % convert from 20s to 30s fasst data 
    name = Df.CRC.score(2,:);
    sw = cell2mat(Df.CRC.score(3,:));
    pc = name(sw == 30);

    if numel(pc) == 1
        % convert files
        sc2 = find(strcmp('Aseega',name));
        Do = crc_convert20to30s(Df, flags);
        Df = spm_eeg_load(Do);
        sc1 = size(Df.CRC.score,2);
    else 
        sc2 = find(strcmp('Aseega',name));
        sc1 = find((sw == 30)-strcmp('Aseega',name));
        if numel(sc1) >1
            [isc1 ok] = listdlg('PromptString','Select a file:',...
                      'SelectionMode','single',...
                      'ListString',Df.CRC.score(2,sc1));
                  sc1 = sc1(isc1);
        end    
    end
    % if more than one Aseega scorer, take the last one!
    if numel(sc2)>1
        sc2 = sc2(end);
    end
    % compare the data 
    scorer1 = cell2mat(Df.CRC.score(1,sc1));
    scorer2 = cell2mat(Df.CRC.score(1,sc2));
    % Check if the recording has been cut 
    % nb: Aseega analyses only windows with a duration of 30s!
    Nwi = nsamples(Df)/(30*fsample(Df));
    Nw = floor(nsamples(Df)/(30*fsample(Df)));
    if Nw<Nwi
        scorer1 = scorer1(1:end-1);
    end
    if numel(scorer1)~=numel(scorer2)
        %Check the fpl opl
        st_ed = Df.CRC.score{4,sc1};
        deb = st_ed(1)/30;
        deb_i = floor(deb);
        rest = deb-deb_i;
        if rest < 0.5
            startp = deb_i+1;
        else 
            startp = deb_i+2;
        end
        La = numel(scorer2);
        scorer1 = scorer1(startp:startp+La-1);
    end
    % Unscorable windows merged with MT windows
    scorer1(scorer1 == 7) = 6;
    % In cas of there are undefined windows in FASST
    scorer1 = scorer1(~isnan(scorer1));
    scorer2 = scorer2(~isnan(scorer1));
    
    flags.name = fname(Df);
    CMu = crc_compareScore(scorer1,scorer2,flags);
    if flags.addxl
        comp_res = crc_stat_aseega(CMu,flags);
        stat = [comp_res.Po comp_res.kappa comp_res.tpr comp_res.ppv comp_res.tpr_tot comp_res.ppv_tot];
        save_xl(flags.pathname, fname(Df), stat)
    end
    CM(:,:,ifi) = CMu;
end

% union of severals files
if and(flags.addxl, size(CM,3)>1)
    flags.addxl = 'union_files';
    if glob
        flags.addpie = 1;
    end
    comp_res = crc_stat_aseega(CM,flags);
    stat = [comp_res.Po comp_res.kappa comp_res.tpr comp_res.ppv comp_res.tpr_tot comp_res.ppv_tot];
    save_xl(flags.pathname, flags.addxl, stat)
end


if flags.display
    sprintf('Cohen''kappa orders: \n * < 0 : poor \n * 0.00-0.40 : fear \n * 0.41-0.60: good \n * 0.61-0.80: excellent \n * 0.81-0.10: almost perfect \n')
    fprintf(['\n * 2 conditions Sleep (S1 to S5) VS Wake (W) and no MT : Pouracy ', num2str(comp_res.Po(1)*100),' and kappa ', num2str(comp_res.kappa(1)),...
        '\n * 2 conditions NREM (S2 to S4) VS REM and no MT : Pouracy ', num2str(comp_res.Po(2)*100),' and kappa ', num2str(comp_res.kappa(2)),...
        '\n * 5 conditions (w - s1 - s2 - s3,4 - rem) : Pouracy ', num2str(comp_res.Po(3)*100),' and kappa ', num2str(comp_res.kappa(3)),...
        '\n * MT : Pouracy ', num2str(comp_res.Po(4)*100),'\n'])
    sprintf('CM matrix: lines from 1 to 6 correspond to sleep stages 0 (awake) to 5 (rem) with s3 and s4 merged.\n rows are sleep stages seen in Fasst and columns are sleep stages seen in Aseega \n')
    comp_res.CM
end
end

%% SUBFUNCTIONS
function [dat,Ndat] = subfct_loaddata(fn)

Ndat = size(fn,1);
for ii=1:Ndat
    tmp = load(deblank(fn(ii,:)));
    dat(ii) = tmp.aseega_analysis;
end
end

function save_xl(pathname, name, stat)
    icd = cd;
    cd(pathname)
    L = ls('result_comparison.xls');
    loc = fullfile(pathname,'result_comparison.xls');
    if isempty(L)
        init_excel(loc);
        s = 2;
    else 
        [num text raw]=xlsread(loc,'compare');
        s = size(text,1);
    end
    namefile = spm_str_manip(name,'s');
    xlswrite(loc,{char(namefile)},'compare',['A' num2str(s+2)])
    xlswrite(loc,stat,'compare',['B' num2str(s+2)]);
    cd(icd);
end

function init_excel(loc)
    classes = {'Name','Po S-W','Po NREM-REM','Po 5 SS','MT','k S-W','k NREM-REM','k 5SS','TPR NREM','TPR REM','PPV NREM','PPV REM','TPR w','TPR s1','TPR S2','TPR S3,4','TPR REM','PPV w','PPV S1','PPV S2','PPV S3,4','PPV REM'};        
    xlswrite(loc,classes,'compare','A2') 
end