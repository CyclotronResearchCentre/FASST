function data = COF_data(pathtodata)

% Modifications:
% - 20171222    MV merged and reorganized all the existing COF_DATA.m into one complete data structure

% Defining required resources
data = struct('DIR',pathtodata, ... %dir of the study
    'num',[],'weight',[],'height',[],'id',[],'dir',[],'THKID',[],'THKDose',[],...
    'PTname',[],'MRI',[],'PTfile',[]','sessname',[],'sspname',[],'EEGfile',[],...
    'Sname',[],'channel',[],'closestE',[],'EEGsession',[],'LongMRI',[],'T1w',[],'PDw',[],'MTw',[],...
    'B1',[],'B0magn',[],'B0phase',[],'DCMfile',[],'FPL_OPL',[],'Sleep_cycles',[],...
    'wkTime',[],'SleepTime',[],'dlmo',[],'twoback',[],'threeback',[],'sart',[],...
    'pvt',[],'kss_alone',[],'kss_tasks',[],'kss_tms',[],'vas_alone',[],'vas_tasks',[],...
    'vas_tms',[],'pupil_cal',[],'pupil_ER',[]);
num = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% This is a template. Please use it to create new subjects %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% 0 - COF000
% num                     = num + 1;
% data(num).num           = num;
% 
% % Informations
% data(num).weight        = ;     % Kg
% data(num).height        = ;    % Cm
% data(num).id            = '';
% data(num).dir           = fullfile(pathtodata,data(num).id);
% % PET Scans
% data(num).THKID         = '';
% data(num).THKDose       = ;
% % TMS Pretests
% data(num).PTname        = ['P1' ; 'P2'];
% data(num).MRI           = {',1'};
% data(num).PTfile        = {{'.mat'} ...   %P1
%     {'.mat'} ...                        %P2
%     };
% % TMS Sessions
% data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
% data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
% data(num).EEGfile     = {{'.mat'} ... %S1
%     {'.mat'} ... %S2
%     {'.mat'} ... %S3
%     {'.mat'} ... %S4
%     {'.mat'} ... %S5
%     };
% data(num).Sname         = {'LSMA'};
% data(num).channel       = '';
% data(num).closestE      = ;
% data(num).EEGsession    = [];
% % MRI MPM
% data(num).LongMRI       = '';
% data(num).T1w           = '';
% data(num).PDw           = '';
% data(num).MTw           = '';
% data(num).B1            = '';
% data(num).B0magn        = '';
% data(num).B0phase       = '';
% % DCM
% data(num).DCMfile       = {'.mat'};
% % BL
% data(num).FPL_OPL = [];
% data(num).Sleep_cycles  = {};
% % Circadian Phase
% data(num).wkTime        = '';
% data(num).SleepTime     = '';
% data(num).dlmo          = '';
% % Cognitive tasks
% data(num).twoback       = [];
% data(num).threeback     = [];
% data(num).sart          = [];
% data(num).pvt           = [];
% data(num).kss_alone     = [];
% data(num).kss_tasks     = [];
% data(num).kss_tms       = [];
% data(num).vas_alone     = [];
% data(num).vas_tasks     = [];
% data(num).vas_tms       = [];
% % Pupillometry
% data(num).pupil_cal     = [];
% data(num).pupil_ER      = []; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1 - COF002
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 60;     % Kg
data(num).height        = 160;    % Cm
data(num).id            = 'COF002';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = 'E0709';
data(num).THKDose       = 190;
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s845-0002-00001-000224-01.nii,1'};
data(num).PTfile        = {{'spmeeg_COF002_PT1_20160609_150353_4.mat'} ...  %P1
    {'spmeeg_COF002_PT2_20160609_150353_10.mat'} ...                        %P2
    };
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5' ; 'S6'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5' ; 'S6'];
data(num).EEGfile       = {{'spmeeg_COF002_S1_EEGTMS_20160615_103442_4.mat'} ...    %S1
    {'spmeeg_COF002_S2_EEGTMS_20160615_143347_2.mat'} ...                           %S2
    {'spmeeg_COF002_S3_EEGTMS_20160615_174210_3.mat'} ...                           %S3
    {'spmeeg_COF002_S4_TMSEEG_0160615_201630_2.mat'} ...                            %S4
    {'spmeeg_COF002_S5_TMSEEG_20160615_223919_2.mat'} ...                           %S5
    {'spmeeg_COF002_S6_TMSEEG_20160616_024056_2.mat'} ...                           %S6
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 28;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1337';
data(num).T1w           = '8';
data(num).PDw           = '10';
data(num).MTw           = '12';
data(num).B1            = '14';
data(num).B0magn        = '15';
data(num).B0phase       = '16';
% DCM
data(num).DCMfile       = {'merged_fmbLAPefdfMspmeeg_COF002_S1_EEGTMS_20160615_103442_4.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 8;
data(num).SleepTime     = 1;
data(num).dlmo          = 19.89;
% Cognitive tasks
data(num).twoback       = [4 8];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [5 7];
data(num).kss_tasks     = [2];
data(num).kss_tms       = [];
data(num).vas_alone     = [5 7];
data(num).vas_tasks 	= [2];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [2];
data(num).pupil_ER      = [2 7]; 

%% 2 - COF004
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 92;     % Kg
data(num).height        = 181;    % Cm
data(num).id            = 'COF004';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = 'E0573';
data(num).THKDose       = 190;
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s939-0005-00001-000224-01.nii,1'};
data(num).PTfile        = {{'spmeeg_COF004_P1_20160705_090823_2.mat'} ...   %P1
    {'spmeeg_COF004_P2_20160705_090823_3_2.mat'} ...                        %P2
    };
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF004_S1_TMSEEG_20160713_071709_4.mat'} ... %S1
    {'spmeeg_COF004_S2_TMSEEG_20160713_115942_2.mat'} ... %S2
    {'spmeeg_COF004_S3_TMSEEG_20160713_154558_2.mat'} ... %S3
    {'spmeeg_COF004_S4_TMSEEG_20160713_192459_3.mat'} ... %S4
    {'spmeeg_COF004_S5_TMSEEG_20160713_222122_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 9;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '993';
data(num).T1w           = '2';
data(num).PDw           = '4';
data(num).MTw           = '6';
data(num).B1            = '8';
data(num).B0magn        = '9';
data(num).B0phase       = '10';
% DCM
data(num).DCMfile       = {'merged_fmbLAPefdfMspmeeg_COF004_S1_TMSEEG_20160713_071709_4.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 5.5;
data(num).SleepTime     = 22.5;
data(num).dlmo          = 22.34;
% Cognitive tasks
data(num).twoback       = [7 8 9];      % D' below 0 - to check
data(num).threeback     = [6 7 8 9];    % D' below 0 - to check
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks 	= [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [2];
data(num).pupil_ER      = [2]; 

%% 3 - COF008
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 74;     % Kg
data(num).height        = 169;    % Cm
data(num).id            = 'COF008';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = 'E0586';
data(num).THKDose       = 173;
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s914-0002-00001-000224-01.nii,1'};
data(num).PTfile     = {{'spmeeg_COF008_P1_20160624_103503_3.mat'} ... %P1
    {'spmeeg_COF008_P2_0160624_103503_5.mat'} ...%P2
    };
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF008_S1_TMSEEG_20160825_082122_3.mat'} ... %S1
    {'spmeeg_COF008_S2_TMSEEG_20160825_114617_2.mat'} ... %S2
    {'spmeeg_COF008_S3_TMSEEG_20160825_174518_2.mat'} ... %S3
    {'spmeeg_COF008_S4_TMSEEG_20160825_211614_2.mat'} ... %S4
    {'spmeeg_COF008_S5_TMSEEG_20160825_233422_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 19;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '916';
data(num).T1w           = '2';
data(num).PDw           = '4';
data(num).MTw           = '6';
data(num).B1            = '8';
data(num).B0magn        = '9';
data(num).B0phase       = '10';
% DCM
data(num).DCMfile       = {'merged_fmbLAPefdfMspmeeg_COF008_S1_TMSEEG_20160825_082122_3.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 7;
data(num).SleepTime     = 23;
data(num).dlmo          = 20.94;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = []; 

%% 4 - COF009
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 58;     % Kg
data(num).height        = 144;    % Cm
data(num).id            = 'COF009';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = '';
data(num).THKDose       = '';
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s1006-0005-00001-000224-01.nii,1'};
data(num).PTfile     = {{'spmeeg_COF009_P1_20160920_094402_5.mat'} ... %P1
{'spmeeg_COF009_P2_20160920_094402_8.mat'} ...%P2
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF009_S1_TMSEEG_20161122_081554_3.mat'} ... %S1
    {'spmeeg_COF009_S2_TMSEEG_20161122_150037_2.mat'} ... %S2
    {'spmeeg_COF009_S3_TMSEEG_20161122_161826_2.mat'} ... %S3
    {'spmeeg_COF009_S4_TMSEEG_20161122_200114_2.mat'} ... %S4
    {'spmeeg_COF009_S5_TMSEEG_20161122_222154_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 10;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1466';
data(num).T1w           = '2';
data(num).PDw           = '4';
data(num).MTw           = '6';
data(num).B1            = '8';
data(num).B0magn        = '9';
data(num).B0phase       = '10';
% DCM
data(num).DCMfile       = {'merged_fmbLAPefdfMspmeeg_COF009_S1_TMSEEG_20161122_081554_3.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 7;
data(num).SleepTime     = 22.5;
data(num).dlmo          = 18.52;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = []; 


%% 5 - COF012
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 59;     % Kg
data(num).height        = 165;    % Cm
data(num).id            = 'COF012';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = 'E0591';
data(num).THKDose       = '';
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s940-0005-00001-000224-01.nii,1'};
data(num).PTfile     = {{'spmeeg_COF012_20160705_124512_4.mat'} ... %P1
{'spmeeg_COF012_P2_20160705_124512_8.mat'} ...%P2
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF012_S1_TMSEEG_20160713_091521_3.mat'} ... %S1
    {'spmeeg_COF012_S2_TMSEEG_20160713_144852_2.mat'} ... %S2
    {'spmeeg_COF012_S3_TMSEEG_20160713_184024_2.mat'} ... %S3
    {'spmeeg_COF012_S4_TMSEEG_20160713_210803_3.mat'} ... %S4
    {'spmeeg_COF012_S5_TMSEEG_20160714_010232_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 17;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '940';
data(num).T1w           = '8';
data(num).PDw           = '10';
data(num).MTw           = '12';
data(num).B1            = '14';
data(num).B0magn        = '15';
data(num).B0phase       = '16';
% DCM
data(num).DCMfile       = {'merged_fmbLAPefdfMspmeeg_COF012_S1_TMSEEG_20160713_091521_3.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 6.5;
data(num).SleepTime     = 22.5;
data(num).dlmo          = 20.06;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = []; 

%% 6 - COF013
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 64;     % Kg
data(num).height        = 170;    % Cm
data(num).id            = 'COF013';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = 'E0716';
data(num).THKDose       = 185;
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s994-0014-00001-000224-01.nii,1'};
data(num).PTfile     = {{'spmeeg_COF013_P1_20160822_140742_6.mat'} ... %P1
{'spmeeg_COF013_P2_20160822_140742_9.mat'} ...%P2
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF013_S1_TMSEEG_20160825_101432_4.mat'} ... %S1
    {'spmeeg_COF013_S2_TMSEEG_20160825_162111_4.mat'} ... %S2
    {'spmeeg_COF013_S3_TMSEEG_20160825_200654_2.mat'} ... %S3
    {'spmeeg_COF013_S4_TMSEEG_20160825_222953_2.mat'} ... %S4
    {'spmeeg_COF013_S5_TMSEEG_20160826_021744_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 10;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1009';
data(num).T1w           = '8';
data(num).PDw           = '10';
data(num).MTw           = '12';
data(num).B1            = '14';
data(num).B0magn        = '15';
data(num).B0phase       = '16';
% DCM
data(num).DCMfile       = {'merged_fmbLAPefdfMspmeeg_COF013_S1_TMSEEG_20160825_101432_4.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 8;
data(num).SleepTime     = 23.25;
data(num).dlmo          = 20.83;
% Cognitive tasks
data(num).twoback       = [9];      % D' below 0 - To check
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [3];
data(num).pupil_ER      = [3]; 

%% 7 - COF014
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 81;     % Kg
data(num).height        = 181;    % Cm
data(num).id            = 'COF014';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = 'E0710';
data(num).THKDose       = 190;
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s995-0005-00001-000224-01.nii,1'};
data(num).PTfile     = {{'spmeeg_COF014_P1_20160817_122626_6.mat'} ... %P1
{'spmeeg_COF014_P2_20160817_122626_8.mat'} ...%P2
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF014_S1_TMSEEG_20160907_100508_4.mat'} ... %S1
    {'spmeeg_COF014_S2_TMSEEG_20160907_163015_2.mat'} ... %S2
    {'spmeeg_COF014_S4_TMSEEG_20160907_223114_2.mat'} ... %S4
    {'spmeeg_COF014_S5_TMSEEG_20160908_023028_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 28;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '997';
data(num).T1w           = '2';
data(num).PDw           = '4';
data(num).MTw           = '6';
data(num).B1            = '8';
data(num).B0magn        = '9';
data(num).B0phase       = '10';
% DCM
data(num).DCMfile       = {'merged_fmbLAPefdfMspmeeg_COF014_S1_TMSEEG_20160907_100508_4.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 8;
data(num).SleepTime     = 0;
data(num).dlmo          = 22.18;
% Cognitive tasks
data(num).twoback       = [9];      % D' below 0 - To check
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [3];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [3];
% Pupillometry
data(num).pupil_cal     = [8 9];
data(num).pupil_ER      = [8 9]; 

%% 8 - COF015
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 71;     % Kg
data(num).height        = 168;    % Cm
data(num).id            = 'COF015';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = 'E0722';
data(num).THKDose       = 199;
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s1007-0005-00001-000224-01.nii,1'};
data(num).PTfile     = {{'spmeeg_COF015_P1_20160920_120500_8.mat'} ... %P1
{'spmeeg_COF015_P2_20160920_120500_11.mat'} ...%P2
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF015_S1_TMSEEG_20170307_080820_6.mat'} ... %S1
    {'spmeeg_COF015_S2_20170307_111547_3.mat'} ... %S2
    {'spmeeg_COF015_S3_TMSEEG_20170307_185755_2.mat'} ... %S3
    {'spmeeg_COF015_S4_TMSEEG_20170307_212048_3.mat'} ... %S4
    {'spmeeg_COF015_S5_20170308_012700_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 28;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1007';
data(num).T1w           = '8';
data(num).PDw           = '10';
data(num).MTw           = '12';
data(num).B1            = '14';
data(num).B0magn        = '15';
data(num).B0phase       = '16';
% DCM
data(num).DCMfile       = {'merged_fmbLAPefdfMspmeeg_COF015_S1_TMSEEG_20170307_080820_6.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 7;
data(num).SleepTime     = 22.5;
data(num).dlmo          = 20.83;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [6];
data(num).pupil_ER      = [1 2 6 7];        % 6 is missing. Others are of poor quality - To check 

%% 9 - COF016
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 64;     % Kg
data(num).height        = 161;    % Cm
data(num).id            = 'COF016';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = 'E0577';
data(num).THKDose       = 191;
% TMS Pretests
data(num).PTname      = ['P1' ; 'P2' ; 'P3'];
data(num).MRI           = {'s989-0005-00001-000224-01.nii,1'};
data(num).PTfile     = {{'spmeeg_COF016_P1_20160803_124831_6.mat'} ... %P1
{'spmeeg_COF016_P2_20160803_124831_9.mat'} ...%P2
{'spmeeg_COF016_P3_20160803_124831_14.mat'} ...%P3
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF016_S1_TMSEEG_20160811_101158_3.mat'} ... %S1
    {'spmeeg_COF016_S2_TMSEEG_20160811_163013_2.mat'} ... %S2
    {'spmeeg_COF016_S3_TMSEEG_20160811_200454_2.mat'} ... %S3
    {'spmeeg_COF016_S4_TMSEEG_20160811_223726_2.mat'} ... %S4
    {'spmeeg_COF016_S5_TMSEEG_20160812_023633_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 10;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '989';
data(num).T1w           = '8';
data(num).PDw           = '10';
data(num).MTw           = '12';
data(num).B1            = '14';
data(num).B0magn        = '15';
data(num).B0phase       = '16';
% DCM
data(num).DCMfile       = {'merged_fmbLAPefdfMspmeeg_COF016_S1_TMSEEG_20160811_101158_3.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 8;
data(num).SleepTime     = 0;
data(num).dlmo          = 19.79;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];       % To check
data(num).kss_tasks     = [];       % To check
data(num).kss_tms       = [];       % To check
data(num).vas_alone     = [];       % To check
data(num).vas_tasks     = [];       % To check
data(num).vas_tms       = [];       % To check
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = []; 

%% 10 - COF017
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 83;     % Kg
data(num).height        = 182;    % Cm
data(num).id            = 'COF017';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = 'E0578';
data(num).THKDose       = 194;
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s988-0006-00001-000224-01.nii,1'};
data(num).PTfile     = {{'spmeeg_COF017_P1_20160803_095649_7.mat'} ... %P1
{'spmeeg_COF017_P2_20160803_095649_10.mat'} ...%P2
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF017_S1_TMSEEG_20160811_081922_3.mat'} ... %S1
    {'spmeeg_COF017_S2_TMSEEG_20160811_110607_2.mat'} ... %S2
    {'spmeeg_COF017_S3_TMSEEG_20160811_174536_2.mat'} ... %S3
    {'spmeeg_COF017_S4_TMSEEG_20160811_210107_2.mat'} ... %S4
    {'spmeeg_COF017_S5_TMSEEG_20160811_233922_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 28;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '990';
data(num).T1w           = '2';
data(num).PDw           = '4';
data(num).MTw           = '6';
data(num).B1            = '8';
data(num).B0magn        = '9';
data(num).B0phase       = '10';
% DCM
data(num).DCMfile       = {'merged_fmbLAPefdfMspmeeg_COF017_S1_TMSEEG_20160811_081922_3.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 7;
data(num).SleepTime     = 22.5;
data(num).dlmo          = 21.41;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = []; 

%% 11 - COF020
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 66;     % Kg
data(num).height        = 166;    % Cm
data(num).id            = 'COF020';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = 'E0588';
data(num).THKDose       = 196;
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s999-0005-00001-000224-01.nii,1'};
data(num).PTfile     = {{'spmeeg_COF020_P1_20160830_103128_8.mat'} ... %P1
{'spmeeg_COF020_P2_20160830_103128_12.mat'} ...%P2
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF020_S1_TMSEEG_20160907_004057_3.mat'} ... %S1
    {'spmeeg_COF020_S2_TMSEEG_20160907_111422_2.mat'} ... %S2
    {'spmeeg_COF020_S3_TMSEEG_20160907_172546_1.mat'} ... %S3
    {'spmeeg_COF020_S4_TMSEEG_20160907_212919_2.mat'} ... %S4
    {'spmeeg_COF020_S5_TMSEEG_20160907_232617_6.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 19;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '999';
data(num).T1w           = '8';
data(num).PDw           = '10';
data(num).MTw           = '12';
data(num).B1            = '14';
data(num).B0magn        = '15';
data(num).B0phase       = '16';
% DCM
data(num).DCMfile       = {'merged_fmbLAPefdfMspmeeg_COF020_S1_TMSEEG_20160907_004057_3.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 7;
data(num).SleepTime     = 23.5;
data(num).dlmo          = 20.17;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [7];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [7];
data(num).vas_tasks     = [6];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = []; 

%% 12 - COF023
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 76;     % Kg
data(num).height        = 174;    % Cm
data(num).id            = 'COF023';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = 'E0715';
data(num).THKDose       = 202;
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s1040-0005-00001-000224-01.nii,1'};
data(num).PTfile     = {{'spmeeg_COF023_P1_20161026_104129_4.mat'} ... %P1
{'spmeeg_COF023_P2_20161026_104129_6.mat'} ...%P2
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['1' ; '2' ; '3' ; '4' ; '5'];
data(num).EEGfile     = {{'spmeeg_COF023_S1_TMSEEG_20161116_082955_5.mat'} ... %S1
    {'spmeeg_COF023_S2_TMSEEG_20161116_111907_2.mat'} ... %S2
    {'spmeeg_COF023_S3_TMSEEG_20161116_180514_4.mat'} ... %S3
    {'spmeeg_COF023_S4_TMSEEG_20161116_212304_2.mat'} ... %S4
    {'spmeeg_COF023_S5_TMSEEG_20161117_023202_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 19;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1040';
data(num).T1w           = '8';
data(num).PDw           = '10';
data(num).MTw           = '13';
data(num).B1            = '15';
data(num).B0magn        = '16';
data(num).B0phase       = '17';
% DCM
data(num).DCMfile       = {'merged_fmbLAPefdfMspmeeg_COF023_S1_TMSEEG_20161116_082955_5.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 8;
data(num).SleepTime     = 22.5;
data(num).dlmo          = 20.44;        % Not that clear
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [1 2 3];
data(num).pupil_ER      = [1 2 3]; 

%% 13 - COF024
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 100;     % Kg
data(num).height        = 188;    % Cm
data(num).id            = 'COF024';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = 'E0707';
data(num).THKDose       = 195;
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s1020-0005-00001-000224-01.nii,1'};
data(num).PTfile     = {{'spmeeg_COF024_P1_20160930_132228_7.mat'} ... %P1
{'spmeeg_COF024_P2_20160930_132228_12.mat'} ...%P2
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF024_S1_TMSEEG_20161007_075347_3.mat'} ... %S1
    {'spmeeg_COF024_S2_TMSEEG_20161007_104250_4.mat'} ... %S2
    {'spmeeg_COF024_S3_TMSEEG_20161007_171812_2.mat'} ... %S3
    {'spmeeg_COF024_S4_TMSEEG_20161007_204209_2.mat'} ... %S4
    {'spmeeg_COF024_S5_TMSEEG_20161007_225748_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 18;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1020';
data(num).T1w           = '8';
data(num).PDw           = '10';
data(num).MTw           = '12';
data(num).B1            = '14';
data(num).B0magn        = '15';
data(num).B0phase       = '16';
% DCM
data(num).DCMfile       = {'merged_fmbLAPefdfMspmeeg_COF024_S1_TMSEEG_20161007_075347_3.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 6.5;
data(num).SleepTime     = 23.5;
data(num).dlmo          = 21.23;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [4];
data(num).pupil_ER      = [4]; 

%% 14 - COF025
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 54;     % Kg
data(num).height        = 163;    % Cm
data(num).id            = 'COF025';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = 'E0735';
data(num).THKDose       = '';
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s1026-0005-00001-000224-01.nii,1'};
data(num).PTfile     = {{'spmeeg_COF025_P1_20161006_101426_9.mat'} ... %P1
{'spmeeg_COF025_P2_20161006_101426_13.mat'} ...%P2
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF025_S1_TMSEEG_20161021_072544_6.mat'} ... %S1
    {'spmeeg_COF025_S2_TMSEEEG_20161021_144027_2.mat'} ... %S2
    {'spmeeg_COF025_S3_TMSEEG_20161021_153749_2.mat'} ... %S3
    {'spmeeg_COF025_S4_TMSEEG_20161021_202327_2.mat'} ... %S4
    {'spmeeg_COF025_S5_TMSEEG_20161022_002540_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 9;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1026';
data(num).T1w           = '8';
data(num).PDw           = '10';
data(num).MTw           = '12';
data(num).B1            = '14';
data(num).B0magn        = '15';
data(num).B0phase       = '16';
% DCM
data(num).DCMfile       = {'.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 6;
data(num).SleepTime     = 22.75;
data(num).dlmo          = 20.6;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = []; 

%% 15 - COF026
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 85;     % Kg
data(num).height        = 170;    % Cm
data(num).id            = 'COF026';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = '';
data(num).THKDose       = '';
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s1035-0005-00001-000224-01.nii,1'};
data(num).PTfile     = {{'spmeeg_COF026_P1_20161020_094201_14.mat'} ... %P1
{'spmeeg_COF026_P2_20161020_094201_16.mat'} ...%P2
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF026_S1_TMSEEG_20170316_090832_6.mat'} ... %S1
    {'spmeeg_COF026_S2_TMSEEG_20170316_153236_3.mat'} ... %S2
    {'spmeeg_COF026_S3_TMSEEG_20170316_190532_2.mat'} ... %S3
    {'spmeeg_COF026_S4_TMSEEG_20170316_214132_2.mat'} ... %S4
    {'spmeeg_COF206_S5_TMSEEG_170317_015941_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 9;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1046';
data(num).T1w           = '2';
data(num).PDw           = '4';
data(num).MTw           = '6';
data(num).B1            = '8';
data(num).B0magn        = '9';
data(num).B0phase       = '10';
% DCM
data(num).DCMfile       = {'merged_fmbLAPefdfMspmeeg_COF026_S1_TMSEEG_20170316_090832_6.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 7.25;
data(num).SleepTime     = 23;
data(num).dlmo          = 19.92;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [8];
data(num).pupil_ER      = [8]; 

%% 16 - COF027
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 83;     % Kg
data(num).height        = 179;    % Cm
data(num).id            = 'COF027';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = '';
data(num).THKDose       = '';
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s160930090910DST131221107524366021-0005-00001-000224-01.nii,1'};
data(num).PTfile     = {{'spmeeg_COF027_P1_20160930_094435_6.mat'} ... %P1
{'spmeeg_COF027_P2_20160930_094435_14.mat'} ...%P2
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF027_S1_TMSEEG_20161007_094328_3.mat'} ... %S1
    {'spmeeg_COF027_S2_TMSEEG_20161007_155810_5.mat'} ... %S2
    {'spmeeg_COF027_S3_TMSEEG_20161007_193622_3.mat'} ... %S3
    {'spmeeg_COF027_S4_TMSEEG_20161007_214919_2.mat'} ... %S4
    {'spmeeg_COF027_S5_TMSEEG_20161008_014227_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 18;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1422';
data(num).T1w           = '2';
data(num).PDw           = '4';
data(num).MTw           = '6';
data(num).B1            = '8';
data(num).B0magn        = '9';
data(num).B0phase       = '10';
% DCM
data(num).DCMfile       = {'merged_fmbLAPefdfMspmeeg_COF027_S1_TMSEEG_20161007_094328_3.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 7.5;
data(num).SleepTime     = 23;
data(num).dlmo          = 20.13;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = []; 

%% 17 - COF028
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 53;     % Kg
data(num).height        = 158;    % Cm
data(num).id            = 'COF028';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = 'E0711';
data(num).THKDose       = 191;
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s1049-0005-00001-000224-01.nii,1'};
data(num).PTfile     = {{'spmeeg_COF028_P1_20161111_103234_6.mat'} ... %P1
{'spmeeg_COF028_P2_20161111_103234_9.mat'} ...%P2
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname     = ['1' ; '2' ; '3' ; '4' ; '5'];
data(num).EEGfile     = {{'spmeeg_COF028_S1_TMSEEG_20161126_092529_3.mat'} ... %S1
    {'spmeeg_COF028_S2_TMSEEG_20161126_155038_2.mat'} ... %S2
    {'spmeeg_COF028_S3_TMSEEG_20161126_193801_2.mat'} ... %S3
    {'spmeeg_COF028_S4_TMSEEG_20161126_220622_2.mat'} ... %S4
    {'spmeeg_COF028_S5_TMSEEG_20161127_014102_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 28;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1049';
data(num).T1w           = '8';
data(num).PDw           = '10';
data(num).MTw           = '12';
data(num).B1            = '14';
data(num).B0magn        = '15';
data(num).B0phase       = '16';
% DCM
data(num).DCMfile       = {'merged_fmbLAPefdfMspmeeg_COF028_S1_TMSEEG_20161126_092529_3.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% % Circadian Phase
data(num).wkTime        = 7.5;
data(num).SleepTime     = 22.5;
data(num).dlmo          = 20.39;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = []; 

%% 18 - COF034
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 72;     % Kg
data(num).height        = 171;    % Cm
data(num).id            = 'COF034';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = 'E0712';
data(num).THKDose       = 189;
% TMS Pretests
data(num).PTname      = ['P1' ; 'P2' ; 'P3'];
data(num).MRI           = {'s1213-0005-00001-000224-01.nii,1'};
data(num).PTfile     = {{'spmeeg_COF034_P1_20170309_091607_10.mat'} ... %P1
{'spmeeg_COF034_P2_20170309_091607_14.mat'} ...%P2
{'spmeeg_COF034_P3_20170309_091607_15.mat'} ...%P3
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF034_S1_TMSEEG_20170315_165742_8.mat'} ... %S1
    {'spmeeg_COF034_S2_TMSEEG_20170316_132050_2.mat'} ... %S2
    {'spmeeg_COF034_S3_TMSEEG_20170316_164556_3.mat'} ... %S3
    {'spmeeg_COF034_S4_TMSEEG_20170316_202706_2.mat'} ... %S4
    {'spmeeg_COF034_S5_TMSEEG_20170317_003357_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 18;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1437';
data(num).T1w           = '2';
data(num).PDw           = '4';
data(num).MTw           = '6';
data(num).B1            = '8';
data(num).B0magn        = '9';
data(num).B0phase       = '10';
% DCM
data(num).DCMfile       = {'merged_fmbLAPefdfMspmeeg_COF034_S1_TMSEEG_20170315_165742_8.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 6.25;
data(num).SleepTime     = 22.5;
data(num).dlmo          = 20.56;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = [1 3]; 

%% 19 - COF035
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 64;     % Kg
data(num).height        = 160;    % Cm
data(num).id            = 'COF035';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = 'E0733';
data(num).THKDose       = '';
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s1274-0005-00001-000224-01.nii,1'};
data(num).PTfile     = {{'spmeeg_COF035_P1_20170425_165227_8.mat'} ... %P1
{'spmeeg_COF035_P2_20170425_165227_10.mat'} ...%P2
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF035_S1_TMSEEG_20170428_074759_7.mat'} ... %S1
    {'spmeeg_COF035_S2_TMSEEG_20170428_133440_2.mat'} ... %S2
    {'spmeeg_COF035_S3_TMSEEG_20170428_184956_3.mat'} ... %S3
    {'spmeeg_COF035_S4_TMSEEG_20170428_203227_2.mat'} ... %S4
    {'spmeeg_COF035_S5_TMSEEG_20170428_232602_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 9;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1480';
data(num).T1w           = '5';
data(num).PDw           = '7';
data(num).MTw           = '9';
data(num).B1            = '11';
data(num).B0magn        = '12';
data(num).B0phase       = '13';
% DCM
data(num).DCMfile       = {'merged_fmbLAPefdfMspmeeg_COF035_S1_TMSEEG_20170428_074759_7.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 6.5;
data(num).SleepTime     = 23.5;
data(num).dlmo          = 22.58;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [8];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [1 2];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = [1]; 

%% 20 - COF036
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 53;     % Kg
data(num).height        = 160;    % Cm
data(num).id            = 'COF036';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = 'E0708';
data(num).THKDose       = 186;
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s1044-0005-00001-000224-01.nii,1'};
data(num).PTfile     = {{'spmeeg_COF036_P1_20161108_093501_3.mat'} ... %P1
{'spmeeg_COF036_P2_20161108_093501_6.mat'} ...%P2
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF036_S1_EEGTMS_20170307_101536_3.mat'} ... %S1
    {'spmeeg_COF036_S2_20170307_172807_2.mat'} ... %S2
    {'spmeeg_COF036_S3_TMSEEG_20170307_200004_2.mat'} ... %S3
    {'spmeeg_COF036_S4_20170307_222611_2.mat'} ... %S4
    {'spmeeg_COF036_S5_0170308_015944_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 17;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1465';
data(num).T1w           = '2';
data(num).PDw           = '4';
data(num).MTw           = '6';
data(num).B1            = '8';
data(num).B0magn        = '9';
data(num).B0phase       = '10';
% DCM
data(num).DCMfile       = {'merged_fmbLAPefdfMspmeeg_COF036_S1_EEGTMS_20170307_101536_3.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 8;
data(num).SleepTime     = 23.5;
data(num).dlmo          = 20.15;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [6];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [3 4 8];
data(num).pupil_ER      = [2 3 4 5 7 8 9];

%% 21 - COF039
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 58;     % Kg
data(num).height        = 170;    % Cm
data(num).id            = 'COF039';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = 'E0574';
data(num).THKDose       = 197;
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2' ; 'P3'];
data(num).MRI           = {'s1060-0005-00001-000224-01.nii,1'};
data(num).PTfile        = {{'spmeeg_COF039_P1_20161118_103439_2.mat'} ... %P1
 {'spmeeg_COF039_P2_P220161118_103439_5.mat'} ...%P2
 {'spmeeg_COF039_P3_20161118_103439_10.mat'} ...%P3
 };
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname     = ['1' ; '2' ; '3' ; '4' ; '5'];
data(num).EEGfile     = {{'spmeeg_COF039_S1_TMSEEG_20161126_075844_4.mat'} ... %S1
    {'spmeeg_COF039_S2_TMSEEG_20161126_102721_2.mat'} ... %S2
    {'spmeeg_COF039_S3_TMSEEG_20161126_165945_2.mat'} ... %S3
    {'spmeeg_COF039_S4_TMSEEG_0161126_203328_2.mat'} ... %S4
    {'spmeeg_COF039_S5_TMSEEG_20161126_225339_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 17;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1133';
data(num).T1w           = '2';
data(num).PDw           = '4';
data(num).MTw           = '6';
data(num).B1            = '8';
data(num).B0magn        = '9';
data(num).B0phase       = '10';
% DCM
data(num).DCMfile       = {'merged_fmbLAPefdfMspmeeg_COF039_S1_TMSEEG_20161126_075844_4.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 6.5;
data(num).SleepTime     = 22.5;
data(num).dlmo          = 21.36;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [5];
data(num).pupil_ER      = [5 6]; 

%% 22 - COF040
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 62.5;     % Kg
data(num).height        = 166;    % Cm
data(num).id            = 'COF040';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = 'E0752';
data(num).THKDose       = '';
% TMS Pretests
data(num).PTname        = ['P1'];
data(num).MRI           = {'s1053-0005-00001-000224-01,1'};
data(num).PTfile        = {{'spmeeg_COF040_P1_20170905_150821_3.mat'} ... %P1
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF040_S1_TMSEEG_20170914_071847_3.mat'} ... %S1
    {'spmeeg_COF040_S2_TMSEEG_20170914_100419_2.mat'} ... %S2
    {'spmeeg_COF040_S3_TMSEEG_20170914_173024_2.mat'} ... %S3
    {'spmeeg_COF040_S4_TMSEEG_20170914_200314_2.mat'} ... %S4
    {'spmeeg_COF040_S5_EEGTMS_20170914_222323_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 28;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1472';
data(num).T1w           = '2';
data(num).PDw           = '4';
data(num).MTw           = '6';
data(num).B1            = '8';
data(num).B0magn        = '9';
data(num).B0phase       = '10';
% DCM
data(num).DCMfile       = {'.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 6;
data(num).SleepTime     = 22;
data(num).dlmo          = 19.82;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];       % 1, 2, and 3 are "KSS_VAS_Alone" instead of "KSS_VAS_TMS"
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];       % 1, 2, and 3 are "KSS_VAS_Alone" instead of "KSS_VAS_TMS"
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = []; 


%% 23 - COF043
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 63;     % Kg
data(num).height        = 151;    % Cm
data(num).id            = 'COF043';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = 'E0706';
data(num).THKDose       = 186;
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s1214-0005-00001-000224-01.nii,1'};
data(num).PTfile     = {{'spmeeg_COF043_P1_20170309_125711_3.mat'} ... %P1
{'spmeeg_COF043_P2_20170309_125711_6.mat'} ...%P2
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF043_S1_TMSEEG_20170407_065740_5.mat'} ... %S1
    {'spmeeg_COF043_S2_TMSEEG_20170407_101814_2.mat'} ... %S2
    {'spmeeg_COF043_S3_TMSEEG_20170407_161938_2.mat'} ... %S3
    {'spmeeg_COF043_S4_TMSEEG-20170407_201719_2.mat'} ... %S4
    {'spmeeg_COF043_S5_TMSEEG_20170407_223027_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 19;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1214';
data(num).T1w           = '8';
data(num).PDw           = '10';
data(num).MTw           = '12';
data(num).B1            = '14';
data(num).B0magn        = '15';
data(num).B0phase       = '16';
% DCM
data(num).DCMfile       = {'merged_fmbLAPefdfMspmeeg_COF043_S1_TMSEEG_20170407_065740_5.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 6;
data(num).SleepTime     = 23.5;
data(num).dlmo          = 20.38;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];       
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];       
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = []; 


%% 24 - COF044
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 83.5;     % Kg
data(num).height        = 185;    % Cm
data(num).id            = 'COF044';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = '';
data(num).THKDose       = '';
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2' ; 'P3'];
data(num).MRI           = {'s1509-0002-00001-000224-01,1'};
data(num).PTfile     = {{'spmeeg_COF044_P1_20171018_101005_8.mat'} ... %P1
{'spmeeg_COF044_P2_20171018_101005_13.mat'} ... %P2
{'spmeeg_COF044_P3_20171018_101005_16.mat'}... %P3
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF044_S1_TMSEEG_20171101_064627_6.mat'} ... %S1
    {'spmeeg_COF044_S2_TMSEEG_20171101_145408_2.mat'} ... %S2
    {'spmeeg_COF044_S3_TMSEEG_20171101_161518_2.mat'} ... %S3
    {'spmeeg_COF044_S4_TMSEEG_20171101_185553_2.mat'} ... %S4
    {'spmeeg_COF044_S5_TMSEEG_20171102_000632_2.mat'} ... %S5
       };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 28;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1509';
data(num).T1w           = '5';
data(num).PDw           = '7';
data(num).MTw           = '9';
data(num).B1            = '11';
data(num).B0magn        = '12';
data(num).B0phase       = '13';
% DCM
data(num).DCMfile       = {'.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 6;
data(num).SleepTime     = 22;
data(num).dlmo          = 18.58;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = []; 


%% 25 - COF045
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 70;     % Kg
data(num).height        = 168;    % Cm
data(num).id            = 'COF045';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = 'E0734';
data(num).THKDose       = '';
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s1252-0005-00001-000224-01.nii,1'};
data(num).PTfile     = {{'spmeeg_COF045_P1_20170327_092625_10.mat'} ... %P1
{'spmeeg_COF045_P2_20170327_092625_13.mat'} ...%P2
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF045_S1_TMSEEG_20170407_090318_6.mat'} ... %S1
    {'spmeeg_COF045_S2_TMSEEG_20170407_152650_5.mat'} ... %S2
    {'spmeeg_COF045_S3_TMSEEG_20170407_185938_2.mat'} ... %S3
    {'spmeeg_COF045_S4_TMSEEG_20170407_212026_2.mat'} ... %S4
    {'spmeeg_COF045_S5_TMSEEG_20170408_011648_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 28;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1542';
data(num).T1w           = '2';
data(num).PDw           = '4';
data(num).MTw           = '6';
data(num).B1            = '8';
data(num).B0magn        = '9';
data(num).B0phase       = '10';
% DCM
data(num).DCMfile       = {'merged_fmbLAPefdfMspmeeg_COF045_S1_TMSEEG_20170407_090318_6.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 7;
data(num).SleepTime     = 23;
data(num).dlmo          = 20.24;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = [3 5 6 7 8 9]; 


%% 26 - COF048
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 70;     % Kg
data(num).height        = 174;    % Cm
data(num).id            = 'COF048';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = '';
data(num).THKDose       = '';
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s1511-0002-00001-000224-01,1'};
data(num).PTfile     = {{'spmeeg_COF048_P1_20171019_095336_3.mat'} ... %P1
{'spmeeg_COF048_P2_20171019_095336_9.mat'} ...%P2
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF048_S1_TMSEEG_20171122_091530_6.mat'} ... %S1
    {'spmeeg_COF048_S2_TMSEEG_20171122_153735_2.mat'} ... %S2
    {'spmeeg_COF048_S3_TMSEEG_20171122_185205_3.mat'} ... %S3
    {'spmeeg_COF048_S4_TMSEEG_20171122_211916_2.mat'} ... %S4
    {'spmeeg_COF048_S5_TMSEEG_20171123_010640_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 9;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1511';
data(num).T1w           = '5';
data(num).PDw           = '7';
data(num).MTw           = '9';
data(num).B1            = '11';
data(num).B0magn        = '12';
data(num).B0phase       = '13';
% DCM
data(num).DCMfile       = {'.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 7;
data(num).SleepTime     = 23.5;
data(num).dlmo          = 22.66;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = []; 


%% 27 - COF049
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 70;     % Kg
data(num).height        = 168;    % Cm
data(num).id            = 'COF049';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = 'E0590';
data(num).THKDose       = 192;
% TMS Pretests
data(num).PTname      = ['P1' ; 'P2' ; 'P3'];
data(num).MRI           = {'s1261-0005-00001-000224-01.nii,1'};
data(num).PTfile     = {{'spmeeg_COF049_P1_20170330_091620_4.mat'} ... %P1
{'spmeeg_COF049_P2_20170330_091620_7.mat'} ...%P2
{'spmeeg_COF049_P3_20170330_091620_11.mat'} ...%P3
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF049_S1_TMSEEG_20170411_092315_6.mat'} ... %S1
    {'spmeeg_COF049_S2_TMSEEG_20170411_160731_2.mat'} ... %S2
    {'spmeeg_COF049_S3_TMSEEG_20170411_193020_2.mat'} ... %S3
    {'spmeeg_COF049_S4_TMSEEG_20170411_215748_2.mat'} ... %S4
    {'spmeeg_COF049_S5_TMSEEG_20170412_013755_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 27;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1357';
data(num).T1w           = '2';
data(num).PDw           = '4';
data(num).MTw           = '6';
data(num).B1            = '8';
data(num).B0magn        = '9';
data(num).B0phase       = '10';
% DCM
data(num).DCMfile       = {'merged_fmbLAPefdfMspmeeg_COF049_S1_TMSEEG_20170411_092315_6.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 7.5;
data(num).SleepTime     = 22.5;
data(num).dlmo          = 20.09;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = [1 4 5 7 8]; 


%% 28 - COF050
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 74;     % Kg
data(num).height        = 162;    % Cm
data(num).id            = 'COF050';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = 'E0732';
data(num).THKDose       = '';
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s1262-0005-00001-000224-01.nii,1'};
data(num).PTfile     = {{'spmeeg_COF050_P1_20170330_113515_3.mat'} ... %P1
{'spmeeg_COF050_P2_20170330_144015_7.mat'} ...%P2
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF050_S1_TMSEEG_20170411_073613_4.mat'} ... %S1
    {'spmeeg_COF050_S2_TMSEEG_20170411_145452.mat'} ... %S2
    {'spmeeg_COF050_S3_TMSEEG_20170411_185441_2.mat'} ... %S3
    {'spmeeg_COF050_S4_TMSEEG_20170411_203320_3.mat'} ... %S4
    {'spmeeg_COF050_S5_TMSEEG_20170411_225344_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 39;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1262';
data(num).T1w           = '8';
data(num).PDw           = '10';
data(num).MTw           = '12';
data(num).B1            = '14';
data(num).B0magn        = '15';
data(num).B0phase       = '16';
% DCM
data(num).DCMfile       = {'merged_fmbLAPefdfMspmeeg_COF050_S1_TMSEEG_20170411_073613_4.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% DLMO
data(num).wkTime        = 6.5;
data(num).SleepTime     = 23.5;
data(num).dlmo          = 20.93;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = [2 3 4 7 8 9]; 


%% 29 - COF055
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 81;     % Kg
data(num).height        = 176;    % Cm
data(num).id            = 'COF055';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = 'E0721';
data(num).THKDose       = 194;
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s1307-0005-00001-000224-01.nii,1'};
data(num).PTfile        = {{'spmeeg_COF055_P1_20170511_092326_3.mat'} ...   %P1
    {'spmeeg_COF055_P2_20170511_092326_6.mat'} ...                        %P2
    };
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF055_S1_TMSEEG_20170525_090642_4.mat'} ... %S1
    {'spmeeg_COF055_S2_TMSEEG_20170525_151510_2.mat'} ... %S2
    {'spmeeg_COF055_S3_TMSEEG_20170525_184249_2.mat'} ... %S3
    {'spmeeg_COF055_S4_TMSEEG_20170525_210427_3.mat'} ... %S4
    {'spmeeg_COF055_S5_TMSEEG_20170526_005521_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 9;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1307';
data(num).T1w           = '8';
data(num).PDw           = '10';
data(num).MTw           = '12';
data(num).B1            = '14';
data(num).B0magn        = '15';
data(num).B0phase       = '16';
% DCM
data(num).DCMfile       = {'merged_fmbLAPefdfMspmeeg_COF055_S1_TMSEEG_20170525_090642_4.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% DLMO
data(num).wkTime        = 6.75;
data(num).SleepTime     = 23.75;
data(num).dlmo          = 19.57;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = []; 


%% 30 - COF057
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 96;     % Kg
data(num).height        = 185;    % Cm
data(num).id            = 'COF057';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = 'E0705';
data(num).THKDose       = 188;
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s1298-0005-00001-000224-01.nii,1'};
data(num).PTfile        = {{'spmeeg_COF057_P1_20170509_084900_4.mat'} ...   %P1
    {'spmeeg_COF057_P2_20170509_084900_7.mat'} ...                        %P2
    };
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF057_S1_TMSEEG_20170525_074437_3.mat'} ... %S1
    {'spmeeg_COF057_S2_TMSEEG_20170525_134900_2.mat'} ... %S2
    {'spmeeg_COF057_S3_TMSEEG_20170525_160503_2.mat'} ... %S3
    {'spmeeg_COF057_S4_TMSEEG_20170525_194739_2.mat'} ... %S4
    {'spmeeg_COF057_S5_TMSEEG_20170526_000043_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 28;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1338';
data(num).T1w           = '2';
data(num).PDw           = '4';
data(num).MTw           = '6';
data(num).B1            = '8';
data(num).B0magn        = '9';
data(num).B0phase       = '10';
% DCM
data(num).DCMfile       = {'merged_fmbLAPefdfMspmeeg_COF057_S1_TMSEEG_20170525_074437_3.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% DLMO
data(num).wkTime        = 5.75;
data(num).SleepTime     = 22.5;
data(num).dlmo          = 21.6;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = []; 


%% 31 - COF059
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 70;     % Kg
data(num).height        = 172;    % Cm
data(num).id            = 'COF059';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = 'E0726';
data(num).THKDose       = 185;
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s1316-0005-00001-000224-01.nii,1'};
data(num).PTfile        = {{'spmeeg_COF059_P1_20170519_095739_8.mat'} ...   %P1
    {'spmeeg_COF059_P2_20170519_095739_11.mat'} ...                        %P2
    };
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF059_S1_TMSEEG_20170601_075520_3.mat'} ... %S1
    {'spmeeg_COF059_S2_TMSEEG_20170601_094849_2.mat'} ... %S2
    {'spmeeg_COF059_S3_TMSEEG_20170601_162022_2.mat'} ... %S3
    {'spmeeg_COF059_S4_TMSEEG_20170601_194712_2.mat'} ... %S4
    {'spmeeg_COF059_S5_TMSEEG_20170601_221107_3.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 5;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1316';
data(num).T1w           = '8';
data(num).PDw           = '10';
data(num).MTw           = '12';
data(num).B1            = '14';
data(num).B0magn        = '15';
data(num).B0phase       = '16';
% DCM
data(num).DCMfile       = {'merged_fmbLAPefdfMspmeeg_COF059_S1_TMSEEG_20170601_075520_3.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% DLMO
data(num).wkTime        = 7;
data(num).SleepTime     = 23;
data(num).dlmo          = 20.84;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = []; 


%% 32 - COF060
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 75;     % Kg
data(num).height        = 172;    % Cm
data(num).id            = 'COF060';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = 'E0717';
data(num).THKDose       = '';
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s1273-0005-00001-000224-01.nii,1'};
data(num).PTfile     = {{'spmeeg_COF060_P1_20170421_100614_3.mat'} ... %P1
{'spmeeg_COF060_P2_20170421_100614_6.mat'} ...%P2
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF060_S1_TMSEEG_20170428_094156_4.mat'} ... %S1
    {'spmeeg_COF060_S2_TMSEEG_20170428_160611_2.mat'} ... %S2
    {'spmeeg_COF060_S3_TMSEEG_20170428_192237_2.mat'} ... %S3
    {'spmeeg_COF060_S4_TMSEEG_20170428_214453_2.mat'} ... %S4
    {'spmeeg_COF060_S5_TMSEEG_20170429_014006_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 28;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1490';
data(num).T1w           = '2';
data(num).PDw           = '4';
data(num).MTw           = '6';
data(num).B1            = '8';
data(num).B0magn        = '9';
data(num).B0phase       = '10';
% DCM
data(num).DCMfile       = {'merged_fmbLAPefdfMspmeeg_COF060_S1_TMSEEG_20170428_094156_4.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% DLMO
data(num).wkTime        = 7.5;
data(num).SleepTime     = 23;
data(num).dlmo          = 20.91;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = []; 


%% 33 - COF061
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 62;     % Kg
data(num).height        = 175;    % Cm
data(num).id            = 'COF061';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = '';
data(num).THKDose       = '';
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s1504-0002-00001-000224-01,1'};
data(num).PTfile        = {{'spmeeg_COF061_P1_20171012_202107_9.mat'} ...   %P1
    {'spmeeg_COF061_P2_20171012_202107_13.mat'} ...                        %P2
    };
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF061_S1_TMSEEG_171118_094527_2.mat'} ... %S1
    {'spmeeg_COF061_S2_TMSEEG_20171118_151121_2.mat'} ... %S2
    {'spmeeg_COF061_S3_TMSEEG_20171118_185145_2.mat'} ... %S3
    {'spmeeg_COF061_S4_TMSEEG_20171118_212645_2.mat'} ... %S4
    {'spmeeg_COF061_S5_TMSEEG_20171119_011103_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 28;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1537';
data(num).T1w           = '3';
data(num).PDw           = '5';
data(num).MTw           = '7';
data(num).B1            = '9';
data(num).B0magn        = '10';
data(num).B0phase       = '11';
% DCM
data(num).DCMfile       = {'.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% DLMO
data(num).wkTime        = 7;
data(num).SleepTime     = 22;
data(num).dlmo          = 19.86;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = [];


%% 34 - COF062
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 73;     % Kg
data(num).height        = 170;    % Cm
data(num).id            = 'COF062';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = 'E0750';
data(num).THKDose       = '';
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s1451-0002-00001-000224-01,1'};
data(num).PTfile     = {{'spmeeg_COF062_P1_20170830_092928_5.mat'} ... %P1
{'spmeeg_COF062_P2_20170830_092928_7.mat'} ...%P2
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF062_S1_EEGTMS_20170926_064855_8.mat'} ... %S1
    {'spmeeg_COF062_S2_TMSEEG_20170926_100443_3.mat'} ... %S2
    {'spmeeg_COF062_S3_TMSEEG_20170926_181522_2.mat'} ... %S3
    {'spmeeg_COF062_S4_TMSEEG_20170926_202013_2.mat'} ... %S4
    {'spmeeg_COF062_S5_TMSEEG_20170926_222738_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 19;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1451';
data(num).T1w           = '5';
data(num).PDw           = '7';
data(num).MTw           = '9';
data(num).B1            = '11';
data(num).B0magn        = '12';
data(num).B0phase       = '13';
% DCM
data(num).DCMfile       = {'.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% DLMO
data(num).wkTime        = 6;
data(num).SleepTime     = 23;
data(num).dlmo          = 19.64;        % Computed with ascending level set to 5 because 2.3 did not seem correct
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [9];
data(num).pupil_ER      = [9];


%% 35 - COF063
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 75;     % Kg
data(num).height        = 165;    % Cm
data(num).id            = 'COF063';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = 'E0725';
data(num).THKDose       = 190;
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s1311-0005-00001-000224-01.nii,1'};
data(num).PTfile        = {{'spmeeg_COF063_P1_20170516_170533_21.mat'} ...   %P1
    {'spmeeg_COF063_P2_20170516_170533_25.mat'} ...                        %P2
    };
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF063_S1_TMSEEG_20170520_065515_5.mat'} ... %S1
    {'spmeeg_COF063_S2_TMSEEG-20170520_093205_2.mat'} ... %S2
    {'spmeeg_COF063_S3_TMSEEG_20170520_172301_3.mat'} ... %S3
    {'spmeeg_COF063_S4_TMSEEG_20170520_192821_2.mat'} ... %S4
    {'spmeeg_COF063_S5_TMSEEG_20170520_220508_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 5;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1507';
data(num).T1w           = '2';
data(num).PDw           = '4';
data(num).MTw           = '6';
data(num).B1            = '8';
data(num).B0magn        = '9';
data(num).B0phase       = '10';
% DCM
data(num).DCMfile       = {'merged_fmbLAPefdfMspmeeg_COF063_S1_TMSEEG_20170520_065515_5.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% DLMO
data(num).wkTime        = 6.5;
data(num).SleepTime     = 22;
data(num).dlmo          = 18.61;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = []; 


%% 36 - COF065
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 66;     % Kg
data(num).height        = 171;    % Cm
data(num).id            = 'COF065';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = 'E0742';
data(num).THKDose       = '';
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s1461-0002-00001-000224-01,1'};
data(num).PTfile        = {{'spmeeg_COF065_P1_20170906_095610_4.mat'} ...   %P1
    {'spmeeg_COF065_P2_20170906_095610_7.mat'} ...                        %P2
    };
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF065_S1_TMSEEG_20170919_080303_5.mat'} ... %S1
    {'spmeeg_COF065_S2_EEGTMS_20170919_102854_2.mat'} ... %S2
    {'spmeeg_COF065_S3_TMSEEG_20170919_170818_2.mat'} ... %S3
    {'spmeeg_COF065_S4_EEGTMS_20170919_213223_2.mat'} ... %S4
    {'spmeeg_COF065_S5_TMEEEG_20170919_230916_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 5;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1461';
data(num).T1w           = '5';
data(num).PDw           = '7';
data(num).MTw           = '9';
data(num).B1            = '11';
data(num).B0magn        = '12';
data(num).B0phase       = '13';
% DCM
data(num).DCMfile       = {'.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% DLMO
data(num).wkTime        = 6.5;
data(num).SleepTime     = 22.75;
data(num).dlmo          = 20.01;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = [];


%% 37 - COF068
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 60;     % Kg
data(num).height        = 176;    % Cm
data(num).id            = 'COF068';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = '';
data(num).THKDose       = '';
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s1520-0002-00001-000224-01,1'};
data(num).PTfile     = {{'spmeeg_COF068_P1_20171023_182821_4.mat'} ... %P1
{'spmeeg_COF068_P2_20171023_182821_8.mat'} ...%P2
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF068_S1_TMSEEG_20171028_093218_5.mat'} ... %S1
    {'spmeeg_COF068_S2_TMSEEG_20171028_161127_3.mat'} ... %S2
    {'spmeeg_COF068_S3_EEGTMS_20171028_193629_4.mat'} ... %S3
    {'spmeeg_COF068_S4_TMSEEG_20171028_215946_2.mat'} ... %S4
    {'spmeeg_COF068_S5_TMSEEG_20171029_015800_2.mat'} ... %S5
       };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 5;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1541';
data(num).T1w           = '5';
data(num).PDw           = '7';
data(num).MTw           = '9';
data(num).B1            = '11';
data(num).B0magn        = '12';
data(num).B0phase       = '13';
% DCM
data(num).DCMfile       = {'.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 7.5;
data(num).SleepTime     = 22.5;
data(num).dlmo          = 22.31;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = []; 


%% 38 - COF070
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 60;     % Kg
data(num).height        = 151;    % Cm
data(num).id            = 'COF070';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = 'E0743';
data(num).THKDose       = '';
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s1430-0002-00001-000224-01,1'};
data(num).PTfile        = {{'spmeeg_COF070_P1_20170817_153127_3.mat'} ...   %P1
    {'spmeeg_COF070_P2_20170817_153127_8.mat'} ...                        %P2
    };
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF070_S1_TMSEEG_20170914_090234_3.mat'} ... %S1
    {'spmeeg_COF070_S2_TMSEEG_20170914_152014_2.mat'} ... %S2
    {'spmeeg_COF070_S3_TMSEEG_20170914_190055_2.mat'} ... %S3
    {'spmeeg_COF070_EEGTMS_S4_20170914_211249_2.mat'} ... %S4
    {'spmeeg_COF070_S5_EEGTMS_20170915_012309_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 29;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1527';
data(num).T1w           = '2';
data(num).PDw           = '4';
data(num).MTw           = '6';
data(num).B1            = '8';
data(num).B0magn        = '9';
data(num).B0phase       = '10';
% DCM
data(num).DCMfile       = {'.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 7;
data(num).SleepTime     = 22.75;
data(num).dlmo          = 20.93;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [4];
data(num).pupil_ER      = []; 


%% 39 - COF073
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 80;     % Kg
data(num).height        = 172;    % Cm
data(num).id            = 'COF073';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = '';
data(num).THKDose       = '';
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2' ; 'P3'];
data(num).MRI           = {'s1524-0002-00001-000224-01,1'};
data(num).PTfile     = {{'spmeeg_COF073_P1_20171025_091403_6.mat'} ... %P1
{'spmeeg_COF073_P2_20171025_091403_11.mat'} ... %P2
{'spmeeg_COF073_P3_20171025_091403_14.mat'}... %P3
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF073_S1_TMSEEG_20171108_074606_6.mat'} ... %S1
    {'spmeeg_COF073_S2_EEGTMS_20171108_110908_2.mat'} ... %S2
    {'spmeeg_COF073_S3_TMSEEG_20171108_164904_2.mat'} ... %S3
    {'spmeeg_COF073_S4_TMSEEG_20171108_201950_2.mat'} ... %S4
    {'spmeeg_COF073_S5_TMSEEG_20171108_225931_3.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 29;     % 19 second best, 4 third best
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1524';
data(num).T1w           = '5';
data(num).PDw           = '7';
data(num).MTw           = '9';
data(num).B1            = '11';
data(num).B0magn        = '12';
data(num).B0phase       = '13';
% DCM
data(num).DCMfile       = {'.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 6.5;
data(num).SleepTime     = 23;
data(num).dlmo          = 20.90;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = []; 


%% 40 - COF075
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 65;     % Kg
data(num).height        = 160;    % Cm
data(num).id            = 'COF075';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = '';
data(num).THKDose       = '';
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s1464-0002-00001-000224-01,1'};
data(num).PTfile     = {{'spmeeg_COF075_P1_20170908_141801_5.mat'} ... %P1
{'spmeeg_COF075_P2_20170908_141801_10.mat'} ...%P2
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF075_S1_TMSEEG_20170919_093441_4.mat'} ... %S1
    {'spmeeg_COF075_S2_EEGTMS_20170919_154133_2.mat'} ... %S2
    {'spmeeg_COF075_S3_TMSEEG_20170919_194706_2.mat'} ... %S3
    {'spmeeg_COF075_S4_EEGTMS_20170919_222312_2.mat'} ... %S4
    {'spmeeg_COF075_S5_EEGTMS_20170920_014128_2.mat'} ... %S5
       };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 5;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1464';
data(num).T1w           = '5';
data(num).PDw           = '7';
data(num).MTw           = '9';
data(num).B1            = '11';
data(num).B0magn        = '12';
data(num).B0phase       = '13';
% DCM
data(num).DCMfile       = {'.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 7.5;
data(num).SleepTime     = 22.5;
data(num).dlmo          = 19.60;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = []; 


%% 41 - COF080
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 51;     % Kg
data(num).height        = 160;    % Cm
data(num).id            = 'COF080';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = '';
data(num).THKDose       = '';
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2' ; 'P3'];
data(num).MRI           = {'s1481-0002-00001-000224-01,1'};
data(num).PTfile     = {{'spmeeg_COF080_P1_20170923_101308_4.mat'} ... %P1
{'spmeeg_COF080_P2_20170923_101308_7.mat'} ...%P2
{'spmeeg_COF080_P3_20170923_101308_9.mat'} ...%P3
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF080_S1_EEGTMS_20170926_091928_3.mat'} ... %S1
    {'spmeeg_COF080_S2_EEGTMS_20170926_151532_3.mat'} ... %S2
    {'spmeeg_COF080_S3_TMSEEG_20170926_191402_4.mat'} ... %S3
    {'spmeeg_COF080_S4_TMSEEG_20170926_211310_2.mat'} ... %S4
    {'spmeeg_COF080_S5_TMSEEG_20170927_010803_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 28;           % 17/29 - To check
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1481';
data(num).T1w           = '5';
data(num).PDw           = '7';
data(num).MTw           = '9';
data(num).B1            = '11';
data(num).B0magn        = '12';
data(num).B0phase       = '13';
% DCM
data(num).DCMfile       = {'.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 7;
data(num).SleepTime     = 22.75;
data(num).dlmo          = 21.79;        % Computed with ascending level set to 5 because 2.3 did not seem correct
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = []; 


%% 42 - COF083
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 88;     % Kg
data(num).height        = 182;    % Cm
data(num).id            = 'COF083';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = '';
data(num).THKDose       = '';
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s1525-0002-00001-000224-01,1'};
data(num).PTfile        = {{'spmeeg_COF083_P1_20171025_125825_14.mat'} ...   %P1
    {'spmeeg_COF083_P2_20171025_125825_19.mat'} ...                        %P2
    };
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF083_S1_TMSEEG_20171108_103918_6.mat'} ... %S1
    {'spmeeg_COF083_S2_EEGTMS_20171108_153114_2.mat'} ... %S2
    {'spmeeg_COF083_S3_TMSEEG_20171108_192931_2.mat'} ... %S3
    {'spmeeg_COF083_S4_TMSEEG_20171108_215823_2.mat'} ... %S4
    {'spmeeg_COF083_S5_TMSEEG_20171109_015349_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 27; % But really far... 4 second best except for S1
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1525';
data(num).T1w           = '5';
data(num).PDw           = '7';
data(num).MTw           = '9';
data(num).B1            = '11';
data(num).B0magn        = '12';
data(num).B0phase       = '13';
% DCM
data(num).DCMfile       = {'.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 6.5;
data(num).SleepTime     = 22.5;
data(num).dlmo          = 19.48;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = []; 


%% 43 - COF085
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 63;     % Kg
data(num).height        = 158;    % Cm
data(num).id            = 'COF085';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = 'E0753';
data(num).THKDose       = '';
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2' ; 'P3'];
data(num).MRI           = {'s1505-0019-00001-000224-01,1'};
data(num).PTfile     = {{'spmeeg_COF085_P1_20171013_102448_5.mat'} ... %P1
{'spmeeg_COF085_P2_20171013_102448_10.mat'} ...%P2
{'spmeeg_COF085_P3_20171020_101252_9.mat'}
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF085_S1_TMSEEG_20171028_071421_7.mat'} ... %S1
    {'spmeeg_COF085_S2_TMSEEG_20171028_103922_2.mat'} ... %S2
    {'spmeeg_COF085_S3_EEGTMS_20171028_173415_2.mat'} ... %S3
    {'spmeeg_COF085_S4_TMSEEG_20171028_193609_2.mat'} ... %S4
    {'spmeeg_COF085_S5_TMSEEG_20171028_231020_2.mat'} ... %S5
           };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 5;        % But bad for 4th session - To check
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1505';
data(num).T1w           = '5';
data(num).PDw           = '7';
data(num).MTw           = '9';
data(num).B1            = '11';
data(num).B0magn        = '12';
data(num).B0phase       = '13';
% DCM
data(num).DCMfile       = {'.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 6.5;
data(num).SleepTime     = 22;
data(num).dlmo          = 21.28;        % WARNING: Lots of 0 values for melatonin concentration!
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [7];
data(num).kss_tasks     = [];
data(num).kss_tms       = [4];
data(num).vas_alone     = [7];
data(num).vas_tasks     = [];
data(num).vas_tms       = [4];
% Pupillometry
data(num).pupil_cal     = [5];
data(num).pupil_ER      = [5];


%% 44 - COF086
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 65;     % Kg
data(num).height        = 153;    % Cm
data(num).id            = 'COF086';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = '';
data(num).THKDose       = '';
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s1528-0002-00001-000224-01,1'};
data(num).PTfile     = {{'spmeeg_COF086_P1_20171027_095822_3.mat'} ... %P1
{'spmeeg_COF086_P2_20171027_095822_5.mat'} ...%P2
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF086_S1_TMSEEG_20171118_071637_5.mat'} ... %S1
    {'spmeeg_COF086_S2_TMSEEG_20171118_141354_2.mat'} ... %S2
    {'spmeeg_COF086_S3_TMSEEG_20171118_163610_2.mat'} ... %S3
    {'spmeeg_COF086_S4_TMSEEG_20171118_200426_2.mat'} ... %S4
    {'spmeeg_COF086_S5_TMSEEG_20171118_222346_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 28;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1528';
data(num).T1w           = '5';
data(num).PDw           = '7';
data(num).MTw           = '9';
data(num).B1            = '11';
data(num).B0magn        = '12';
data(num).B0phase       = '13';
% DCM
data(num).DCMfile       = {'.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 6;
data(num).SleepTime     = 22.5;
data(num).dlmo          = 18.17;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = []; 


%% 45 - COF088
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 54;     % Kg
data(num).height        = 163;    % Cm
data(num).id            = 'COF088';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = '';
data(num).THKDose       = '';
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s1557-0002-00001-000224-01,1'};
data(num).PTfile     = {{'spmeeg_COF088_P1_20171111_123937_8.mat'} ... %P1
{'spmeeg_COF088_P2_20171111_123937_14.mat'} ...%P2
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF088_S1_TMSEEG_20171128_085306_4.mat'} ... %S1
    {'spmeeg_COF088_S2_TMSEEG_20171128_154237_3.mat'} ... %S2
    {'spmeeg_COF088_S3_TMSEEG_20171128_190245_2.mat'} ... %S3
    {'spmeeg_COF088_S4_TMSEEG_20171128_212413_2.mat'} ... %S4
    {'spmeeg_COF088_S5_TMSEEG_20171129_012304_3.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 4;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1557';
data(num).T1w           = '5';
data(num).PDw           = '7';
data(num).MTw           = '9';
data(num).B1            = '11';
data(num).B0magn        = '12';
data(num).B0phase       = '13';
% DCM
data(num).DCMfile       = {'.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 7;
data(num).SleepTime     = 22;
data(num).dlmo          = 18.88;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = []; 


%% 46 - COF089
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 73;     % Kg
data(num).height        = 186;    % Cm
data(num).id            = 'COF089';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = '';
data(num).THKDose       = '';
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2' ; 'P3'];
data(num).MRI           = {'s1556-0002-00001-000224-01,1'};
data(num).PTfile     = {{'spmeeg_COF089_P1_20171111_093652_4.mat'} ... %P1
{'spmeeg_COF089_P2_20171111_093652_8.mat'} ...%P2
{'spmeeg_COF089_P3_20171111_093652_12.mat'} ...%P2
};
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF089_S1_TMSEEG_20171128_071032_4.mat'} ... %S1
    {'spmeeg_COF089_S2_TMSEEG_20171128_145556_2.mat'} ... %S2
    {'spmeeg_COF089_S3_TMSEEG_20171128_174151_2.mat'} ... %S3
    {'spmeeg_COF089_S4_TMSEEG_20171128_200656_2.mat'} ... %S4
    {'spmeeg_COF089_S5_TMSEEG_20171128_221251_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 18;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1556';
data(num).T1w           = '5';
data(num).PDw           = '7';
data(num).MTw           = '9';
data(num).B1            = '11';
data(num).B0magn        = '12';
data(num).B0phase       = '13';
% DCM
data(num).DCMfile       = {'.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 6;
data(num).SleepTime     = 22;
data(num).dlmo          = 18.97;
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = []; 


%% 47 - COF091
num                     = num + 1;
data(num).num           = num;

% Informations
data(num).weight        = 66;     % Kg
data(num).height        = 172;    % Cm
data(num).id            = 'COF091';
data(num).dir           = fullfile(pathtodata,data(num).id);
% PET Scans
data(num).THKID         = '';
data(num).THKDose       = '';
% TMS Pretests
data(num).PTname        = ['P1' ; 'P2'];
data(num).MRI           = {'s1554-0002-00001-000224-01,1'};
data(num).PTfile     = {{'spmeeg_COF091_P1_20171110_102804_4.mat'} ... %P1
{'spmeeg_COF091_P2_20171110_102804_6.mat'} ...%P2
};  
% TMS Sessions
data(num).sessname      = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).sspname       = ['S1' ; 'S2' ; 'S3' ; 'S4' ; 'S5'];
data(num).EEGfile     = {{'spmeeg_COF091_S1_TMSEEG_20171122_072716_4.mat'} ... %S1
    {'spmeeg_COF091_S2_TMSEEG_20171122_103419_2.mat'} ... %S2
    {'spmeeg_COF091_S3_EEGTMS_20171122_163049_2.mat'} ... %S3
    {'spmeeg_COF091_S4_TMSEEG_20171122_194735_2.mat'} ... %S4
    {'spmeeg_COF091_S5_TMSEEG_20171122_221948_2.mat'} ... %S5
    };
data(num).Sname         = {'LSMA'};
data(num).channel       = '';
data(num).closestE      = 17;
data(num).EEGsession    = [];
% MRI MPM
data(num).LongMRI       = '1554';
data(num).T1w           = '5';
data(num).PDw           = '7';
data(num).MTw           = '9';
data(num).B1            = '11';
data(num).B0magn        = '12';
data(num).B0phase       = '13';
% DCM
data(num).DCMfile       = {'.mat'};
% BL
data(num).FPL_OPL = [];
data(num).Sleep_cycles  = {};
% Circadian Phase
data(num).wkTime        = 6;
data(num).SleepTime     = 22.75;
data(num).dlmo          = 19.84;        % Computed with ascending value set to 5 because 2.3 did not seem correct
% Cognitive tasks
data(num).twoback       = [];
data(num).threeback     = [];
data(num).sart          = [];
data(num).pvt           = [];
data(num).kss_alone     = [];
data(num).kss_tasks     = [];
data(num).kss_tms       = [];
data(num).vas_alone     = [];
data(num).vas_tasks     = [];
data(num).vas_tms       = [];
% Pupillometry
data(num).pupil_cal     = [];
data(num).pupil_ER      = []; 


%% ********************************************************************* %%


return