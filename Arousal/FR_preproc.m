% Preprocessing of data for using of the AASM Toolbox from Doro
% Written by Franziska Rudzik, 2017.
% Cyclotron Research Centre, University of Liege, Belgium

% INPUT
mainpath  = 'C:\Users\Franziska\Documents\FRANZISKA\SiRENE\FHR\Cooperation\Liege_Toolbox\Data\';
pathin_1    = [mainpath, 'Datasets\ana00\'];                 
pathin_2   = [mainpath, 'Datasets\ana01\']; 
% OUTPUT
outpath  = 'C:\Users\Franziska\Documents\FRANZISKA\SiRENE\FHR\Cooperation\Liege_Toolbox\Data\Output\All\';
nom_start = 7;
nom_end = 4;

% Set relevant PARAMETER
% preprocessing FASST: ./.

%% Preprocessing (created by FR)
% Open all files with Fasst so that a .mat is created from the edf
filePattern = fullfile(pathin_1, '*.edf');
matFiles = dir(filePattern);

for k = 1:length(matFiles)
  baseFileName = matFiles(k).name;
  fullFileName = fullfile(pathin_1, baseFileName);
  crc_eeg_load(fullFileName);
end

% Chunking files
filePattern = fullfile(pathin_2, '*.mat');
matFiles = dir(filePattern);

% Load lights off-text-file
fid = fopen([mainpath, 'lightsOff.csv']);
M = textscan(fid,'%s %d %d','Delimiter',',');
fclose(fid);

% adjust parameter for chunking function of FASST
flags = struct('overwr',[],'numchunk',[]);
flags.numchunk = 1;

% Relative times
for k = 1:length(matFiles)
  baseFileName = matFiles(k).name;
  fullFileName = fullfile(pathin_2, baseFileName);
  D = spm_eeg_load(fullFileName);
  toC = baseFileName(1:length(baseFileName)-4); %depends on file naming: COULD NEED ADJUSTMENT
  argh = find(ismember(M{1},toC));
  start_x = M{2}(argh) * fsample(D);
  stop_x = M{3}(argh) * fsample(D);
  crc_process_chunk(D,start_x,stop_x,flags);
end

