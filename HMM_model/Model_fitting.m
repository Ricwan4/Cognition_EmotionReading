%% Read file, do settings for EyeHMM

addpath(genpath('src\'))
addpath(genpath('images\'))

tStart = tic;


OPT = struct;
OPT.fixationfile = 'Fixations_Raw.xlsx';
OPT.imgdir = 'images/';
OPT.imgsize = [1024; 768];
OPT.DEBUGMODE = 0;
OPT.S = 1:8;
OPT.imginfo = 'stimuliid_boundaries';
OPT.IMAGE_EXT = 'jpg';

[alldata, SubjNames, TrialNames, StimuliNames] = read_xls_fixations2(OPT.fixationfile);

tmp = unique(cat(1,StimuliNames{:}));
StimuliNamesImages = cell(length(tmp),2);
StimuliNamesImages(:,1) = tmp;
for j=1:length(tmp)
    StimuliNamesImages{j,2} = imread([OPT.imgdir tmp{j}]);
end

% pre-process
alldata = preprocess_fixations(alldata, StimuliNames, OPT.imginfo, StimuliNamesImages);

% NOTE: the order of Stimuli will be alphabetized
[alldataC, TrialNamesC, StimuliNamesC] = sort_coclustering(alldata, TrialNames, StimuliNames);

% remove to avoid confusion
clear alldata TrialNames StimuliNames
close all
[Nstimuli, Nsubjects] = size(alldataC);

% map from StimuliNamesC to entry in StimuliNamesImages
mapC = zeros(1,Nstimuli);
for j=1:Nstimuli
    mapC(j) = find(strcmp(StimuliNamesC{j}, StimuliNamesImages(:,1)));
end
  
% reorder to match StimuliNamesC
StimuliNamesImagesC = StimuliNamesImages(mapC,:);
clear StimuliNamesImages

% remove trial that are bad

for i = 1:52
    for ii = 1:40
        if ~isempty(alldataC{i,ii})
        if length(alldataC{i,ii}{1})<10 % find that trial where only a few fixation
            alldataC{i,ii} = [];
        end
        end
    end
end
clear i ii j tmp tempdir_path

% vbhmm parameters
vbopt.alpha0 = 0.1;
vbopt.mu0    = OPT.imgsize/2;
vbopt.W0     = 0.005;
vbopt.beta0  = 1;
vbopt.v0     = 5;
vbopt.epsilon0 = 0.1;
vbopt.showplot = 0; 
vbopt.bgimage  = '';
vbopt.learn_hyps = 1;
vbopt.seed = 1000; % random seed
vbopt.numtrials = 300;
vbopt.verbose = 1;

emhmm_toggle_color(1);

% Co-clustering parameters
hemopt.tau = round(get_median_length(alldataC(:))); 
hemopt.seed = 1000;  
hemopt.trials = 300;  
hemopt.initmode = 'gmmNew'; 

disp('====================== Settings Complete ======================')
%% HMM FITTING

hmm_subjects = cell(52, 40);

for j = 1:Nstimuli
    fprintf('********************** Stimuli %d *****************************\n', j);
   [tmp] = vbhmm_learn_batch(alldataC(j,:), OPT.S, vbopt);
    hmm_subjects(j,:) = tmp;
end

disp('====================== HMM Complete ======================')

%% Reorganzie individual hmm in question type and emotion type

mapping_table = readtable('mapping_table.xlsx'); % Mapping conditions
names = mapping_table{:,:}; % Extracts all entries into a cell array
names(:,1) = [];
names(1,:) = [];

clear mapping_table

% Create logical masks for cells containing "Happy" and "Sad"

maskGlobal_Happy = contains(names, 'globalHappy');
maskGlobal_Sad = contains(names, 'globalSad');
masklocal_Happy = contains(names, 'localHappy');
masklocal_Sad = contains(names, 'localSad');

GH_hmm = hmm_subjects;
GS_hmm = hmm_subjects;
LH_hmm = hmm_subjects;
LS_hmm = hmm_subjects;

GH_hmm(~maskGlobal_Happy) = {[]};
GS_hmm(~maskGlobal_Sad) = {[]};
LH_hmm(~masklocal_Happy) = {[]};
LS_hmm(~masklocal_Sad) = {[]};

new_hmm = [GH_hmm,GS_hmm,LH_hmm,LS_hmm];

disp('====================== Co-clustering Settings Complete ======================')

clear GH_data GS_data LH_data LS_data GH_hmm GS_hmm LH_hmm LS_hmm maskGlobal_Happy maskGlobal_Sad masklocal_Happy masklocal_Sad

%% Group level (specify to 2, using vhem)

[cogroup_hmms, trials] = vhem_cocluster(new_hmm, 2, [], hemopt);

disp('====================== Co-clustering Complete ======================')

totalTime = toc(tStart);
fprintf('Total running time: %.2f seconds\n', totalTime);

save("HMM_output/Trained_Model.mat", 'StimuliNamesC','cogroup_hmms');
