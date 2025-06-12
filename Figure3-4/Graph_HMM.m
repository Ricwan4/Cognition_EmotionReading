%% get model
addpath(genpath('..\HMM_model\src\'))
addpath(genpath('..\HMM_model\images\'))

load("..\HMM_model\HMM_output\Trained_Model.mat")
emhmm_toggle_color(1);

%% draw group hmm - LOCAL (Figure3)

g = 26;
vhem_plot(cogroup_hmms{g},StimuliNamesC{g},'c',[],{'Global','Local'})

% heatmap - continued at heatmap.ipynb

%% draw group hmm - GLOBAL (Figure4)

g =50;
vhem_plot(cogroup_hmms{g},StimuliNamesC{g},'c',[],{'Global','Local'})

% heatmap - continued at heatmap.ipynb
