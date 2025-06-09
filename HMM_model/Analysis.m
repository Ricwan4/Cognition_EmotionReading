%% get significant
addpath(genpath('src\'))

addpath(genpath('images\'))

% load("HMM_output\Trained_Model.mat")
Nsubjects = length(new_data);
Nstimuli = height(new_data);
HEM_K = length(cogroup_hmms{1}.hmms);
cg_hmms = cogroup_hmms;
grps = cogroup_hmms{1}.groups;

for i=1:Nsubjects
  for j=1:Nstimuli
    for k=1:HEM_K
      if ~isempty(new_data{j,i})
        % compute LL: if more than one sample for this stimuli, then average.
        % the LL for each subject of each situli under each model
        LL(i,j,k) = mean(vbhmm_ll(cg_hmms{j}.hmms{k}, new_data{j,i}));
      else
        LL(i,j,k) = 0;
      end

    end
  end
end

% (assumes K=2)
% LL of subject under model 1 and 2
LL1 = sum(LL(:,:,1),2);  % sum over stimuli
LL2 = sum(LL(:,:,2),2);

% compute AB scale
% since some stimuli are missing, then normalize
AB = (LL1 - LL2) ./ (abs(LL1) + abs(LL2));

stim_LL1 = [];
stim_LL2 = [];
  
for j=1:Nstimuli
    % get LL for group 1 under models 1 and 2
    stim_LL1(:,j,1) = LL(grps{1}, j, 1);
    stim_LL1(:,j,2) = LL(grps{1}, j, 2);
    
    % get LL for group 2 under models 2 and 1
    stim_LL2(:,j,1) = LL(grps{2}, j, 2);
    stim_LL2(:,j,2) = LL(grps{2}, j, 1);
end

fprintf('=== group difference (overall) ===\n');
  
% average LL over stimulus for group 1 under models 1 and 2
LL1_grp1 = mean(stim_LL1(:,:,1),1); % subjects in group under group 1 model
LL1_grp2 = mean(stim_LL1(:,:,2),1); % subjects in group under group 2 model
  
[h1,p1,ci1,stats1] = ttest(LL1_grp1, LL1_grp2, 'tail', 'both');
deffect = computeCohen_d(LL1_grp1,LL1_grp2, 'paired');
fprintf('model1 vs model2 using grp1 data: p=%0.4g, t(%d)=%0.4g, cohen_d=%0.4g\n', ...
p1, stats1.df, stats1.tstat, deffect);
  
% average LL over stimulus for group 2 under models 1 and 2
LL2_grp1 = mean(stim_LL2(:,:,1),1); % subjects in group under group 1 model
LL2_grp2 = mean(stim_LL2(:,:,2),1);% subjects in group under group 2 model
   
[h2,p2,ci2,stats2] = ttest(LL2_grp1, LL2_grp2, 'tail', 'both');
deffect = computeCohen_d(LL2_grp1,LL2_grp2, 'paired');
fprintf('model1 vs model2 using grp2 data: p=%0.4g, t(%d)=%0.4g, cohen_d=%0.4g\n', ...
    p2, stats2.df, stats2.tstat, deffect);

clear LL1 LL2 AB stim_LL1 stim_LL2 LL1_grp1 LL1_grp2 h1 p1 ci1 stats1 deffect LL2_grp1 LL2_grp2 h2 p2 ci2 stat2 cg_hmms HEM_K

%% A-B scale
% Initialize variables
Nsubjects = length(new_data);
Nstimuli = height(new_data);
HEM_K = length(cogroup_hmms{1}.hmms); % Should be 2 for this setup
LL = zeros(Nsubjects, Nstimuli, HEM_K);

% Compute log-likelihood for each subject under each model
for i = 1:Nsubjects
    for j = 1:Nstimuli
        if ~isempty(new_data{j, i})
            for k = 1:HEM_K 
                LL(i, j, k) = mean(vbhmm_ll(cogroup_hmms{j}.hmms{k}, new_data{j, i}));
            end
        else
            LL(i, j, :) = NaN; % Handle missing data
        end
    end
end

% Compute LL1 and LL2 (sum across stimuli for each subject)
LL1 = sum(LL(:, :, 1), 2, 'omitnan');
LL2 = sum(LL(:, :, 2), 2, 'omitnan');

% Compute A-B scale
AB_scale = (LL1 - LL2) ./ (abs(LL1) + abs(LL2));
AB_scale=reshape(AB_scale,40,4);
disp(AB_scale);
clear trials LL1 LL2 k j i HEM_K

%% Entropy

for i=1:size(hmm_subjects,2)
    for j = 1:size(hmm_subjects,1)
      % select the ith HMM
      myhmm = hmm_subjects{j,i};
      % overall entropy
      if ~isempty(myhmm)
      ent(j,i) = vbhmm_entropy(myhmm, alldataC{j,i});
      end
   end
end

%% get Ent
% Initialize categories
categories = ["global", "local"];

% Initialize result
result_ent = [];

% Process each category in order
for i = 1:length(categories)
    % Find indices where mapping_table matches the current category
    indices = strcmp(names, categories(i));
    
    % Extract corresponding values from ent and concatenate
    result_ent = [result_ent; ent(indices)];
end

% mean -----------------------------------------------------

num_subjects = size(ent, 2); 
ent_means = zeros(2, num_subjects); 

% Loop through each subject
for subj = 1:num_subjects
    for cond = 1:length(categories)
        % Find indices where names match current condition (global/local)
        indices = contains(names(:, subj), categories(cond));
        
        % Get values and compute mean
        current_values = ent(indices, subj);
        if ~isempty(current_values)
            ent_means(cond, subj) = mean(current_values);
        else
            ent_means(cond, subj) = NaN;
        end
    end
end

% Transpose so each row is a subject
ent_means = ent_means';
