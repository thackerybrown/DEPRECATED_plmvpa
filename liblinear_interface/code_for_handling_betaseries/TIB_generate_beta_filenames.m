function [ ] = TIB_generate_beta_filenames(S)
%TIB_generate_beta_filenames generates indices and lists of betas to be
%used for MVPA analysis using PM_run_mvpa_general
%names must match those of the target condition names identified in the
%_params file. Order must also be the same.
%
%FOR BINARY CLASSIFICATION: assign condition index for 1st + 2nd classes to
%"betaidx{1,1}" + "betaidx{1,2}" - as opposed to whatever cell they might normally be
%associated with (i.e., if you are classifying between classes 2 and 5, we
%don't want to store the indices in cells 2 and 5 for binary.
%
%last modified: 11/13/2014 - add train on one, test on another
%functionality


% load in onsets and regs
onsetsmat = [S.mvpa_dir S.onsets_filename];
load(onsetsmat);

if strcmp(S.trainTask,S.testTask) %if we are doing x-validation classification
%set up a beta index file (equivalent to the onsets element of a normal
%model file). ***IMPORTANT***The index numbers should make sense in light
%of the original model files - i.e. if there were 29 goal_1 onsets, we
%should have 29 1s in the goal_1 idx cell

% betaidx{1,1} = strcmpi('plan1',names);
% betaidx{1,2} = strcmpi('plan2',names);
% betaidx{1,3} = strcmpi('plan3',names);
% betaidx{1,4} = strcmpi('plan4',names);
% betaidx{1,5} = strcmpi('plan5',names);

%% binary
% betaidx{1,1} = strcmpi('cues1',names);
% betaidx{1,2} = strcmpi('cues3',names);

%% 5-way
betaidx{1,1} = strcmpi('cues1',names);
betaidx{1,2} = strcmpi('cues2',names);
betaidx{1,3} = strcmpi('cues3',names);
betaidx{1,4} = strcmpi('cues4',names);
betaidx{1,5} = strcmpi('cues5',names);

% betaidx{1,1} = strcmpi('goal1',names);
% betaidx{1,2} = strcmpi('goal2',names);
% betaidx{1,3} = strcmpi('goal3',names);
% betaidx{1,4} = strcmpi('goal4',names);
% betaidx{1,5} = strcmpi('goal5',names);

%% 10-way
% betaidx{1,1} = strcmpi('cues1_1',names);
% betaidx{1,2} = strcmpi('cues1_2',names);
% betaidx{1,3} = strcmpi('cues2_1',names);
% betaidx{1,4} = strcmpi('cues2_2',names);
% betaidx{1,5} = strcmpi('cues3_1',names);
% betaidx{1,6} = strcmpi('cues3_2',names);
% betaidx{1,7} = strcmpi('cues4_1',names);
% betaidx{1,8} = strcmpi('cues4_2',names);
% betaidx{1,9} = strcmpi('cues5_1',names);
% betaidx{1,10} = strcmpi('cues5_2',names);

%% binary
% bnames{1,1} = 'cues1';
% bnames{1,2} = 'cues3';

%% 5-way
% bnames{1,1} = 'plan1';
% bnames{1,2} = 'plan2';
% bnames{1,3} = 'plan3';
% bnames{1,4} = 'plan4';
% bnames{1,5} = 'plan5';

bnames{1,1} = 'cues1';
bnames{1,2} = 'cues2';
bnames{1,3} = 'cues3';
bnames{1,4} = 'cues4';
bnames{1,5} = 'cues5';

% bnames{1,1} = 'goal1';
% bnames{1,2} = 'goal2';
% bnames{1,3} = 'goal3';
% bnames{1,4} = 'goal4';
% bnames{1,5} = 'goal5';

%% 10-way
% bnames{1,1} = 'cues1_1';
% bnames{1,2} = 'cues1_2';
% bnames{1,3} = 'cues2_1';
% bnames{1,4} = 'cues2_2';
% bnames{1,5} = 'cues3_1';
% bnames{1,6} = 'cues3_2';
% bnames{1,7} = 'cues4_1';
% bnames{1,8} = 'cues4_2';
% bnames{1,9} = 'cues5_1';
% bnames{1,10} = 'cues5_2';

allbetafilenames = dir(fullfile(S.mvpa_dir, 'beta*.img'));

for idx = 1:length(betaidx{1,1})%-1%note, we are filling in the beta file names based on how many betas OF INTEREST we have (length(betaidx)). We don't care about the error reg betas for this analysis
    beta_filenames{idx,1} = [S.mvpa_dir allbetafilenames(idx).name]; %create the analog to "raw_filenames.mat" - i.e. a list of all filenames including the path
end

  


cd(S.mvpa_dir);

savename_betaidx=[S.subj_id '_betas_idx.mat'];
save(savename_betaidx, 'bnames', 'betaidx');

savename_betafnms=['beta_filenames.mat'];
save(savename_betafnms, 'beta_filenames');

else

%% binary
% betaidx{1,1} = strcmpi('cues1',names);
% betaidx{1,2} = strcmpi('cues3',names);

%% 5-way
betaidx_tr{1,1} = strcmpi('goal1',names);%training set beta indices
betaidx_tr{1,2} = strcmpi('goal2',names);
betaidx_tr{1,3} = strcmpi('goal3',names);
betaidx_tr{1,4} = strcmpi('goal4',names);
betaidx_tr{1,5} = strcmpi('goal5',names);

betaidx_te{1,1} = strcmpi('cues1',names);%testing set beta indices
betaidx_te{1,2} = strcmpi('cues2',names);
betaidx_te{1,3} = strcmpi('cues3',names);
betaidx_te{1,4} = strcmpi('cues4',names);
betaidx_te{1,5} = strcmpi('cues5',names);

%betaidx{1,6} = strcmpi('goal_1',names);
%betaidx{1,7} = strcmpi('goal_2',names);
%betaidx{1,8} = strcmpi('goal_3',names);
%betaidx{1,9} = strcmpi('goal_4',names);
%betaidx{1,10} = strcmpi('goal_5',names);

%% 10-way
% betaidx_tr{1,1} = strcmpi('cues1_1',names);
% betaidx_te{1,1} = strcmpi('cues1_2',names);
% betaidx_tr{1,2} = strcmpi('cues2_1',names);
% betaidx_te{1,2} = strcmpi('cues2_2',names);
% betaidx_tr{1,3} = strcmpi('cues3_1',names);
% betaidx_te{1,3} = strcmpi('cues3_2',names);
% betaidx_tr{1,4} = strcmpi('cues4_1',names);
% betaidx_te{1,4} = strcmpi('cues4_2',names);
% betaidx_tr{1,5} = strcmpi('cues5_1',names);
% betaidx_te{1,5} = strcmpi('cues5_2',names);

%% binary
% bnames{1,1} = 'cues1';
% bnames{1,2} = 'cues3';

%% 5-way
bnames_tr{1,1} = 'goal1';
bnames_tr{1,2} = 'goal2';
bnames_tr{1,3} = 'goal3';
bnames_tr{1,4} = 'goal4';
bnames_tr{1,5} = 'goal5';

bnames_te{1,1} = 'cues1';
bnames_te{1,2} = 'cues2';
bnames_te{1,3} = 'cues3';
bnames_te{1,4} = 'cues4';
bnames_te{1,5} = 'cues5';

%bnames{1,6} = 'goal_1';
%bnames{1,7} = 'goal_2';
%bnames{1,8} = 'goal_3';
%bnames{1,9} = 'goal_4';
%bnames{1,10} = 'goal_5';

%% 10-way
% bnames_tr{1,1} = 'cues1_1';
% bnames_te{1,1} = 'cues1_2';
% bnames_tr{1,2} = 'cues2_1';
% bnames_te{1,2} = 'cues2_2';
% bnames_tr{1,3} = 'cues3_1';
% bnames_te{1,3} = 'cues3_2';
% bnames_tr{1,4} = 'cues4_1';
% bnames_te{1,4} = 'cues4_2';
% bnames_tr{1,5} = 'cues5_1';
% bnames_te{1,5} = 'cues5_2';

allbetafilenames = dir(fullfile(S.mvpa_dir, 'beta*.img'));

for idx = 1:length([betaidx_tr{1,1} | betaidx_te{1,1}])%-1%note, we are filling in the beta file names based on how many betas OF INTEREST we have (length(betaidx)). We don't care about the error reg betas for this analysis
    beta_filenames{idx,1} = [S.mvpa_dir allbetafilenames(idx).name]; %create the analog to "raw_filenames.mat" - i.e. a list of all filenames including the path
end

  


cd(S.mvpa_dir);

savename_betaidx_tr=[S.subj_id '_betas_idx_tr.mat'];
savename_betaidx_te=[S.subj_id '_betas_idx_te.mat'];
save(savename_betaidx_tr, 'bnames_tr', 'betaidx_tr');
save(savename_betaidx_te, 'bnames_te', 'betaidx_te');

savename_betafnms=['beta_filenames.mat'];
save(savename_betafnms, 'beta_filenames');    
    
end


end

