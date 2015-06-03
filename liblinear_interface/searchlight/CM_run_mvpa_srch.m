function [Subj res] = CM_run_mvpa_srch(subjArray, classification_name, proclus)
% Runs a searchlight classification analysis on the countermeasures dataset
% Args:
%   subjArray: array of subject numbers to include in the analysis
%
% Returns:
%   Map of highest classification performance

% ID numbers for all subjects in this study
ALL_SUBS = [1,3:10,12:26];
% this selector will be set up by cm_create_sl_selectors and its
% name will be passed to feature_select
SL_SELECTOR_TO_USE = 'trEx_teCm_sl';

% if specific subjects are not specified to classify
if isempty(subjArray)
    subjArray=ALL_SUBS;
end

res = [];
% initialize some classifier and searchlight specific parameter variables
class_args.train_funct_name = 'train_pLR';
class_args.test_funct_name = 'test_pLR';
class_args.penalty=1;

scratch.class_args = class_args;
scratch.perfmet_funct = 'perfmet_auc';
scratch.perfmet_args = struct([]);

statmap_srch_arg.obj_funct = 'statmap_classify';
statmap_srch_arg.scratch = scratch;

% iterate through each subject performing classification
for iSub = subjArray
    % load relevant mvpa and experiment parameters
    [Expt, classArgs, par] = CM_mvpa_params(iSub, 'ret', proclus);
    % create necessary strings
    Expt.sl_selector_to_use = SL_SELECTOR_TO_USE;
    subjId = sprintf(Expt.subjIdFormat,iSub);
    Expt.saveName = cm_make_sl_savename(Expt);
    Expt.impMapDirStr = fullfile(Expt.dir, '/mvpa_results/importance_maps/', Expt.saveName);
    Expt.roiFname = fullfile(Expt.dir,'Masks',Expt.roiName);
    thisMvpaDir = fullfile(Expt.dir, subjId, Expt.mvpaDirStr);
    Expt.subjFname = fullfile(thisMvpaDir, [subjId '_' Expt.roiName '_s8mm_wa.mat']);
    
    
    % training and testing on the same TR, whether you like it or not
    assert(isequal(Expt.trWeights_train, Expt.trWeights_test));
    Expt.trWeights = Expt.trWeights_train;
    
    % load or create the subj structure and output a summary
    if ~exist(Expt.subjFname,'file')
        Subj = CM_mvpa_load_and_preprocess_raw_data(subjId, Expt, nRuns, 1);
    else
        Subj = load(Expt.subjFname,'subj');
        Subj=Subj.subj;
    end
    summarize(Subj);
    
    % load onsets and names of scan files
    onsetsFile = fullfile(thisMvpaDir, Expt.onsetsFname);
    Expt.ons = load(onsetsFile);
    Expt.scanfiles = vertcat(par.swascanfiles.(par.task));
    
    % condense onsets and patterns
    fprintf('\n\nPrepping onsets and patterns for CM%03d...',iSub);
    Subj = cm_condense_onsets_and_patterns(Subj, Expt);
    summarize(Subj);
    
    active_trials = find(sum(get_mat(Subj,'regressors','conds')));
    actives_selector = zeros(1,size(get_mat(Subj,'regressors','conds'),2));
    actives_selector(active_trials) =1;
    Subj = initset_object(Subj,'selector','conditions_of_interest',actives_selector);
    
    % create selectors that will determine training and testing scheme
    fprintf('\n\nPrepping selectors for CM%03d...',iSub);
    Subj = cm_create_sl_selectors(Subj);
    
    summarize(Subj);
    
    Subj_prebalancing = Subj;
    
    for iResIteration = 1:Expt.num_results_iter
        
        Subj = Subj_prebalancing;
        new_active_trials = [];
        new_actives_selector = zeros(1,size(Subj.regressors{1}.mat,2));
        
        if Expt.equate_number_of_trials_in_cond_1_and_2
            Subj = create_balanced_xvalid_selectors_searchlight(Subj, 'conds',SL_SELECTOR_TO_USE);
        end
        
        nSelectors = length(find_group(Subj,'selector',[SL_SELECTOR_TO_USE '_bal']));
        % assert only one selector, since when we have multiple selectors
        % to use it becomes difficult to figure out how to set up selectors
        % ... although, if the active_trials is just relevant for zscoring,
        % then why don't we just use the balanced selectors
        assert(nSelectors==1);
        for iSelector = 1:nSelectors
            new_active_trials = horzcat(new_active_trials, ...
                find(ismember(Subj.selectors{end-nSelectors+iSelector}.mat,[1 2 3])));
        end
        new_active_trials = unique(new_active_trials);
        new_actives_selector(new_active_trials)=1;
        Subj = initset_object(Subj,'selector',...
            'conditions_of_interest_bal_within_runs', new_actives_selector);
        active_trials = new_active_trials;
        
        if Expt.perform_second_round_of_zscoring
            display('Performing second round of z-scoring');
            pat_ix = get_number(Subj,'pattern','epi_d_hp_z_condensed');
            Subj.patterns{pat_ix}.mat(:,active_trials) = zscore(Subj.patterns{pat_ix}.mat(:,active_trials)')';
        end
        
             
        % create spheres to use for classification. given some radius,
        % produce a matrix of nvox_in_mask x nvox_in_sphere where row i
        % contains the voxels in the sphere centered at voxel i
        Subj.adj_sphere = create_adj_list(Subj,Expt.roiName,'radius',Expt.srch_radius);
        statmap_srch_arg.adj_list = Subj.adj_sphere;
        
        
        
        Subj = JR_scramble_regressors(Subj,'conds',[SL_SELECTOR_TO_USE '_bal_1'],'conditions_of_interest_bal_within_runs','conds_scrambled');
       
        
        % run searchlight classification
        Subj = feature_select(Subj, ...
            'epi_d_hp_z_condensed', ... %data
            'conds', ... % binary regs
            SL_SELECTOR_TO_USE, ... % selector
            'statmap_funct', 'statmap_searchlight', ... %function
            'statmap_arg',statmap_srch_arg, ...
            'new_map_patname', 'epi_d_hp_z_condensed_srch', ...
            'thresh', []);
        
        % rerun with the scrambled regressors
        Subj = feature_select(Subj, ...
            'epi_d_hp_z_condensed', ... %data
            'conds_scrambled', ... % binary regs
            SL_SELECTOR_TO_USE, ... % selector
            'statmap_funct', 'statmap_searchlight', ... %function
            'statmap_arg',statmap_srch_arg, ...
            'new_map_patname', 'epi_d_hp_z_condensed_srch_scrambled', ...
            'thresh', []);
        
               
        % visualize 1st resulting searchlight pattern
        figure;
        Subj = load_spm_mask(Subj,'wholebrain',fullfile(Expt.dir, '/Masks/SEPT09_MVPA_MASK_resliced4mm.nii'));
        if ~proclus
		view_pattern_overlay(Subj,'wholebrain','epi_d_hp_z_condensed_srch_1',...
            'over_cmap','redblue','autoslice',true)
    	end
        % populate res structure to save
        if iResIteration == 1
            res = Subj.patterns{end};
            res.masks =  Subj.masks;
            res.auc = nan(Subj.patterns{end}.matsize(1), Expt.num_results_iter);
            res.auc_scram_regs = nan(Subj.patterns{end}.matsize(1), Expt.num_results_iter);
            res.parameters = Expt;
        end
        res.auc(:,iResIteration) = get_mat(Subj,'pattern','epi_d_hp_z_condensed_srch_1');
        res.auc_scram_regs(:,iResIteration) = get_mat(Subj,'pattern','epi_d_hp_z_condensed_srch_scrambled_1');
        
%         % hack together header information for the new pattern, so that we
%         % can save out a Nifti corresponding to the new info
%         assert(length(Subj.patterns) == 6) % I'm being lazy in how the header gets defined, so let's make sure to assert that I change it in the future
%         Subj.patterns{6}.header.vol = Subj.patterns{5}.header.vol{1}
%         Subj.patterns{6}.header.vol.fname = ...
%             sprintf(fullfile(Expt.dir, '/CM%03d/test_searchlight_output_%s.nii'),iSub,Expt.roiName);
%         write_to_spm(Subj,'pattern','epi_d_hp_z_condensed_srch_1');
        
        
    end
    
    saveDir = fullfile(Expt.dir, 'sl_mvpa', classification_name);
    if ~exist(saveDir,'dir')
        mkdir(saveDir);
    end
    save(fullfile(saveDir, sprintf('CM%03d.mat',iSub)),'res');
    close all;
end

% if ~exist(Expt.groupMvpaDir)
%     mkdir(Expt.groupMvpaDir);
% end
% save(fullfile(Expt.group_mvpa_dir, Expt.saveName),'res', 'Subj');


end

function name = cm_make_sl_savename(Expt)
condNamesStr = strrep(strjoin(Expt.condNames),',','V');
trTeStr = strrep(num2str(Expt.which_traintest),'  ','_');
trWStr_tr = strrep(num2str(Expt.trWeights_train),' ','_');
roiStr = regexprep(Expt.roiName,'(\w*).*','$1');
name = sprintf('srchlight_conds_%s_trTe%s_trW%s_roi%s',condNamesStr,trTeStr,trWStr_tr,roiStr);

end

function Subj = cm_condense_onsets_and_patterns(Subj, Expt)
% Condenses the loaded onsets and patterns to remove timepoints of
% non-interest
condensed_runs = [];
nExamples = [];

names = Expt.ons.names;
onsets = Expt.ons.onsets;

nScanFiles = length(Expt.scanfiles);
nRuns = nScanFiles/Expt.numTpPerRun;
nPatts = size(Subj.patterns{end}.mat,2)
nCondsTotal = size(onsets,2);
assert(nScanFiles == nPatts);

% find the
Expt.condCols = makeCondCols(Expt, names);
nCondsToClassify = length(Expt.condNames);
condOnsets = zeros(nCondsToClassify, nPatts);
i = 1;
for iCond = Expt.condCols;
    % record number of examples available for each condition
    nExamples(i) = length(onsets{iCond});
    % turn the onsets into a value that corresponds with volume acquisition
    % number, rather than time passed
    time_idx = onsets{iCond}/Expt.trSecs + 1;
    % put 1's in the condOnsets to create regressors
    condOnsets(i, time_idx) = 1;
    i = i + 1;
end
% visualize the trials associated with each condition
subplot(1,2,1);
imagesc(condOnsets');

% identify rest timepoints (defined here as timepoints where neither
% condition of interest is relevant
restTp = (sum(condOnsets,1) == 0);
% create condensed regressors, which exlude timepoints of noninterest
condensedCondRegs = condOnsets(:,~restTp);
subplot(1,2,2);
imagesc(condensedCondRegs');

all_trials = sum(condOnsets,1);
for trNum = Expt.trsToAverageOver
    data_by_tr(trNum,:,:) = Expt.trWeights(trNum)*Subj.patterns{end}.mat(:,find(all_trials)+trNum-1);
end
temporally_condensed_data = squeeze(sum(data_by_tr(Expt.trsToAverageOver,:,:),1));
clear data_by_tr;

if Expt.remove_outlier_trials
    mean_across_voxels = mean(temporally_condensed_data,1);
    z_mean_across_voxels = zscore(mean_across_voxels);
    upper_outliers = find(z_mean_across_voxels> Expt.remove_outlier_trials);
    lower_outliers = find(z_mean_across_voxels< -1 * Expt.remove_outlier_trials);
    all_outliers = union(upper_outliers,lower_outliers)
    condensedCondRegs(:,all_outliers) = 0;
    fprintf('removing outliers >%i standard deviations from the mean (trials %s)',Expt.remove_outlier_trials,all_outliers)
end

% create condensed runs be removing rest timepoints
condensed_runs = Subj.selectors{1}.mat(~restTp);

% add condensed images to Subj
Subj = initset_object(Subj, 'regressors', 'conds', condensedCondRegs, 'condnames', Expt.condNames);
Subj = duplicate_object(Subj, 'pattern', 'epi_d_hp_z','epi_d_hp_z_condensed');
Subj = set_mat(Subj, 'pattern', 'epi_d_hp_z_condensed', temporally_condensed_data,'ignore_diff_size',true);
Subj = add_history(Subj,'pattern','epi_d_hp_z_condensed','Pattern created created by JR custom code');
Subj = remove_mat(Subj,'pattern','epi_d_hp_z');
summarize(Subj);

% add condensed runs variable to Subj
Subj.selectors{1}.mat = condensed_runs;
Subj.selectors{1}.matsize = size(condensed_runs);

end


function Subj = cm_create_sl_selectors(Subj)
FIRST_CM_RUN = 5;
CM_RUNS = 5:8;
runs = get_mat(Subj,'selector','runs');
nRuns = max(runs);
nTimepoints = length(runs);

runs_xval_sl = ones(nRuns, nTimepoints);
% % set up train explicit test countermeasures selector by setting all explicit trials to be training examples (1s), 3 of 4 CM trials are set as searchlight training runs (3s), and the remaining CM trial is set as a final generalization test trial (2)
% trEx_teCm_sl = ones(4, nTimepoints);
% cm_run_ix = find(ismember(runs,CM_RUNS));
% trEx_teCm_sl(:, cm_run_ix) = 3;
% for r = FIRST_CM_RUN:nRuns
%     r_ix = r-FIRST_CM_RUN+1;
%     cur_final_test_run = find(runs == nRuns-r_ix+1);
%     trEx_teCm_sl(r_ix, cur_final_test_run) = 2;
%     cur_name = sprintf('trEx_teCm_sl_%i',r_ix);
%     Subj = initset_object(Subj, 'selector', cur_name, ...
%         trEx_teCm_sl(r_ix,:),'group_name','trEx_teCm_sl');
% end
% 
% imagesc(trEx_teCm_sl);
% set(gca, 'CLim', [1 3]);
% colorbar;

trEx_teCm_sl = ones(1,nTimepoints);
cm_run_ix = find(ismember(runs,CM_RUNS));
trEx_teCm_sl(:, cm_run_ix) = 3;
trEx_teCm_sl(:,Subj.selectors{end}.mat==0) = 0;
Subj = initset_object(Subj, 'selector', 'trEx_teCm_sl_1', ...
         trEx_teCm_sl,'group_name','trEx_teCm_sl');

% set up xval selector
for r = 1:nRuns
    cur_final_test_run = find(runs==r);
    runs_xval_sl(r, cur_final_test_run) = 2;
end

%imagesc(runs_xval_sl);
%set(gca, 'CLim', [1 3]);
%colorbar

for r = 1:nRuns
    cur_searchlight_gen_run = find(runs == nRuns-r+1);
    runs_xval_sl(r, cur_searchlight_gen_run) = 3;
    runs_xval_sl(r,Subj.selectors{end}.mat==0) = 0;
    cur_name = sprintf('runs_xval_sl_%i',r);
    Subj = initset_object(Subj, 'selector', cur_name, ...
        runs_xval_sl(r,:), ...
        'group_name', 'runs_xval_sl' );
end

%imagesc(runs_xval_sl)
%colorbar



end


