function [res, results]= TIB_run_mvpa_repclustering(subj_array, task, TRsperRun, studyName, portion)%PM_run_mvpa_general(subj_array, task, TRsperRun, saveName, portion)
%% modification of searchlight code to quantify the clustering of representational information; i.e., to what degree do voxels with same class preference cluster together?
%% Dependencies (For TIB Circmaze Study)
%*_mvpa_params.m
%TIB_generate_betafilenames
%TB_mvpa_onsets_and_images
%PM_mvpa_load_and_preprocess_raw_data
%....

for b=(1:length(subj_array))
    
    %% load general parameter information
    tic;
    %[S idxTr idxTe par] = PM_mvpa_loc_params(subj_array(b), task, TRsperRun); %PM_mvpa_params(subj_array(b), task, TRsperRun);%
    [S idxTr idxTe par] = TIB_mvpa_params_betas(subj_array(b), task, TRsperRun);
    %[S idxTr idxTe par] = TIB_surrogate_mvpa_params(subj_array(b), task, TRsperRun);%runs with melina's localizer data
    
    if nargin > 4%3 %AG had 3 here... set to 4 just to keep things rolling for now with testing the script
        S.portion = portion;%******QUESTION for AG****what was this used for?
    else
        S.portion = [];
    end
    S.idxTr = idxTr;
    S.idxTe = idxTe;
    S.saveName = [studyName '_' S.nwayclass 'way_' S.xvaltype '_' S.subj_id]%set name for the .mat results and data log file. Will contain all the goodies for analysis.
    S.saveName2 = [studyName '_' S.nwayclass 'way_MeanActivity' S.subj_id]
    S.subj_array = subj_array; %subjects, input to function at the "call". put in as strings of subject numbers - e.g. '12'.
    
    for trset = 1:length(S.TR_weights_set) %for each set of TRs in the TR weights set (used to test multiple TR weightings)
        %% information about which TRs to include in classification
        %which weighted combination of post-stimulus TRs should be used to train the classifier?
        S.TR_weights_train = S.TR_weights_set{trset}{1}; % should sum to 1
        S.TRs_to_average_over_train = 1:length(S.TR_weights_train);
        
        %which weighted combination of post-stimulus TRs should be used to test the classifier?
        S.TR_weights_test = S.TR_weights_set{trset}{2}; % should sum to 1
        S.TRs_to_average_over_test = 1:length(S.TR_weights_test);
        
        S.TR_weights = S.TR_weights_set{trset};
        S.TRs_to_average_over = 1:length(S.TR_weights);
        
        %% Onsets
        S = TB_mvpa_onsets_and_images(S);%PM_mvpa_onsets_and_images(S);
        S.num_conds = size(S.onsets,2);
        
        
        %% Workspace stuff
        existWorkspace = exist(S.workspace);
        
        for n = 1: S.num_results_iter
            % load workspace
            if (S.use_premade_workspace&&existWorkspace) % if we are supposed to use premade workspace, and one with the correct name exists
                load(S.workspace, 'subj');
            else
                [subj] = PM_mvpa_load_and_preprocess_raw_data(S);
            end
            
            %% mask a workspace mask with another mask.
            if ~isempty(S.secondaryMask)
                subj = load_spm_mask(subj, 'secondaryMask', S.secondaryMask);
                subj = intersect_masks(subj,S.roi_name,'secondaryMask');
                subj = create_pattern_from_mask(subj, S.preprocPatName, subj.masks{end}.name , [S.preprocPatName '_masked']);
            end
            
            %% begin classifier loops
            
            if strcmp(S.patternType, 'raw')
                all_regs = zeros(S.num_conds,S.num_vols); % initialize regs matrix as conditions x timepoints
                
                % convert from seconds to TRs
                for cond = 1:S.num_conds
                    for trial = 1: length(S.onsets{cond})
                        time_idx = round(S.onsets{cond}(trial)/S.TR) + 1; % divide by 2 and add 1 to convert back from sec to TRs (first timepoint = 0 sec; first TR = 1)
                        all_regs(cond, round(time_idx)) = 1;
                    end
                end
                
                % condense regs by removing zeros
                condensed_runs = [];
                condensed_regs_of_interest = [];
                trial_counter = 1;
                for i = 1: size(all_regs,2)
                    if ~isempty(find(all_regs(:,i))) % if not a rest timepoint
                        condensed_regs_of_interest(:,trial_counter) = all_regs(:,i);
                        condensed_runs(1,trial_counter) = subj.selectors{1}.mat(i);
                        trial_counter = trial_counter + 1;
                    end
                end
                idx_condense =find(sum(all_regs));
                
                % condense meta_runs using idx_condense - create a "run" label
                % corresponding to each onset?
                trial_idx = 1;
                m_runs = 0;
                for r = 1:length(S.meta_runs)
                    m_runs(trial_idx:trial_idx+S.meta_runs(r)-1)=r;
                    trial_idx = trial_idx+S.meta_runs(r);
                end
                meta_runs_condensed = m_runs(idx_condense);
                
                
                %% select active trials
                S.actives = ones(size(meta_runs_condensed));
                
                % index active training trials
                allTrainOns = sort([S.onsets_train_in_classifier{:}]);
                allOns = sort([S.onsets{:}]);
                S.trainActives = ismember(allOns, allTrainOns);
                subj = init_object(subj,'selector','trainActives');
                subj = set_mat(subj,'selector','trainActives', S.trainActives);
                
                subj = init_object(subj,'selector','actives');
                subj = set_mat(subj,'selector','actives', S.actives);
                
                
                %% load data
                
                % create meta_runs structure
                all_trials = sum(all_regs,1);
                if strcmp(S.trainTask, S.testTask)
                    meta_runs_train = find(all_trials);
                    meta_runs_test = [];
                    if strcmp(S.xvaltype, 'loo')%leave-one-out based on run number
                        subj = init_object(subj,'selector','leave_one_out');
                        subj = set_mat(subj,'selector','leave_one_out', meta_runs_condensed);
                        subj = set_objfield(subj, 'selector', 'leave_one_out', 'group_name', 'leave_one_outGroup');
                        subj = PM_create_xvalid_indices_trainActivesOnly(subj,'leave_one_out', 'actives_selname', 'trainActives');
                    elseif strcmp(S.xvaltype, 'nf')%nfold based on param set in PM_mvpa_params
                        randomNFold = ceil(shuffle(1:length(meta_runs_condensed))/(length(meta_runs_condensed)/S.nFolds));%this is NOT leave-one-out xval. This is for random NFold xval
                        %You specify "I want to run 10 xvalidation iterations" (in PM_mvpa_params) and then we are are going to leave 1/10th of the data out during xval, randomly distributed across runs.
                        %This ideal in cases where you have few runs, so you are able to squeeze more iterations out of the data (and 10 is just kind of standard).
                        %BUT if you have lots of runs (say 16), 1) you will be doing fewer iterations than you could (unless you set the nfold to 16) and 2) you are not taking advantage of the independence afforded by testing on data from a run untouched in training.
                        subj = init_object(subj,'selector','randomNFold');
                        subj = set_mat(subj,'selector','randomNFold', randomNFold);
                        subj = set_objfield(subj, 'selector', 'randomNFold', 'group_name', 'randomNFoldGroup');
                        subj = PM_create_xvalid_indices_trainActivesOnly(subj,'randomNFold', 'actives_selname', 'trainActives');
                    end
                else
                    %applies only when train and test data are different
                    [x, trainind] = ismember(allTrainOns, allOns); %since we have doubled the number of trial instances, find indices that correspond to training trials
                    TrainTestOneIter = 2*ismember(meta_runs_condensed, 1:length(S.runs_vector));%first set all values = 2 (test)
                    TrainTestOneIter(trainind) = 1;%now replace 2s with 1s for the patterns that are training pats
                    %TrainTestOneIter = 1*ismember(meta_runs_condensed, 1:length(S.runs_vector))+1*~S.trainActives%currently assigns a 1 to every trial, then replaces it with a 2 for every test trial (defined as the inverse of trainActives logical index)
                    %TrainTestOneIter = 1*ismember(meta_runs_condensed, 1:length(S.TrainRuns)) + 2*ismember(meta_runs_condensed, length(S.TrainRuns)+1: length(S.runs_vector));% *this code only seems to work where we train on 1 (first) run and test on rest
                    TrainTestOneIter(TrainTestOneIter==1) = S.trainActives(TrainTestOneIter==1);% *this code may cause problems in some circumstances - keep your eye on the output of this line.
                    meta_runs_train = idx_condense(find(TrainTestOneIter==1));
                    meta_runs_test = idx_condense(find(TrainTestOneIter==2));
                    subj = init_object(subj,'selector','TrainTestOneIter');
                    subj = set_mat(subj,'selector','TrainTestOneIter', TrainTestOneIter);
                    subj = set_objfield(subj, 'selector', 'TrainTestOneIter', 'group_name', 'TrainTestOneIterGroup');
                end
                
                % load training data
                data_by_TR_train = [];
                for dt = 1:length(S.TR_weights_set{trset}{1})
                    data_by_TR_train(dt,:,:) = S.TR_weights_train(dt)*subj.patterns{end}.mat(:,meta_runs_train+(dt-1));
                end
                temporally_condensed_data_train = squeeze(sum(data_by_TR_train(S.TRs_to_average_over_train,:,:),1));
                clear data_by_TR_train
                
                % load testing data
                data_by_TR_test = [];
                for dt = 1:length(S.TR_weights_set{trset}{2})
                    data_by_TR_test(dt,:,:) = S.TR_weights_test(dt)*subj.patterns{end}.mat(:,meta_runs_test+(dt-1));
                end
                temporally_condensed_data_test = squeeze(sum(data_by_TR_test(S.TRs_to_average_over_test,:,:),1));%this will be blank, anyway, in xval situation, based on specification of meta_runs_test as [] above...
                clear data_by_TR_test
                
                % combine training and testing data
                temporally_condensed_data = horzcat(temporally_condensed_data_train, temporally_condensed_data_test);
                clear temporally_condensed_data_train;
                clear temporally_condensed_data_test;
                
                % Important note
                % train patterns and onsets are always first, followed
                % by test patterns and onsets.
                
                %create meta_runs_condensed selector
                subj = init_object(subj,'selector','meta_runs_condensed');
                subj = set_mat(subj,'selector','meta_runs_condensed', meta_runs_condensed);
                subj = create_xvalid_indices(subj,'meta_runs_condensed');
                
                % only include 'active' patterns in selector
                grp = find_group(subj, 'selector', S.thisSelector);
                for g = 1:length(grp)
                    this_mat = get_mat(subj,'selector',grp{g});
                    this_mat(this_mat==1) = this_mat(this_mat==1) .* S.actives(this_mat==1);
                    subj = set_mat(subj,'selector',grp{g},this_mat);
                end
                
                % add conditions
                subj = init_object(subj,'regressors','conds');
                subj = set_mat(subj,'regressors','conds',condensed_regs_of_interest);
                subj = set_objfield(subj,'regressors','conds','condnames',S.condnames);
                
                % add new condensed activation pattern
                subj = duplicate_object(subj,'pattern',S.preprocPatNameFinalMask,S.preprocPatCondensedName);
                subj = set_mat(subj,'pattern',S.preprocPatCondensedName,temporally_condensed_data,'ignore_diff_size',true);
                zhist = sprintf('Pattern ''%s'' created by AG custom code',S.preprocPatCondensedName);
                subj = add_history(subj,'pattern',S.preprocPatCondensedName,zhist,true);
                
                % clean up workspace to save RAM
                subj = remove_mat(subj,'pattern',S.preprocPatNameFinalMask);
                subj = remove_mat(subj,'pattern',S.preprocPatName);
                subj.selectors{1}.mat = condensed_runs;
                subj.selectors{1}.matsize = size(condensed_runs);
                
                S.classSelector = S.thisSelector;
                
            elseif strcmp(S.patternType, 'betas')
                [subj S] = PM_organizeBetasForClassification(subj,S);
            end
            
            %equate the training set.
            if S.equate_number_of_trials_in_groups
                subj = PM_balanceTrainPats(S, subj);%TB_balanceTrainPats(S, subj); %
                S.classSelector = [S.thisSelector 'balanced'];
            end
            
            S.classifier_pattern = S.preprocPatCondensedName; % data to use for classification.
            S.classifier_mask = subj.masks{end}.name; % mask to use for classification.
            
            %zscore the patterns prior to classification
            if S.perform_second_round_of_zscoring
                display('Performing second round of z-scoring')
                if strcmp(S.patternType, 'raw')
                    subj = zscore_runs_TIB(subj,S.preprocPatCondensedName,'runs'); % Z-score the data
                    S.classifier_pattern = [S.preprocPatCondensedName '_z']; % update the classifier data of interest
                elseif strcmp(S.patternType, 'betas')
                    subj = zscore_runs(subj,S.preprocPatName,'runs'); % Z-score the data
                    S.classifier_pattern = [S.preprocPatName '_z']; % update the classifier data of interest
                end
            end
            
            % run feature selection ANOVA: specify #of voxels (if desired)
            if S.class_args.nVox>0
                display('Performing feature selection')
                statmap_arg.use_mvpa_ver = 1;%statmap_arg = []; %%TIB - edited this so that we draw on Princeton MVPA ANOVA func instead of needing stats toolbox
                subj = TIB_feature_select_top_N_vox(subj,S.preprocPatCondensedName,'conds',S.classSelector,'nVox_thresh',S.class_args.nVox, 'statmap_funct', S.statmap_funct, 'statmap_arg',statmap_arg);
                S.classifier_mask = subj.masks{end}.name; % use group of masks created by ANOVA
                S.classifier_mask_group = subj.masks{end}.group_name;
            end
            
            %S.class_args.penalty = S.penaltyParams(pnl);
            
            %if S.extractMeanSignal
            [subj results] = TIB_extractMeanSignal(subj,S,'storevoxlevel',1,'genscrambledmeans',1);
            
            %                 %store the results
            %                 subjnum = str2num(subj_array{b})%convert subject number input from string to number format
            %                 res.subj{subjnum}.nVox(1).weights(1).activity{1} = results;
            %                 res.subj{subjnum}.nVox(1).weights(1).S = S;
            %                 res.subjArray{subjnum} = S.subj_id;
            %
            %                 savepath = '/Users/thackery/Work/Circmaze/';
            %                 save (fullfile(savepath, S.saveName2), 'res');
            %                 clear subj
            %else
            
            
            S.srch_radius = 3.2;%4.8; %radius in mm, might make sense to use a multiple of voxel dim
            vxrad = num2str(S.srch_radius/1.6);
            % create spheres to use for classification. given some radius,
            % produce a matrix of nvox_in_mask x nvox_in_sphere where row i
            % contains the voxels in the sphere centered at voxel i
            subj.adj_sphere = create_adj_list(subj,S.classifier_mask,'radius',S.srch_radius);
            statmap_srch_arg.adj_list = subj.adj_sphere;
            
            %compile params for searchlight
            scratch.class_args = S.class_args;
            scratch.perfmet_funct = S.perfmet_functs;
            scratch.perfmet_args = struct([]);
            statmap_srch_arg.obj_funct = 'statmap_classify';
            statmap_srch_arg.scratch = scratch;
            
            %employ strong memory-saving measures?
            %   1 = yes, move patterns from searchlight to separate files.
            %   2 = yes, move patterns from searchlight to separate
            %   files AND remove scratch.class_scratch after it is no
            %longer immediately needed. NOTE: you may still want this
            %for some functions or debugging scenarios.
            statmap_srch_arg.memsave = 2;
            
            statmap_srch_arg.parallel = 0; %run searchlights in parallel? 1 = yes
            
            %         subj = JR_scramble_regressors(subj,'conds',[SL_SELECTOR_TO_USE '_bal_1'],'conditions_of_interest_bal_within_runs','conds_scrambled');
            
            
            
            %set testing selector indices from 2s to 3s - this is a quirk of
            %the searchlight code
            %(https://code.google.com/p/princeton-mvpa-toolbox/wiki/TutorialSpheres)
            %where 2s are assumed to be an additional left-out dataset for
            %testing generalization
            
            selnames = find_group(subj,'selector',S.classSelector);
            for sn = 1:length(selnames)
                sel  = get_mat(subj,'selector',selnames{sn});
                xval3_idx = find(sel==2);
                sel(xval3_idx) = 3
                subj = set_mat(subj,'selector',selnames{sn},sel);
            end
            
         
            
            masked_by = get_objfield(subj,'pattern',S.preprocPatCondensedName,'masked_by');
            subj = initset_object(subj,'pattern','vox_pref',results.meanvoxcondpref, ...
                'masked_by',masked_by);
            
            
            %add randomization patterns
            for s = 1:100
                scramidx = num2str(s);
                scrampatname = ['vox_pref' scramidx];
                subj = initset_object(subj,'pattern',scrampatname,results.meanvoxrandcondpref{s}, ...
                    'masked_by',masked_by);
                
            end
            
            vox_prefsel = [1];
            subj = init_object(subj,'selector','vox_prefsel');
            subj = set_mat(subj,'selector','vox_prefsel', vox_prefsel);
            subj = set_objfield(subj, 'selector', 'vox_prefsel', 'group_name', 'vox_prefsel');
            
            %elseif S.scrambleregs == 0
            subj = feature_select(subj, ...
                'vox_pref', ... %data
                'conds', ... % binary regs
                'vox_prefsel', ... %SL_SELECTOR_TO_USE, ... % selector, typically runs, or balanced run iterations (runs*balancing its)
                'statmap_funct', 'TIB_statmap_repclustering_searchlight', ... %function
                'statmap_arg',statmap_srch_arg, ...
                'new_map_patname', 'epi_d_hp_z_condensed_srch', ...
                'thresh', []);
            tmptrueres = subj.Pweighted;
            %end
            
            %% generate patterns with scrambled condition labels
            tmpscramres = []; %initialize vector to store preference scores
            for s = 1:100
                scramidx = num2str(s);
                scrampatname = ['vox_pref' scramidx];
                newscrampatname = ['epi_d_hp_z_condensed_srch_' scramidx 'rand'];
                
                subj = feature_select(subj, ...
                    scrampatname, ... %data
                    'conds', ... % binary regs
                    'vox_prefsel', ... %SL_SELECTOR_TO_USE, ... % selector, typically runs, or balanced run iterations (runs*balancing its)
                    'statmap_funct', 'TIB_statmap_repclustering_searchlight', ... %function
                    'statmap_arg',statmap_srch_arg, ...
                    'new_map_patname', newscrampatname, ...
                    'thresh', []);
                
                tmpscramres = [tmpscramres subj.Pweighted];
                
            end
            
            
            
            %% WRITE OUT Cluster SEARCHLIGHT MAP TO .IMG FILE
         
                        
            vol_info = subj.patterns{1,3}.header.vol{1,1}{1,1}
            
            vol_info.fname = [S.group_mvpa_dir '/' S.saveName '_' vxrad 'vox_radius_clustering.img'];%critical - otherwise you will overwrite your beta!
            
            sl_map = zeros(vol_info.dim);
            
            
            included_voxels = find(subj.masks{1,3}.mat);
            sl_map(included_voxels) = subj.patterns{1,4}.mat% sloppy hardcode, change to flexibly find correct pattern in future
            
            spm_write_vol(vol_info, sl_map);
            
            %save out subj for later stats
            subjfpath = [S.group_mvpa_dir '/'  S.saveName 'clusterstats'];
            save (subjfpath, 'tmptrueres', 'tmpscramres', '-v7.3')
            
            
        end
        clear subj
    end
    
    

    close all;
    
end
end





function name = cm_make_sl_savename(Expt)
condNamesStr = strrep(strjoin(Expt.condNames),',','V');
trTeStr = strrep(num2str(Expt.which_traintest),'  ','_');
trWStr_tr = strrep(num2str(Expt.trWeights_train),' ','_');
roiStr = regexprep(Expt.roiName,'(\w*).*','$1');
name = sprintf('srchlight_conds_%s_trTe%s_trW%s_roi%s',condNamesStr,trTeStr,trWStr_tr,roiStr);

end

function subj = cm_condense_onsets_and_patterns(subj, Expt)
% Condenses the loaded onsets and patterns to remove timepoints of
% non-interest
condensed_runs = [];
nExamples = [];

names = Expt.ons.names;
onsets = Expt.ons.onsets;

nScanFiles = length(Expt.scanfiles);
nRuns = nScanFiles/Expt.numTpPerRun;
nPatts = size(subj.patterns{end}.mat,2)
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
    data_by_tr(trNum,:,:) = Expt.trWeights(trNum)*subj.patterns{end}.mat(:,find(all_trials)+trNum-1);
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
condensed_runs = subj.selectors{1}.mat(~restTp);

% add condensed images to subj
subj = initset_object(subj, 'regressors', 'conds', condensedCondRegs, 'condnames', Expt.condNames);
subj = duplicate_object(subj, 'pattern', 'epi_d_hp_z','epi_d_hp_z_condensed');
subj = set_mat(subj, 'pattern', 'epi_d_hp_z_condensed', temporally_condensed_data,'ignore_diff_size',true);
subj = add_history(subj,'pattern','epi_d_hp_z_condensed','Pattern created created by JR custom code');
subj = remove_mat(subj,'pattern','epi_d_hp_z');
summarize(subj);

% add condensed runs variable to subj
subj.selectors{1}.mat = condensed_runs;
subj.selectors{1}.matsize = size(condensed_runs);

end


function subj = cm_create_sl_selectors(subj)
FIRST_CM_RUN = 5;
CM_RUNS = 5:8;
runs = get_mat(subj,'selector','runs');
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
%     subj = initset_object(subj, 'selector', cur_name, ...
%         trEx_teCm_sl(r_ix,:),'group_name','trEx_teCm_sl');
% end
%
% imagesc(trEx_teCm_sl);
% set(gca, 'CLim', [1 3]);
% colorbar;

trEx_teCm_sl = ones(1,nTimepoints);
cm_run_ix = find(ismember(runs,CM_RUNS));
trEx_teCm_sl(:, cm_run_ix) = 3;
trEx_teCm_sl(:,subj.selectors{end}.mat==0) = 0;
subj = initset_object(subj, 'selector', 'trEx_teCm_sl_1', ...
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
    runs_xval_sl(r,subj.selectors{end}.mat==0) = 0;
    cur_name = sprintf('runs_xval_sl_%i',r);
    subj = initset_object(subj, 'selector', cur_name, ...
        runs_xval_sl(r,:), ...
        'group_name', 'runs_xval_sl' );
end

%imagesc(runs_xval_sl)
%colorbar



end


