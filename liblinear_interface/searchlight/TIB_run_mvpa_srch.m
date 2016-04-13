function [res, results]= TIB_run_mvpa_srch(subj_array, task, TRsperRun, studyName, portion)%PM_run_mvpa_general(subj_array, task, TRsperRun, saveName, portion)
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
            
            if S.extractMeanSignal
                [subj results] = TIB_extractMeanSignal(subj,S);
                
                %store the results
                subjnum = str2num(subj_array{b})%convert subject number input from string to number format
                res.subj{subjnum}.nVox(1).weights(1).activity{1} = results;
                res.subj{subjnum}.nVox(1).weights(1).S = S;
                res.subjArray{subjnum} = S.subj_id;
                
                savepath = '/Users/thackery/Work/Circmaze/';
                save (fullfile(savepath, S.saveName2), 'res');
                clear subj
            else
                
                %            [subj] =  JR_scramble_regressors(subj,'conds','randomNFold','trainActives','conds_scrambled')%here we feed in 'randomNFold' as a surrogate for "runs" because we want to randomize within each "bin" used for xvalidation. We feed in trainActives for active datapoints because this also reflects all datapoints used for training and testing <- this may change depending on study!
                %[subj] =  JR_scramble_regressors(subj,'conds','runs','trainActives','conds_scrambled')
                %% run the classification.
                %             if S.class_args.nVox>0 %if we are using data-derived feature selection (e.g. top n voxels) we feed in the mask grp name such that each x-val iteration gets its own, non-biased, masked set of data
                %                 [subj results] = cross_validation(subj,S.classifier_pattern,'conds', ...
                %                     S.classSelector, S.classifier_mask_group,S.class_args, 'perfmet_functs', S.perfmet_functs);
                %             else %if we aren't doing data-driven feature selection, we just use the user-specified mask for the data
                %                 [subj results] = cross_validation(subj,S.classifier_pattern,'conds', ...
                %                     S.classSelector, S.classifier_mask,S.class_args, 'perfmet_functs', S.perfmet_functs);
                %             end
                %
                % %                          [subj results] = cross_validation(subj,S.classifier_pattern,'conds_scrambled', ...
                % %                              S.classSelector, S.classifier_mask,S.class_args, 'perfmet_functs', S.perfmet_functs);
                %
                %
                %         %set up importance maps.
                %         if S.generate_importance_maps == 1
                %             for rif = 1:length(results.iterations);
                %                 thisScratch = results.iterations(rif).scratchpad.w(2:end,:)';%liblinear
                %                 %thisScratch = results.iterations(rif).scratchpad.weights(1:end,:)';%pLR
                %                 results_IW{rif}.iterations(1).scratchpad.net.IW{1} = thisScratch;
                %             end
                %         end
                %
                %
                %         % generate importance maps.
                %         if S.generate_importance_maps
                %             TIB_generate_importance_maps(subj, results, results_IW, S)
                %         end
                %
                %         %save results
                %         if ~(exist(S.group_mvpa_dir))
                %             mkdir(S.group_mvpa_dir);
                %         end
                %
                %         if exist([fullfile(S.group_mvpa_dir, S.saveName) '.mat'], 'file')
                %             load(fullfile(S.group_mvpa_dir, S.saveName))
                %         end
                %
                %         %store results
                %         subjnum = str2num(subj_array{b})%convert subject number input from string to number format
                %         res.subj{subjnum}.penalty(1).nVox(1).weights(1).iter{n} = results;
                %         res.subj{subjnum}.penalty(1).nVox(1).weights(1).S = S;
                %         res.subjArray{subjnum} = S.subj_id;
                % %        res.subj{subjnum}.penalty(1).nVox(1).weights(1).iter{n}.confm = multiple_iterations_confusion(results);
                %         %         res.subj{b}.penalty(1).nVox(1).weights(1).iter{n} = results;
                %         %         res.subj{b}.penalty(1).nVox(1).weights(1).S = S;
                %         %         res.subjArray = subj_array;
                %
                %
                %         savepath = '/Users/thackery/Work/Circmaze/'
                %         save (fullfile(savepath, S.saveName), 'res', '-v7.3');
                %
                %         %save (fullfile(S.group_mvpa_dir, S.saveName), 'res', '-v7.3');
                %
                %
                %         %save (fullfile(S.group_mvpa_dir, S.saveName), 'res');
                %
                %         % display time it took.
                %         time2finish = toc/60;
                %         display(['Finished ' S.subj_id ' in ' num2str(time2finish) ' minutes']);
                %
                %
                %
                %         clear subj
                %         end
                %     end
                %         a.iterations = []
                %     for i = 1:length(res.subj{1,7}.penalty.nVox.weights.iter)
                %         a.iterations = [a.iterations res.subj{1,7}.penalty.nVox.weights.iter{1,i}.iterations]
                %     end
                %     testo = multiple_iterations_confusion(a)
                
                
                
                
                
                % Runs a searchlight classification analysis on the countermeasures dataset
                % Args:
                %   subjArray: array of subject numbers to include in the analysis
                %
                % Returns:
                %   Map of highest classification performance
                
                % ID numbers for all subjects in this study
                % ALL_SUBS = [1,3:10,12:26];
                % % this selector will be set up by cm_create_sl_selectors and its
                % % name will be passed to feature_select
                % SL_SELECTOR_TO_USE = 'trEx_teCm_sl';
                %
                % % if specific subjects are not specified to classify
                % if isempty(subjArray)
                %     subjArray=ALL_SUBS;
                % end
                %
                % res = [];
                % % initialize some classifier and searchlight specific parameter variables
                % class_args.train_funct_name = 'train_pLR';
                % class_args.test_funct_name = 'test_pLR';
                % class_args.penalty=1;
                %
                % scratch.class_args = class_args;
                % scratch.perfmet_funct = 'perfmet_auc';
                % scratch.perfmet_args = struct([]);
                %
                % statmap_srch_arg.obj_funct = 'statmap_classify';
                % statmap_srch_arg.scratch = scratch;
                %
                % % iterate through each subject performing classification
                % for iSub = subjArray
                %     % load relevant mvpa and experiment parameters
                %     [Expt, classArgs, par] = CM_mvpa_params(iSub, 'ret', proclus);
                %     % create necessary strings
                %     Expt.sl_selector_to_use = SL_SELECTOR_TO_USE;
                %     subjId = sprintf(Expt.subjIdFormat,iSub);
                %     Expt.saveName = cm_make_sl_savename(Expt);
                %     Expt.impMapDirStr = fullfile(Expt.dir, '/mvpa_results/importance_maps/', Expt.saveName);
                %     Expt.roiFname = fullfile(Expt.dir,'Masks',Expt.roiName);
                %     thisMvpaDir = fullfile(Expt.dir, subjId, Expt.mvpaDirStr);
                %     Expt.subjFname = fullfile(thisMvpaDir, [subjId '_' Expt.roiName '_s8mm_wa.mat']);
                %
                %
                %     % training and testing on the same TR, whether you like it or not
                %     assert(isequal(Expt.trWeights_train, Expt.trWeights_test));
                %     Expt.trWeights = Expt.trWeights_train;
                %
                %     % load or create the subj structure and output a summary
                %     if ~exist(Expt.subjFname,'file')
                %         Subj = CM_mvpa_load_and_preprocess_raw_data(subjId, Expt, nRuns, 1);
                %     else
                %         Subj = load(Expt.subjFname,'subj');
                %         Subj=Subj.subj;
                %     end
                %     summarize(Subj);
                %
                %     % load onsets and names of scan files
                %     onsetsFile = fullfile(thisMvpaDir, Expt.onsetsFname);
                %     Expt.ons = load(onsetsFile);
                %     Expt.scanfiles = vertcat(par.swascanfiles.(par.task));
                %
                %     % condense onsets and patterns
                %     fprintf('\n\nPrepping onsets and patterns for CM%03d...',iSub);
                %     Subj = cm_condense_onsets_and_patterns(Subj, Expt);
                %     summarize(Subj);
                %
                %     active_trials = find(sum(get_mat(Subj,'regressors','conds')));
                %     actives_selector = zeros(1,size(get_mat(Subj,'regressors','conds'),2));
                %     actives_selector(active_trials) =1;
                %     Subj = initset_object(Subj,'selector','conditions_of_interest',actives_selector);
                %
                %     % create selectors that will determine training and testing scheme
                %     fprintf('\n\nPrepping selectors for CM%03d...',iSub);
                %     Subj = cm_create_sl_selectors(Subj);
                %
                %     summarize(Subj);
                %
                %     Subj_prebalancing = Subj;
                %
                %     for iResIteration = 1:Expt.num_results_iter
                %
                %         Subj = Subj_prebalancing;
                %         new_active_trials = [];
                %         new_actives_selector = zeros(1,size(Subj.regressors{1}.mat,2));
                %
                %         if Expt.equate_number_of_trials_in_cond_1_and_2
                %             Subj = create_balanced_xvalid_selectors_searchlight(Subj, 'conds',SL_SELECTOR_TO_USE);
                %         end
                %
                %         nSelectors = length(find_group(Subj,'selector',[SL_SELECTOR_TO_USE '_bal']));
                %         % assert only one selector, since when we have multiple selectors
                %         % to use it becomes difficult to figure out how to set up selectors
                %         % ... although, if the active_trials is just relevant for zscoring,
                %         % then why don't we just use the balanced selectors
                %         assert(nSelectors==1);
                %         for iSelector = 1:nSelectors
                %             new_active_trials = horzcat(new_active_trials, ...
                %                 find(ismember(Subj.selectors{end-nSelectors+iSelector}.mat,[1 2 3])));
                %         end
                %         new_active_trials = unique(new_active_trials);
                %         new_actives_selector(new_active_trials)=1;
                %         Subj = initset_object(Subj,'selector',...
                %             'conditions_of_interest_bal_within_runs', new_actives_selector);
                %         active_trials = new_active_trials;
                %
                %         if Expt.perform_second_round_of_zscoring
                %             display('Performing second round of z-scoring');
                %             pat_ix = get_number(Subj,'pattern','epi_d_hp_z_condensed');
                %             Subj.patterns{pat_ix}.mat(:,active_trials) = zscore(Subj.patterns{pat_ix}.mat(:,active_trials)')';
                %         end
                
                S.srch_radius = 4.8; %radius in mm, might make sense to use a multiple of voxel dim
                
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
                %   files AND store out logits per class per trial (VERY large files).
                %   3 = yes, move patterns from searchlight to separate
                %   files BUT remove scratch.class_scratch after it is no
                %longer immediately needed and DON'T write out logits per trial. NOTE: you may still want this
                %for some functions or debugging scenarios.
                statmap_srch_arg.memsave = 3;
                
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
                
                % run searchlight classification
                %scrambled classification analysis
                if S.scrambleregs == 1
                    if strcmp(S.xvaltype,'nf')%if run labels have been replaced by random nfolds
                        [subj] =  JR_scramble_regressors(subj,'conds','randomNFold','trainActives','conds_scrambled');%here we feed in 'randomNFold' as a surrogate for "runs" because we want to randomize within each "bin" used for xvalidation. We feed in trainActives for active datapoints because this also reflects all datapoints used for training and testing <- this may change depending on study!
                    else %if we are doing 'loo'
                        [subj] =  JR_scramble_regressors(subj,'conds','runs','trainActives','conds_scrambled');%here we feed in 'randomNFold' as a surrogate for "runs" because we want to randomize within each "bin" used for xvalidation. We feed in trainActives for active datapoints because this also reflects all datapoints used for training and testing <- this may change depending on study!
                    end
                    
                    % run with the scrambled regressors
                    subj = feature_select(subj, ...
                        S.classifier_pattern, ... %data
                        'conds_scrambled', ... % binary regs
                        S.classSelector, ... % selector
                        'statmap_funct', 'TIB_statmap_searchlight', ... %function
                        'statmap_arg',statmap_srch_arg, ...
                        'new_map_patname', 'epi_d_hp_z_condensed_srch_scrambled', ...
                        'thresh', []);
                    
                elseif S.scrambleregs == 0
                    subj = feature_select(subj, ...
                        S.classifier_pattern, ... %'epi_d_hp_z_condensed', ... %data
                        'conds', ... % binary regs
                        S.classSelector, ... %SL_SELECTOR_TO_USE, ... % selector, typically runs, or balanced run iterations (runs*balancing its)
                        'statmap_funct', 'TIB_statmap_searchlight', ... %function
                        'statmap_arg',statmap_srch_arg, ...
                        'new_map_patname', 'epi_d_hp_z_condensed_srch', ...
                        'thresh', []);
                    
                end
                
                
                %% WRITE OUT MEAN SEARCHLIGHT MAP TO .IMG FILE
                if statmap_srch_arg.memsave >= 1 %If maps for each iteration were stored separately
                    
                    % average searchlight performance maps across runs
                    for r=1:length(selnames)
                        tmp = load(subj.patterns{end-length(selnames)+r}.movehd.pathfilename);
                        
                        sl_voxel_values(:,r)=tmp.mat;
                    end
                    sl_mean_voxel_values = mean(sl_voxel_values,2);
                    
                    
                    
                    vol_info = subj.patterns{1,3}.header.vol{1,1}{1,1}
                    
                    vol_info.fname = [S.group_mvpa_dir '/' S.saveName '_3vox_radius_searchlight' num2str(n) '.img'];%critical - otherwise you will overwrite your beta!
                    
                    sl_map = zeros(vol_info.dim);
                    
                    
                    included_voxels = find(subj.masks{1,3}.mat);
                    sl_map(included_voxels) = sl_mean_voxel_values.*100; %multiply by 100 to improve scaling for visualization (i.e., proportion to percent)
                    
                    spm_write_vol(vol_info, sl_map);
                    
                    if statmap_srch_arg.memsave < 3
                        %save out subj for later map reconstruction
                        subjfpath = [S.group_mvpa_dir '/'  S.saveName 'subjfile'];
                        save (subjfpath, 'subj', '-v7.3');
                    end
                    
                else
                    
                    proclus = 0; %0 = run locally (not proclus)
                    
                    if proclus == 0
                        % visualize 1st resulting searchlight pattern
                        figure;
                        subj = load_spm_mask(subj,'wholebrain',fullfile(S.anat_dir, '/mtl_ctx/r5705_2_1.nii'));
                        
                        %subj = load_spm_mask(subj,'wholebrain',fullfile(Expt.dir, '/Masks/SEPT09_MVPA_MASK_resliced4mm.nii'));
                        
                        view_pattern_overlay(subj,'wholebrain','epi_d_hp_z_condensed_srch_1',...
                            'over_cmap','jet','autoslice',true)
                    end
                    % populate res structure to save
                    %         if iResIteration == 1
                    %             res = subj.patterns{end};
                    %             res.masks =  subj.masks;
                    %             res.auc = nan(subj.patterns{end}.matsize(1), Expt.num_results_iter);
                    %             res.auc_scram_regs = nan(subj.patterns{end}.matsize(1), Expt.num_results_iter);
                    %             res.parameters = Expt;
                    %         end
                    %         res.auc(:,iResIteration) = get_mat(subj,'pattern','epi_d_hp_z_condensed_srch_1');
                    %         res.auc_scram_regs(:,iResIteration) = get_mat(subj,'pattern','epi_d_hp_z_condensed_srch_scrambled_1');
                    
                    %         % hack together header information for the new pattern, so that we
                    %         % can save out a Nifti corresponding to the new info
                    %         assert(length(subj.patterns) == 6) % I'm being lazy in how the header gets defined, so let's make sure to assert that I change it in the future
                    %         subj.patterns{6}.header.vol = subj.patterns{5}.header.vol{1}
                    %         subj.patterns{6}.header.vol.fname = ...
                    %             sprintf(fullfile(Expt.dir, '/CM%03d/test_searchlight_output_%s.nii'),iSub,Expt.roiName);
                    %         write_to_spm(subj,'pattern','epi_d_hp_z_condensed_srch_1');
                    
                    
                    % average searchlight performance maps across runs
                    for r=1:length(selnames)
                        sl_voxel_values(:,r)=subj.patterns{end-length(selnames)+r}.mat;
                    end
                    sl_mean_voxel_values = mean(sl_voxel_values,2);
                    
                    
                    
                    vol_info = subj.patterns{1,3}.header.vol{1,1}{1,1}
                    
                    vol_info.fname = [S.group_mvpa_dir '/' S.saveName '_3vox_radius_searchlight' num2str(n) '.img'];%critical - otherwise you will overwrite your beta!
                    
                    sl_map = zeros(vol_info.dim);
                    
                    
                    included_voxels = find(subj.masks{1,3}.mat);
                    sl_map(included_voxels) = sl_mean_voxel_values.*100; %multiply by 100 to improve scaling for visualization (i.e., proportion to percent)
                    
                    spm_write_vol(vol_info, sl_map);
                    
                end
                
                   
            end
            clear subj
        end

    end
    
    
    %         saveDir = fullfile(Expt.dir, 'sl_mvpa', classification_name);
    %         if ~exist(saveDir,'dir')
    %             mkdir(saveDir);
    %         end
    %         save(fullfile(saveDir, sprintf('CM%03d.mat',iSub)),'res');
    
    %     save (fullfile(S.group_mvpa_dir, S.saveName), 'res', '-v7.3');
    %
    %     % display time it took.
    %     time2finish = toc/60;
    %     display(['Finished ' S.subj_id ' in ' num2str(time2finish) ' minutes']);
    close all;
    
end
end


% if ~exist(Expt.groupMvpaDir)
%     mkdir(Expt.groupMvpaDir);
% end
% save(fullfile(Expt.group_mvpa_dir, Expt.saveName),'res', 'subj');





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


