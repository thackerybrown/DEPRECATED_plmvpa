function [res, results]= TIB_run_mvpa_general(subj_array, task, TRsperRun, studyName, portion)%PM_run_mvpa_general(subj_array, task, TRsperRun, saveName, portion)
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
    
    %[S idxTr idxTe par] = PM_mvpa_params(subj_array(b), task, TRsperRun);
    [S idxTr idxTe par] = TIB_mvpa_params_betas(subj_array(b), task, TRsperRun);
    
    %[S idxTr idxTe par] = TIB_localizer_mvpa_params(subj_array(b), task, TRsperRun);%runs with melina's localizer data
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
    %S.saveName = [studyName '_' S.nwayclass 'way_' S.xvaltype S.subj_id]%set name for the .mat results and data log file. Will contain all the goodies for analysis.
    
    S.subj_array = subj_array; %subjects, input to function at the "call". put in as strings of subject numbers - e.g. '12'.
    
    %% information about which TRs to include in classification
    %which weighted combination of post-stimulus TRs should be used to train the classifier?
    S.TR_weights_train = S.TR_weights_set{1}; % should sum to 1
    S.TRs_to_average_over_train = 1:length(S.TR_weights_train);
    
    %which weighted combination of post-stimulus TRs should be used to test the classifier?
    S.TR_weights_test = S.TR_weights_set{2}; % should sum to 1
    S.TRs_to_average_over_test = 1:length(S.TR_weights_test);
    
    S.TR_weights = S.TR_weights_set;
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
            for dt = 1:length(S.TR_weights_set{1})
                data_by_TR_train(dt,:,:) = S.TR_weights_train(dt)*subj.patterns{end}.mat(:,meta_runs_train+(dt-1));
            end
            temporally_condensed_data_train = squeeze(sum(data_by_TR_train(S.TRs_to_average_over_train,:,:),1));
            clear data_by_TR_train
            
            % load testing data
            data_by_TR_test = [];
            for dt = 1:length(S.TR_weights_set{2})
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
            subj = TIB_feature_select_top_N_vox(subj,S.preprocPatCondensedName,'conds',S.classSelector,'nVox_thresh',S.class_args.nVox, 'statmap_funct', S.statmap_funct, 'statmap_arg',statmap_arg, 'fseltype', S.class_args.fseltype);
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

            %scrambled classification analysis
            if S.scrambleregs == 1
                if strcmp(S.xvaltype,'nf')%if run labels have been replaced by random nfolds
                    [subj] =  JR_scramble_regressors(subj,'conds','randomNFold','trainActives','conds_scrambled');%here we feed in 'randomNFold' as a surrogate for "runs" because we want to randomize within each "bin" used for xvalidation. We feed in trainActives for active datapoints because this also reflects all datapoints used for training and testing <- this may change depending on study!
                else %if we are doing 'loo'
                    [subj] =  JR_scramble_regressors(subj,'conds','runs','trainActives','conds_scrambled');%here we feed in 'randomNFold' as a surrogate for "runs" because we want to randomize within each "bin" used for xvalidation. We feed in trainActives for active datapoints because this also reflects all datapoints used for training and testing <- this may change depending on study!
                end
                [subj results] = cross_validation(subj,S.classifier_pattern,'conds_scrambled', ...
                    S.classSelector, S.classifier_mask,S.class_args, 'perfmet_functs', S.perfmet_functs);
                
                %classify with correct (unscrambled) class labels
            elseif S.scrambleregs == 0
                if S.class_args.nVox>0 %if we are using data-derived feature selection (e.g. top n voxels) we feed in the mask grp name such that each x-val iteration gets its own, non-biased, masked set of data
                    [subj results] = cross_validation(subj,S.classifier_pattern,'conds', ...
                        S.classSelector, S.classifier_mask_group,S.class_args, 'perfmet_functs', S.perfmet_functs);
                else %if we aren't doing data-driven feature selection, we just use the user-specified mask for the data
                    [subj results] = cross_validation(subj,S.classifier_pattern,'conds', ...
                        S.classSelector, S.classifier_mask,S.class_args, 'perfmet_functs', S.perfmet_functs);
                end
                
            end
            
%                          [subj results] = cross_validation(subj,S.classifier_pattern,'conds_scrambled', ...
%                              S.classSelector, S.classifier_mask,S.class_args, 'perfmet_functs', S.perfmet_functs);

        
        %set up importance maps.
        if S.generate_importance_maps == 1
            for rif = 1:length(results.iterations);
                thisScratch = results.iterations(rif).scratchpad.w(2:end,:)';%liblinear
                %thisScratch = results.iterations(rif).scratchpad.weights(1:end,:)';%pLR
                results_IW{rif}.iterations(1).scratchpad.net.IW{1} = thisScratch;
            end
        end
        
        
        % generate importance maps.
        if S.generate_importance_maps
            TIB_generate_importance_maps(subj, results, results_IW, S)
        end
        
        %save results
        if ~(exist(S.group_mvpa_dir))
            mkdir(S.group_mvpa_dir);
        end
        
        if exist([fullfile(S.group_mvpa_dir, S.saveName) '.mat'], 'file')
            load(fullfile(S.group_mvpa_dir, S.saveName))
        end
        
        %store results
        subjnum = str2num(subj_array{b})%convert subject number input from string to number format
        res.subj{subjnum}.penalty(1).nVox(1).weights(1).iter{n} = results;
        res.subj{subjnum}.penalty(1).nVox(1).weights(1).S = S;
        res.subjArray{subjnum} = S.subj_id;
%        res.subj{subjnum}.penalty(1).nVox(1).weights(1).iter{n}.confm = multiple_iterations_confusion(results);
        %         res.subj{b}.penalty(1).nVox(1).weights(1).iter{n} = results;
        %         res.subj{b}.penalty(1).nVox(1).weights(1).S = S;
        %         res.subjArray = subj_array;
        
        
        savepath = '/Users/thackery/Work/Circmaze/'
        save (fullfile(savepath, S.saveName), 'res', '-v7.3');
        
        %save (fullfile(S.group_mvpa_dir, S.saveName), 'res', '-v7.3');
        
        
        %save (fullfile(S.group_mvpa_dir, S.saveName), 'res');
        
        % display time it took.
        time2finish = toc/60;
        display(['Finished ' S.subj_id ' in ' num2str(time2finish) ' minutes']);
        
        
        
        clear subj
        end
    end
    %         a.iterations = []
    %     for i = 1:length(res.subj{1,7}.penalty.nVox.weights.iter)
    %         a.iterations = [a.iterations res.subj{1,7}.penalty.nVox.weights.iter{1,i}.iterations]
    %     end
    %     testo = multiple_iterations_confusion(a)
end


















