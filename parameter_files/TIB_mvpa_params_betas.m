function [S idxTr idxTe par]= TIB_mvpa_params_betas(subj_id, task, TRsperRun)

% establish parameters for mvpa analysis
% <subj_id> - identifier for the given subject. Can be numerical or a
% string
% <task> 'goals' or 'plan'

%% establish general parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
idxTr = [];
idxTe = [];
%par = PM_Params(subj_id, task);

%Study name
S.exp_name = 'Circmaze'; %change this to flexibly redirect the script to different studies
%S.exp_name = 'marlieketest';

%Subname
par.substr = ['cm' subj_id{1}];
S.subj_id = par.substr;

%Task type
par.task = task; %input at function call. For circmaze, this is 'goals' or 'plan'

par.TR = 2; %TR (s)

%Functional image scan selectors
par.scansSelect.goals.loc = 1:1;%***if ALL FILENAMES corresponding to ALL RUNS OF INTEREST are stored in ONE cell of raw_filenames.mat (i.e., not broken up by run), set index to 1 or 1:1. Otherwise, create indexing for elements of cell array raw_filenames.mat corresponding to task of interest (i.e. if cells runs 1:4 correspond to task 1, we want to reference {1}, {2}... in raw_filenames.mat)
par.scansSelect.plan.loc = 1:1;%***if ALL FILENAMES corresponding to ALL RUNS OF INTEREST are stored in ONE cell of raw_filenames.mat (i.e., not broken up by run), set index to 1 or 1:1. Otherwise, create indexing for elements of cell array raw_filenames.mat corresponding to task of interest (i.e. if cells runs 1:4 correspond to task 1, we want to reference {1}, {2}... in raw_filenames.mat)

%input image info
S.inputformat = 'betas'%'betas' %either 'raw' for raw bold images or 'betas' for beta images. Selection here automatically changes some params below.
if strcmp(S.inputformat, 'raw')
    data_imgs_to_use = 'raw_filenames.mat'; % .mat file containing names of all functional images (must exist for each subject; can be created by running raw_filenames = cellstr(SPM.xY.P) on subject's SPM.mat file)
elseif strcmp(S.inputformat, 'betas')
    data_imgs_to_use = 'beta_filenames.mat'; % analogous to raw_filenames, but with a boolean index for which betas correspond to which conditions. Created by calling TIB_generate_beta_filenames.mat below
end

%% tasks
%set the same if you want to train and test on the same set of trials via
%cross-validation (see section below)
S.trainTask = 'goals';%%%%%%%%%'goals' is just a placeholder for now. Eventually, change to 'goals' vs 'plan' so we can be more flexible in analysis.
S.testTask = 'goals';%'plan';%
%x-validation info
S.xvaltype = 'loo'; %set to 'loo' for leave-one-out x-validation or 'nf' for nfold using the S.nFolds defined below

%%model information
if strcmp(S.inputformat, 'raw')
    S.onsets_filename = ['onsets_' S.subj_id '_allruns']
    S.onsets_filename_tr = ['onsets_' S.subj_id '_allruns']%S.onsets_filename %added for train on 1, test on another - this assumes the data are actually in the same set of files.
    S.onsets_filename_tst = ['onsets_' S.subj_id '_allruns']%['onsets_' S.subj_id '_allruns_test']%S.onsets_filename %added for train on 1, test on another - this assumes the data are actually in the same set of files.
elseif strcmp(S.inputformat, 'betas')
    S.onsets_filename = ['onsets_' S.subj_id '_allruns_cuenew_rearranged'];%_wnoplan_d2']%['cuebycue_onsets_' S.subj_id]%
    S.onsets_filename_tr = ['onsets_' S.subj_id '_allruns_cuenew_rearranged']%['onsets_' S.subj_id '_allruns_cuengoal_rearranged']%['onsets_' S.subj_id '_allruns']%S.onsets_filename %added for train on 1, test on another - this assumes the data are actually in the same set of files.
    S.onsets_filename_tst = ['onsets_' S.subj_id '_allruns_cuenew_rearranged']%['onsets_' S.subj_id '_allruns_cuengoal_rearranged']%['onsets_' S.subj_id '_allruns_test']%S.onsets_filename %added for train on 1, test on another - this assumes the data are actually in the same set of files.

    
%    S.onsets_filename = ['onsets_' S.subj_id '_allruns_cuengoal_rearranged_goald7_wnogoalhold'];%_wnoplan_d2']%['cuebycue_onsets_' S.subj_id]%
%    S.onsets_filename_tr = ['onsets_' S.subj_id '_allruns_cuengoal_rearranged_goald7_wnogoalhold']%['onsets_' S.subj_id '_allruns_cuengoal_rearranged']%['onsets_' S.subj_id '_allruns']%S.onsets_filename %added for train on 1, test on another - this assumes the data are actually in the same set of files.
%    S.onsets_filename_tst = ['onsets_' S.subj_id '_allruns_cuengoal_rearranged_goald7_wnogoalhold']%['onsets_' S.subj_id '_allruns_cuengoal_rearranged']%['onsets_' S.subj_id '_allruns_test']%S.onsets_filename %added for train on 1, test on another - this assumes the data are actually in the same set of files.

    S.betaidx_filename = [S.subj_id '_betas_idx']
    S.betaidx_filename_tr = [S.subj_id '_betas_idx_tr']
    S.betaidx_filename_te = [S.subj_id '_betas_idx_te']
    %S.onsets_filename = ['onsets_' S.subj_id '_allruns_cuenew_rearranged'];%_wnoplan_d2']%['cuebycue_onsets_' S.subj_id]%
    
end

%% directories~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%S.expt_dir = ['/mnt/wagner/thackery/' S.exp_name '/subs_preprocessed_fermi_nospikerepair/'];
S.expt_dir = ['/Volumes/group/awagner/wagner/thackery/' S.exp_name '/subs_preprocessed_fermi_nospikerepair/'];

%windows version of path...
%S.expt_dir = ['Z:/wagner/thackery/' S.exp_name '/subs_preprocessed_fermi_nospikerepair/'];


%model directory (onsets and betas in here)
if strcmp(S.inputformat, 'raw')
    %S.mvpa_dir = [S.expt_dir S.subj_id '/model_v2/'];
    S.mvpa_dir = [S.expt_dir S.subj_id '/model_allruns_includemarginals/'];  
elseif strcmp(S.inputformat, 'betas')
    %S.mvpa_dir = [S.expt_dir S.subj_id '/model_allruns_includemarginals/betaseries_cuengoal_rearranged_goald7_nogoalhold/'];%/'];%bckup/'];
    
    S.mvpa_dir = [S.expt_dir S.subj_id '/model_allruns_includemarginals/betaseries_cue_rearranged/'];%bckup/'];
    %S.mvpa_dir = [S.expt_dir S.subj_id '/model_allruns_includemarginals/betaseries_cue_rearranged/'];%bckup/'];
    
end

S.anat_dir = [S.expt_dir S.subj_id '/Masks'];%masks in here
S.importance_maps_dir=['/Users/thackery/Work/PM_mvpa_results/ImpMaps_' date];%[S.expt_dir 'PM_mvpa_results/ImpMaps_' date  ];%*********************
S.group_mvpa_dir = [S.expt_dir 'Circmaze_mvpa_files'];%results are spit out in here
par.subdir =[S.expt_dir S.subj_id];
par.funcdir =[par.subdir '/bolds/']
%if strcmp(task, 'goals')
%    S.univar_dir = [par.subdir '/' 'analysis_loc_perc_3d'];%*******************************maybe spm con images stored here to be used for featur selection etc
%elseif strcmp(task, 'plan')
%    S.univar_dir = [par.subdir '/' 'analysis_loc_mnem'];
%end
S.workspace_dir = [par.subdir '/' 'mvpa_workspace'];

%% preprocessing~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
S.preprocType = 'spm'; % 'spm' for spm preprocessing, 'knk' for kendrick preprocessing
if strcmp(S.inputformat, 'betas')%if we are running a beta analysis, take a moment to create "model" files for the beta images
    TIB_generate_beta_filenames(S);
end

%create raw_filenames on the fly (as opposed to making it by loading data
%from an SPM.mat file
if strcmp(S.inputformat, 'raw')
runfolds = dir(fullfile(par.funcdir, 'run_*'));%dir(fullfile(par.funcdir, 'localizer*'));%
for idxr = 1:length(runfolds)
    allrawfilenames{idxr,1} = dir(fullfile(par.funcdir, runfolds(idxr).name, '/urun*.nii'));%'/swa*.nii'));%
    for idxf = 1:length(allrawfilenames{idxr})
        allrawfilepaths{idxr,1}{idxf,1} = runfolds(idxr).name;
    end
end
allrawfilenames = vertcat(allrawfilenames{:});
allrawfilepaths = vertcat(allrawfilepaths{:});
for idx = 1:length(allrawfilenames)%-1%note, we are filling in the beta file names based on how many betas OF INTEREST we have (length(betaidx)). We don't care about the error reg betas for this analysis
    raw_filenames{idx,1} = [par.funcdir char(allrawfilepaths(idx)) '/' allrawfilenames(idx).name]; %create the analog to "raw_filenames.mat" - i.e. a list of all filenames including the path
end
for idx = 1:length(raw_filenames)
    if length(raw_filenames{idx,1}) == 100%80
        raw_filenames{idx,2} = str2double(raw_filenames{idx,1}(length(raw_filenames{idx,1})-9:length(raw_filenames{idx,1})-9));
    else
        raw_filenames{idx,2} = str2double(raw_filenames{idx,1}(length(raw_filenames{idx,1})-10:length(raw_filenames{idx,1})-9));
    end
end
a = sortrows(raw_filenames, 2);
raw_filenames = a(:,1);
else
    raw_filenames = [];
end

if strcmp(S.inputformat, 'betas')
load([S.mvpa_dir '/' data_imgs_to_use]); %loads predefined cell array
%called raw_filenames or beta_filenames into memory
end



%% cross-validation scheme 
%this code directs the training/testing toward either a cross validation
%across all trials of a given set (normal cross validation) or training on
%one set of data and testing on the other (for example if you had a
%localizer or encoding task and then wanted to test retrieval patterns
if strcmp(S.trainTask, S.testTask)
    S.xval = 1;
    if strcmp(S.xvaltype, 'nf')
        S.thisSelector =  'randomNFold_xval'; % nfold cross validation selector
    elseif strcmp(S.xvaltype, 'loo')
        S.thisSelector =  'leave_one_out_xval'; % leave-one-out cross validation selector
    end
else
    S.xval = 0;
    S.thisSelector = 'TrainTestOneIterGroup'; % train on one group, test on another
end

%% information specific to the training and testing of specific subsets of data
% parTr = par for the training data = dump all params relevant to a
% specific task - e.g. 'mnem' - into a variable parTr ***don't bother for
% circmaze study***

% S.onsetsTrainDir = location of training onsets

% S.condsTrain = conditions on which to train

% S.dnCondsTrain = conditions which which to denoise, if denoising is used

% S.TrainRuns = runs of data on which to train (this is NOT x-validation -
% this is selecting subsets of runs pertinent to analysis (assuming not all
% runs are. For example if the first half of the runs pertained to a
% different task from the other half of the runs)

% S.durTrain = duration, in SECONDS, of runs used for training - this is
% the total number of usable trs from the runs of interest, multiplied by
% TR time to convert into seconds (e.g. for [165 156 153] + TR =2,
% S.durTrain = 948s)

% S.filenames_train = names of data images to use for training

% idxTr = behavioral indices for training task, used by PM_run_MVPA_general
    
if strcmp(S.trainTask,'goals')
    %parTr = PM_Params(subj_id, 'mnem');
    S.onsetsTrainDir = [S.mvpa_dir];%directory containing onsets.mat or betas_idx.mat file to be loaded in
    
    %S.condsTrain = {{'goal1'}  {'goal2'} {'goal3'} {'goal4'} {'goal5'}} ;%corresponds to the names in the onsets.mat or betas_idx.mat files. This is used to select what is being compared with what.
    S.condsTrain = {{'cues1'}  {'cues2'} {'cues3'} {'cues4'} {'cues5'}} %{{'cues1_1'}  {'cues1_2'} {'cues2_1'} {'cues2_2'} {'cues3_1'} {'cues3_2'}  {'cues4_1'} {'cues4_2'} {'cues5_1'} {'cues5_2'}}%{{'cues1'}  {'cues2'} {'cues3'} {'cues4'} {'cues5'}} %{ {'plan1'} {'plan2'} {'plan3'} {'plan4'} {'plan5'}} ;%%corresponds to the names in the onsets.mat or betas_idx.mat files. This is used to select what is being compared with what.
    %S.condsTrain = {{'cues1_1'} {'cues2_1'} {'cues3_1'} {'cues4_1'} {'cues5_1'}}%
    
    %S.condsTrain = {{'objects'}  {'scenes'} } ;%corresponds to the names in the onsets.mat or betas_idx.mat files. This is used to select what is being compared with what.

    %S.condsTrain = {{'goal_1'}  {'goal_2'} {'goal_3'}} ;%corresponds to the names in the onsets.mat or betas_idx.mat files. This is used to select what is being compared with what.
    %S.dnCondsTrain = {{'face'}  {'house'} {'noise'}};
    S.TrainRuns = par.scansSelect.(par.task).loc;%pull up indexing from PM_Params for RUNS corresponding to task of interest (i.e. if runs 2,4,6 correspond to task 1)
    %S.filenames_train = vertcat(par.swrascanfilesByRun{S.TrainRuns});%create a filename list with all image names relevant to analysis of interest (i.e., skip of runs 1,3,5 if they correspond to task 2)
    if strcmp(S.inputformat, 'raw')
        S.filenames_train = raw_filenames;%
    elseif strcmp(S.inputformat, 'betas')
        S.filenames_train = beta_filenames;%
    end
    
    %S.durTrain = sum(par.(par.task).numvols(S.TrainRuns)) * par.TR; %AG
    %version
    S.durTrain = numel(S.filenames_train) * par.TR;
    %[~, idxTr] = fMRIBehAnalysis_Loc(par);
    
elseif strcmp(S.trainTask,'percloc_2sess')
    %par2 = PM_Params3(subj_id, 'mnem');
    S.onsetsTrainDir = [S.expt_dir S.subj_id '/analysis_mvpa_Loc_3d_2sess/'];
    S.condsTrain = {{'noise'}  {'houseCor'}} ;
    S.TrainRuns1 = par.scansSelect.(par.task).loc;
    S.TrainRuns2 = par2.scansSelect.(par2.task).loc;
    S.TrainRuns = [S.TrainRuns1 S.TrainRuns2 + max(par.scansSelect.perc.all)];
    S.durTrain = sum([par.(par.task).numvols(S.TrainRuns1) par2.(par2.task).numvols(S.TrainRuns2)])* par.TR;
    S.filenames_train_h{1} = char(par.rascanfilesByRun{S.TrainRuns1});
    S.filenames_train_h{2} = char(par2.rascanfilesByRun{S.TrainRuns2});
    S.filenames_train = char(S.filenames_train_h);
end

%% testing
if strcmp(S.testTask,'goals')
    S.onsetsTestDir =[S.mvpa_dir];%directory containing onsets.mat or betas_idx.mat file to be loaded in
    %S.condsTest = {{'goal_1_plan'}  {'goal_2_plan'} {'goal_3_plan'} {'goal_4_plan'} {'goal_5_plan'}};
    S.condsTest = {{'cues1'} {'cues2'} {'cues3'} {'cues4'} {'cues5'}} %{{'cues1_1'}  {'cues1_2'} {'cues2_1'} {'cues2_2'} {'cues3_1'} {'cues3_2'}  {'cues4_1'} {'cues4_2'} {'cues5_1'} {'cues5_2'}}%{ {'plan1'} {'plan2'} {'plan3'} {'plan4'} {'plan5'}}%%{'goalhold1'} {'goalhold2'} {'goalhold3'} {'goalhold4'} {'goalhold5'}};
    %S.condsTrain = {{'cues1_1'} {'cues2_1'} {'cues3_1'} {'cues4_1'} {'cues5_1'}}%
    
    
    %S.condsTest = {{'goal1'} {'goal2'} {'goal3'} {'goal4'} {'goal5'}};
    %S.condsTest = {{'objects'}  {'scenes'} } ;
    
    %S.condsTest = {{'goal_1'}  {'goal_2'}  {'goal_3'} };
    S.nwayclass = num2str(numel(S.condsTest));%stores the number classification dimensions just for reference (i.e. is this a 5-way or a 2-way/binary classification?)
    S.TestRuns = par.scansSelect.(par.task).loc;
    %S.durTest = sum(par.(par.task).numvols(S.TestRuns)) * par.TR;
    %S.filenames_test = vertcat(par.betasByRun{S.TestRuns});
    if strcmp(S.inputformat, 'raw')
        S.filenames_test = raw_filenames;%
    elseif strcmp(S.inputformat, 'betas')
        S.filenames_test = beta_filenames;%
    end
    S.durTest = numel(S.filenames_test) * par.TR;
    %[~, idxTe] = fMRIBehAnalysis_Loc(par);
    
elseif strcmp(S.testTask,'plan')
    S.onsetsTestDir =[S.mvpa_dir];%directory containing onsets.mat or betas_idx.mat file to be loaded in
    %S.condsTest = {{'goal_1_plan'}  {'goal_2_plan'} {'goal_3_plan'} {'goal_4_plan'} {'goal_5_plan'}};
    %S.condsTest = {{'goal1'} {'goal2'} {'goal3'} {'goal4'} {'goal5'}};
    S.condsTest = {{'cues1'} {'cues2'} {'cues3'} {'cues4'} {'cues5'}} %
    %S.condsTest = {{'cues1_2'} {'cues2_2'} {'cues3_2'}  {'cues4_2'} {'cues5_2'}}%
    %S.condsTest = {{'objects'}  {'scenes'} } ;
    
    %S.condsTest = {{'goal_1'}  {'goal_2'}  {'goal_3'} };
    S.nwayclass = num2str(numel(S.condsTest));%stores the number classification dimensions just for reference (i.e. is this a 5-way or a 2-way/binary classification?)
    S.TestRuns = par.scansSelect.(par.task).loc;
    %S.durTest = sum(par.(par.task).numvols(S.TestRuns)) * par.TR;
    %S.filenames_test = vertcat(par.betasByRun{S.TestRuns});
    if strcmp(S.inputformat, 'raw')
        S.filenames_test = raw_filenames;%
    elseif strcmp(S.inputformat, 'betas')
        S.filenames_test = beta_filenames;%
    end
    S.durTest = numel(S.filenames_test) * par.TR;
    
elseif strcmp(S.testTask,'mnemonicloc_2sess') 
    S.onsetsTestDir = [S.expt_dir S.subj_id '/analysis_mvpa_Loc_3d_2sess/'];
    S.condsTest = {{'faceCor'}  {'houseCor'}} ;
    S.TestRuns1 = par.scansSelect.(par.task).loc;
    S.TestRuns2 = par2.scansSelect.(par2.task).loc;
    S.TestRuns = [S.TrainRuns1 S.TrainRuns2 + max(par.scansSelect.mnem.all)];
    S.durTest = sum([par.(par.task).numvols(S.TrainRuns1) par2.(par2.task).numvols(S.TrainRuns2)])* par.TR;  
end

S.condnames = S.condsTrain;
S.regName = 'conds';


%% Smoothing Parameters
S.funcType = 3;
if strcmp(S.inputformat, 'betas')%override hard code here if we're playing with betas
    S.funcType = 4;
end

S.smoothTxt = { 'unsmoothed' 'smoothed' 'native' 'betas'};
switch S.funcType
    case 1
    par.filesForPatterns = par.wrascanfiles.all;
    case 2
    par.filesForPatterns = par.swrascanfiles.all;
    case 3
    %par.filesForPatterns = par.rascanfiles.all;  
    par.filesForPatterns = raw_filenames;
    case 4
    par.filesForPatterns = beta_filenames;    
end

%% specify which files to load for classification
if S.xval%if we are doing nfold xval (set above via a 1)
    S.filenames = S.filenames_train;
else
    S.filenames_h{1} = S.filenames_train;
    S.filenames_h{2} = S.filenames_test;
    S.filenames = S.filenames_h{1}%char(S.filenames_h);
end
S.img_files =  mat2cell(S.filenames, [ones(1,size(S.filenames,1))], [size(S.filenames,2)]);

%% Runs Parameters

S.runs_vector = TRsperRun; % number of TRs per scanning run (coded for vector of values)

% if ismember(S.trainTask,{'percloc_2sess' 'mnemonicloc_2sess'})  
%     if S.xval
%         S.runs_vector = [par.(par.task).numvols(S.TrainRuns1) par2.(par2.task).numvols(S.TrainRuns2)];
%     else
%         S.runs_vector = [par.(par.task).numvols(S.TrainRuns1) parTr.(parTr.task).numvols(S.TrainRuns2) par.(par.task).numvols(S.TestRuns)];
%     end
% else
%     if S.xval
%         S.runs_vector =  [par.(par.task).numvols(S.TrainRuns)];
%     else
%         S.runs_vector =  [parTr.(parTr.task).numvols(S.TrainRuns) par.(par.task).numvols(S.TestRuns)];
%         S.TrainTestOrder = [ones(size(S.TrainRuns)) 2*ones(size(S.TestRuns))];
%     end
% end

S.meta_runs = S.runs_vector;
S.num_runs = length(S.runs_vector);
S.num_vols = sum(S.runs_vector);
S.TR = par.TR;

%create an index of from which run each beta comes
if strcmp(S.inputformat, 'betas')
    onsetsmat = [S.mvpa_dir S.onsets_filename];
    load(onsetsmat);
    idxTr.sess = onsets(1:length(onsets)-1);%NOTE: Specific to circmaze, we are skipping the last onsets entry because this is an array of "arrows period" onsets.
    
    varswecareabout = length(idxTr.sess);
    onsetscount = 0;
    for r = 1:length(S.runs_vector)
        onsetscount = onsetscount + (S.runs_vector(r)*2);%sums the time (s) from the current run with previously iterated through runs to set an increasing threshold
        %s.(sprintf('x%d', r)) = 
        %x1= find(S.idxTr.sess)% < 500)
        runsc.(sprintf('x%d', r)) = find([idxTr.sess{:,1:varswecareabout}] <= onsetscount);

    end
    for r = 1:length(S.runs_vector)%creates bidx(r) which has beta numbers corresponding to run 1, beta numbers corresponding to run 2, etc
        if r == 1
            runs.bidx1 = runsc.x1;
            idxTr.sess(runs.bidx1) = {1};
        else
            runs.(sprintf('bidx%d', r)) = setdiff(runsc.(sprintf('x%d', r)), runsc.(sprintf('x%d', r-1)));
            idxTr.sess(runs.(sprintf('bidx%d', r))) = {r};
        end
    end
    idxTr.sess = cell2mat(idxTr.sess)%convert to matrix format
        
end


%% Volume Parameters
S.vol_info = spm_vol(fullfile(par.funcdir, 'run_03', 'urun3_0006.nii')); %get functional data resolution info for spm .img writing

%S.roiWithNonTaksVoxels = fullfile(par.anatdir, 'tnativeOccTemp.nii');
%S.roiWithNonTaksVoxelsName = 'tnativeOccTemp.nii';

S.roi_name = 'rbilat_hipp_trace.nii';
%S.roi_name = 'rphc_trace.nii';
%S.roi_name = 'rprc_trace.nii';
%S.roi_name = 'rec_trace.nii';

%S.roi_name = 'rcd_trace.nii';
%S.roi_name = 'rcerebell.nii';
%S.roi_name = 'rv1.nii'

%S.roi_name = 'rpfc.nii';%was chance for cm25 w/ nf160
%S.roi_name = 'rvs_t.nii';
%S.roi_name = 'gm_cm12.nii';

S.roi_file = [S.expt_dir S.subj_id '/Masks/' S.roi_name]; %this is the large-scale ROI (could be wholebrain) that workspace info is calculated for. Saves time to have this volume include any sub-volumes you are interested in (e.g. MTL if you plan on looking in hippo and phc separately)

%S.roi_name = 'mask.img';
%S.roi_file = [S.mvpa_dir 'mask.img'];

%S.noiseVoxels_file = fullfile(par.anatdir, 'rnativec2V001.nii');
%S.noiseVoxels_name = 'rnativec2V001.nii';

%S.sigVoxels_file = fullfile(par.anatdir, 'tnativeOccTempGrey.nii');
%S.sigVoxels_name = 'tnativeOccTempGrey.nii';

%mask the primary data loaded in the workspace. [] = no secondary mask.
if strcmp(S.inputformat, 'raw')
    S.secondaryMask = [S.expt_dir S.subj_id '/Masks/' S.roi_name]; % secondary mask (the specific classification mask - e.g. hippocampus within MTL)
elseif strcmp(S.inputformat, 'betas')
    S.secondaryMask = [S.mvpa_dir 'mask.img'];
end


%S.secondaryMask = fullfile(par.anatdir, 'tnativeOccTempGrey.nii');
%S.secondaryMask = ['/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/masks/OTnoHipp.img'];
%S.secondaryMask = ['/Users/gordonam/Studies/AG1/mvpa/OT_2012/rhipp.img'];
%S.secondaryMask = ['/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/masks/rparietal.img'];


%% Workspace Parameters
S.use_premade_workspace = 1;
S.workspace = fullfile(S.workspace_dir, [S.subj_id '_' S.roi_name '_' S.smoothTxt{S.funcType} '_train_' S.trainTask '_test_' S.testTask S.preprocType '.mat']);

%% Pattern names
S.patternType = S.inputformat; %'raw' or 'betas'
% S.preprocPatName = 'spiral_dn';
% S.preprocPatName = 'patsAllVox_z_dn';
if strcmp(S.inputformat, 'raw')
    S.preprocPatName = 'spiral_hp_z';%stands for 'spiral imaging'_'high-pass filtered'_'z-scored'
elseif strcmp(S.inputformat, 'betas')
    S.preprocPatName = 'betas';%'betas_z';%use betas_z if z-scoring betas
end

S.preprocPatCondensedName = [S.preprocPatName '_condensed'];

if isempty(S.secondaryMask)
    S.preprocPatNameFinalMask = S.preprocPatName;
else
    S.preprocPatNameFinalMask = [S.preprocPatName '_masked'];
end

%% Artifacts
%S.artFile = (fullfile(par.artrepdir, ['art_global_modified_' par.substr])); %directory where artifact information is contained
S.inactivateArtifacts = 0; %remove artifact trials? 0 = no.

%% Iteration Parameters
S.num_results_iter = 10; % number of times to run the entire classification process (select subset of the data and train/test classifier)
S.num_iter_with_same_data = 1; % number of times to run the classfication step for a given subset of data

%% Balancing Parameters
S.equate_number_of_trials_in_groups = 0; % equate number of trials in conditions
S.numBalancedIts = 1; % number of iterations to run, with different randomization for the balancing

%% Z-Scoring and outlier detection
S.perform_second_round_of_zscoring = 0;  % z-score data again immediately prior to classification 
S.remove_artdetect_outliers = 0; % 1 = remove trials that exhibited movement or global signal artifacts as determined by ArtDetect
S.artdetect_motion_thresh = 0; % specify ArtDetect bin for motion outliers
S.artdetect_global_signal_thresh = 0; % specify ArtDetect bin for global signal outliers
S.remove_outlier_trials = 3;  % on-the-fly outlier detection/removal; specify how many std dev from whole brain mean to exclude as outliers (0 = don't exclude any trials)

%% Importance Maps
S.generate_importance_maps = 0; %visualize classifier weights
S.generateBetaMaps = 0; %use betas, instead of importance values
S.impType = {'pos' 'neg' 'both' 'raw'}; %importance map types
S.regNames = {'locationA' 'locationB' 'locationC' 'locationD' 'locationE'};

%% Special types of analysis
S.searchlightAnalysis = 0; % run a searchlight analysis
S.linReg = 0; % run an analysis with a continuous outcome variable
S.scrambleregs = 1; % run an anlysis with the class labels scrambled on a run-by-run basis.

%% Subsample %%alan stuff only - don't use
 S.subsampleToMatch = 0; %subsample trials to match quantities across them.
 S.balanceHiAndLowConf = 0;% match N of hi and low confidence trials?
 
%% voxel interactions
S.includeVoxelInteractions = 0; %include interactions among voxels? 0 = no   
S.interactionType = 2;
S.intEst = 1;
S.intConcat = 0;
S.intPThresh = .001;
S.intReportIncrement = 100;
S.intFilePrefix = 'intVox';
S.intMaskName = 'interactionMask';
S.intPatName = 'interactions';
S.intGroupName = 'interactionsGroup';
S.intUseIntsWithUniqueInfo = 1;

%% Signal intensity analysis
S.thisSigIntenseSelector = 'randomNFold_xval'; %which selector to use for signal intensity analysis
S.zscoreIntensityVals = 1; % zscore the intensity values?

%% Denoising
S.denoise = 0; %undergo denoising?
S.denoiseOpt.denoisespec = '10001'; %which parts of the glm output do we want to save?

%% Mean Signal Extraction Params
% parameters for selecting the mean signal from a class-specific ROI for each pattern.

S.extractMeanSignal = 0; %1 - do signal extraction. 0 = don't do this. 
S.defineROIsFromANOVAFS = 0; % define ROIs using ANOVA-based feature selection, instead of pre-defining them. 
S.logreg_2Features = 0; %perform a logistic regression, using the two extracted intensity vectors
%**********************************************************************************come back to this section once up and running......%%%%
% S.ROI1PatName = [S.preprocPatCondensedName '_ROI1'];
% S.ROI1_name = [ 'occipitoTemporal_faceVsScene_500vox.img'];
% S.ROI1_file  = [par.subdir '/analysis_loc_mnem/' S.ROI1_name];
% 
% S.ROI2PatName = [S.preprocPatCondensedName '_ROI2'];
% S.ROI2_name  = ['occipitoTemporal_sceneVsFace_500vox.img'];
% S.ROI2_file   = [par.subdir '/analysis_loc_mnem/' S.ROI2_name];

%% TR Weighting
%which post-stimulus TRs should be used (and if more than one, averaged
%across) before feeding data to the classifier?  S.TR_weights_set{1} -
%training weights, S.TR_weights_set{2} = testing weights

%NOTE: for beta analyses, we don't average over multiple images because
%different images = different events
if strcmp(S.inputformat, 'raw')
    %S.TR_weights_set = {[.0072 .2168 .3781 .2742 .1237] [.0072 .2168 .3781 .2742 .1237]}; %approximates the canonical haemodynamic response
    %S.TR_weights_set = {[0 0.35 0.5 0.15 0] [0 0.35 0.5 0.15 0]};%v1
    %S.TR_weights_set = {[0 0.15 0.5 0.35 0] [0 0.15 0.5 0.35 0]};%v2
    S.TR_weights_set = {{[0 0.25 0.5 0.25] [0 0.25 0.5 0.25]}};%v3
elseif strcmp(S.inputformat, 'betas')
    S.TR_weights_set = {{[1] [1]}};%give full weighting to the 1 and only image corresponding to each event
end


%% classifier parameters
S.class_args.train_funct_name = 'train_liblinear_multiclass';%'train_pLR';%; %training function
S.class_args.test_funct_name = 'test_liblinear_multiclass';%'test_pLR';%% %testing function
S.class_args.classType = 'libLin';
S.perfmet_functs = 'perfmet_maxclass'; % performance metric
S.statmap_funct = 'statmap_anova';%'AG_statmap_anova'; % performance metric
S.nPlsCompsSet = 0; % number of pls components to include. 0 = do not use pls components.
S.nFolds = 160; % number of cross validation iterations - used for nFold as opposed to run-by-run leave-one-out

S.class_args.nVox = 0; % number of voxels to select with feature selection e.g. [1000 5000 10000]
S.class_args.fseltype = 'topn' % feature selection format: top N vox (topn) or random N vox (rand)?
S.class_args.libLin = '-q -s 0 -B 1'; %arguments for liblinear
S.class_args.libsvm = '-q -s 0 -t 2 -d 3'; % arguments for libsvm
S.class_args.constant = true; % include a constant term?
S.class_args.prefitWeights = true; 
S.class_args.chooseOptimalPenalty = 0; % 1 = yes. cycle through cost parameters in the training set, and chose the optimal one?
S.class_args.penaltyRange = [.001 .005 .01 .05 .1 .5 1 5 10 50 100 500 1000 50000]; % a vector "[]" of cost parameters to cycle through
S.class_args.radialBasisSelection = [];%[.00001 .0001 .001 .01 .1 1 10];
S.class_args.nFoldsPenaltySelection = 10; % number of cross validation folds for penalty parameter selection. 
S.class_args.penalty = 1; %uncomment if not using optimal penalty
%establishment
end