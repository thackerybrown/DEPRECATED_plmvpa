function [subj] = PM_mvpa_load_and_preprocess_raw_data(S)

% create runs index
trial_idx = 1;
for r = 1:length(S.runs_vector)
    runs(trial_idx:trial_idx+S.runs_vector(r)-1)=r;
    trial_idx = trial_idx+S.runs_vector(r);
end
S.runs = runs;

if strcmp(S.patternType, 'betas')
    % if we are going to use betas to classify, load the betas
    subj = init_subj(S.exp_name,S.subj_id);
    subj = load_spm_mask(subj,S.roi_name,S.roi_file);
    subj = load_analyze_pattern(subj,'betas',S.roi_name, S.img_files,'single',true);

    % make runs vector
    %if S.xval
        runs = S.idxTr.sess;
%     else %commented out by TIB, 11/13/2014
%         runs = [S.idxTr.sess max(S.idxTr.sess) + S.idxTe.sess];
%     end
    
    runs = runs(1:length(subj.patterns{1,1}.mat(1,:))); %shorten vector to length of usable betas, 10/15/2014
    

    subj = init_object(subj,'selector','runs');
    subj = set_mat(subj,'selector','runs',runs);
    
   % zscore the data - not sure if useful for beta analysis
   %      subj = zscore_runs(subj,'betas','runs');
   %      subj = remove_mat(subj,'pattern','betas');
        

    
elseif strcmp(S.patternType, 'raw')
    
    if S.denoise
        
        %initialize subj structure,
        subj = init_subj(S.exp_name,S.subj_id);
        
        %load masks
        subj = load_spm_mask(subj,S.roi_name,S.roi_file);
        subj = load_spm_mask(subj,S.noiseVoxels_name,S.noiseVoxels_file);
        subj = load_spm_mask(subj,S.sigVoxels_name,S.sigVoxels_file);
        
        %create sub-masks for signal and noise
        subj = intersect_masks(subj,S.roi_name,S.sigVoxels_name, 'new_maskname', 'sigMask');
        subj = intersect_masks(subj,S.roi_name,S.noiseVoxels_name, 'new_maskname', 'noiseMask');
        
        %load functional data
        subj = load_analyze_pattern(subj,'patsAllVox', S.roi_name, S.img_files,'single',true);
        
        % make runs vector
        subj = init_object(subj,'selector','runs');
        subj = set_mat(subj,'selector','runs',runs);
        
        % zscore the data from each run
        subj = zscore_runs(subj,'patsAllVox','runs'); % Z-score the data
        subj = remove_mat(subj,'pattern','patsAllVox');
        
        % make a mask of only voxels considered to be noise (e.g. white
        % matter voxels)
        allMaskVox = find(get_mat(subj,'mask',S.roi_name));
        noiseMask = find(get_mat(subj,'mask','noiseMask'));
        noiseWithinMask = ismember(allMaskVox,noiseMask);
        
        opt = S.denoiseOpt;
        opt.noiseWithinMask = noiseWithinMask;
        
        % denoise data, separately for train and test set
        if S.xval
            [~, denoisedData] = denoiseWrapper(S, subj, 1:length(S.runs_vector), 'patsAllVox_z', opt, [], 0);
        else
            [dnResTrain, denoisedDataTrain] = denoiseWrapper(S, subj, find(S.TrainTestOrder==1), 'patsAllVox_z', opt, 'train');
            [dnResTest, denoisedDataTest] = denoiseWrapper(S, subj, find(S.TrainTestOrder==2), 'patsAllVox_z', opt, 'test');
            denoisedData = [denoisedDataTrain denoisedDataTest];
            S.dnInputsTrain = dnResTrain.inputs;
            S.dnInputsTest = dnResTest.inputs;
            clear denoisedDataTrain denoisedDataTest;
        end
        
        %put denoised data in subj struct
        cat_DND = squeeze(cat(4,denoisedData{:}));
        subj = remove_mat(subj,'pattern','patsAllVox');
        subj = duplicate_object(subj,'pattern','patsAllVox','patsAllVox_z_dn');
        subj = set_mat(subj,'pattern','patsAllVox_z_dn',cat_DND,'ignore_diff_size',true);
        
    else % if we're not going to denoise...
        
        % load patterns
        subj = init_subj(S.exp_name,S.subj_id);
        subj = load_spm_mask(subj,S.roi_name,S.roi_file);
        subj = load_analyze_pattern(subj,'spiral',S.roi_name, S.img_files,'single',true);
        
        % set runs
        subj = init_object(subj,'selector','runs');
        subj = set_mat(subj,'selector','runs',runs);
        
        % hp filter the data
        subj = hpfilter_runs(subj,'spiral','runs',100,2);
        subj = remove_mat(subj,'pattern','spiral');
        
        % zscore the data
        subj = zscore_runs(subj,'spiral_hp','runs');
        subj = remove_mat(subj,'pattern','spiral_hp');
    end
    
    %save the workspace
    if ~exist(S.workspace_dir)
        mkdir(S.workspace_dir);
    end
    cd (S.workspace_dir);
    save (S.workspace);
    
end
end


function [results, denoiseddata] = denoiseWrapper(S, subj, runSubset, patName, opt, phase)
% function to call GLM_denoise on a subj

% make a structure tmpS with information about conditions for denoising purposes.
% These may be different than regular conditions.
tmpS = S;
if S.xval
    tmpS.condsTest = tmpS.dnCondsTest;
elseif strcmp(phase, 'train')
    tmpS.condsTrain = tmpS.dnCondsTrain;
    tmpS.condsTest = tmpS.dnCondsTrain;
elseif strcmp(phase, 'test')
    tmpS.condsTrain = tmpS.dnCondsTest;
    tmpS.condsTest = tmpS.dnCondsTest;
end
tmpS = PM_mvpa_onsets_and_images(tmpS, tmpS.idxTr, tmpS.idxTe);

% get the relevant data and reshape it to GLM_denoise conventions
data_h = get_mat(subj,'pattern',patName);
data_h = reshape(data_h, size(data_h,1), 1, 1, size(data_h,2));

thisRunSubset = unique(runSubset);
colsToRemove = [];

% for each run
for i=1:length(thisRunSubset);
    
    TRsTheseRuns = find(S.runs==thisRunSubset(i)); % TRs that belong in the specified run
    
    % initialize design
    if includeJunkReg
        thisDesign = zeros(length(TRsTheseRuns), length(tmpS.onsets)+1);
    else
        thisDesign = zeros(length(TRsTheseRuns), length(tmpS.onsets));
    end
    
    % for each condition
    for j=1:length(tmpS.onsets)
        
        onsetsThisCond = sort(round(tmpS.onsets{j}/S.TR));
        theseOnsets = intersect(onsetsThisCond,TRsTheseRuns);
        
        % onsets corresponding to the given run and given condition, with a baseline of 0
        % corresponding to the time at which the run began.
        theseOnsets_fromBeginningOfRun{j} = theseOnsets - TRsTheseRuns(1) + 1;
        
        % if there are no onsets for the given run/condition, mark the
        % corresponding column for later removal from the design matrix.
        if isempty(find(theseOnsets_fromBeginningOfRun{j}))
            colsToRemove = [colsToRemove j];
        end
        
        % put onsets in the design.
        thisDesign(theseOnsets_fromBeginningOfRun{j},j) = 1;
    end
    
    data{i} = data_h(:,:,:,TRsTheseRuns); %4d data matrix
    design{i} = thisDesign; % design matrix with time as rows and condition as cols
end

%remove cols that did not have relevant onsets.
for i=1:length(design)
    design{i}(:,colsToRemove) = [];
end

%run GLMdenoise
[results,denoiseddata] = PM_GLMdenoisedata(design,data,.5,S.TR,'assume',[],opt,[]);

end