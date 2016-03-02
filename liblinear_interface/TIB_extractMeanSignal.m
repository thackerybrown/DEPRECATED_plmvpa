function [subj results] = TIB_extractMeanSignal(subj, S, varargin)
%Returns mean activity by timepoint-of-interest from current ROI
% Revision 1.1
% adds functionality to save out condition-specific activity by voxel (as
% opposed to ROI-level mean) and overall class preference of each voxel
% (max signal value)


% Revision 1.0
%   currently a very simple approach - we take the same patterns that would
%   be used as testing patterns for classification, collapse across voxels, and return these
%   values for analysis at the post-processing stage (e.g. correlation with
%   behavior, logits, etc)
% *much of this code draws from "cross_validation.m" to read in the
% patterns appropriately

%specific default arguments
defaults.storevoxlevel = 0;%1 = yes. By default we only store out mean activity for whole ROI per trial
defaults.genscrambledmeans = 0;%1 = yes. By default we only store out mean activity for whole ROI per trial

args = propval(varargin,defaults);


patin = S.classifier_pattern;%pattern name fed in - this is currently THE pattern name used for classification. This may not be what you want in some circumstances.
maskgroup = S.classifier_mask;


nIterations = 1;%hardcode to 1 - we don't want, at this time, multiple iterations functionality. %length(selnames);

[patnames ispatgroup] = find_group_single(subj,'pattern',patin);
if ~ispatgroup
    patnames = cellstr(repmat(patnames{1},[nIterations 1]));
end

[masknames ismaskgroup] = find_group_single(subj,'mask',maskgroup,'repmat_times',nIterations);
if ~ismaskgroup
    disp( sprintf('Using the %s mask each time',maskgroup) );
end


% Set the current masked pattern up
cur_patname = patnames{1};
cur_maskname = masknames{1};
masked_pats = get_masked_pattern(subj,cur_patname,cur_maskname);

%collapse across voxels to generate average level of activity for each
%timepoint
meanpats = nanmean(masked_pats);%mean(subj.patterns{1,5}.mat);
if any(isnan(masked_pats(:)))%*note, there should probably be no NaNs in your mask - if there are, this is a warning sign that your ROIs may not be where you think they are.
    warning( sprintf('There are NaNs in your masked pattern...') );
end

%filter meanpats by active test indices (i.e. we only want patterns of interest that
%would have been used for testing in classification)
actives = logical(S.idxOnsets_test_in_classifier);
meanpats = meanpats(actives);

results.meanacts = meanpats;

%store mean voxel signal values separately on condition-by-condition basis,
%and class preference of each voxel
if args.storevoxlevel == 1
    nonmeanpats = masked_pats(:,actives);
    
    meanvoxpercond = zeros(length(nonmeanpats'),S.num_conds); %initialize mean voxel value*cond matrix
    for icon = 1:S.num_conds
        a = nonmeanpats(:,S.onsets_test_in_classifier{icon}); %filter trials to one condition
        b = mean(a'); %mean per voxel across trials
        meanvoxpercond(:,icon) = b';
        
    end
    
    meanvoxcondpref = zeros(length(nonmeanpats'),1); %initialize mean voxel condition preference matrix
    for ivx = 1:length(nonmeanpats')
        [c meanvoxcondpref(ivx,1)] = max(meanvoxpercond(ivx,:));
    end
    
    if args.genscrambledmeans == 1
        %% create preference maps with scrambled regs
        scramblingnum = 100; % how many scrambling iterations to do?
        
        scramidx = []; %initialize empty idx from scrambling
        for s = 1:scramblingnum
            %oldidx = num2str(scramidx);
            newidx = num2str(sum([scramidx 1]));
            %regname = ['conds' oldidx];
            newregname = ['conds' newidx];
            [subj] =  TB_scramble_regressors(subj,'conds','runs','trainActives',newregname, 'ignore_adjacency', 0);
            
            meanvoxperrandcond{s} = zeros(length(nonmeanpats'),S.num_conds); %initialize mean voxel value*cond matrix
            for icon = 1:S.num_conds
                %a = nonmeanpats(:,S.onsets_test_in_classifier{icon}); %filter trials to one condition
                x = logical(subj.regressors{1,s+1}.mat(icon,:));
                a = nonmeanpats(:,x);
                
                b = mean(a'); %mean per voxel across trials
                meanvoxperrandcond{s}(:,icon) = b';
                
            end
            
            meanvoxrandcondpref{s} = zeros(length(nonmeanpats'),1); %initialize mean voxel condition preference matrix
            for ivx = 1:length(nonmeanpats')
                [c meanvoxrandcondpref{s}(ivx,1)] = max(meanvoxperrandcond{s}(ivx,:));
            end
            
            %update scrambling idx
            scramidx = [scramidx 1];
        end
        
    end
    %% store out results
    
    results.meanvoxpercond = meanvoxpercond; %write out mean signal value per voxel separately per cond
    results.meanvoxcondpref = meanvoxcondpref; %write out mean class preference per voxel (max signal)
    results.meanvoxrandcondpref = meanvoxrandcondpref; % write out randomization "preference" maps
end

end

