function [subj results] = TIB_extractMeanSignal(subj, S)
%Returns mean activity by timepoint-of-interest from current ROI
% Revision 1.0
%   currently a very simple approach - we take the same patterns that would
%   be used as testing patterns for classification, collapse across voxels, and return these
%   values for analysis at the post-processing stage (e.g. correlation with
%   behavior, logits, etc)
% *much of this code draws from "cross_validation.m" to read in the
% patterns appropriately

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



end

