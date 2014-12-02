function [subj results] = TIB_extractMeanSignal(subj, S)
%Returns mean activity by timepoint-of-interest from current ROI
% Revision 1.0
%   currently a very simple approach - we take the same patterns that would
%   be used for classification, collapse across voxels, and return these
%   values for analysis at the post-processing stage (e.g. correlation with
%   behavior, logits, etc)

%collapse across voxels to generate average level of activity for each
%timepoint
meanpats = mean(subj.patterns{1,5}.mat);

%filter meanpats by actives (i.e. we only want patterns of interest that
%would have been used in classification)
actives = logical(subj.selectors{1,3}.mat);
meanpats = meanpats(actives);

results.meanacts = meanpats;



end

