function [subj, S] = PM_organizeBetasForClassification(subj, S)
%Modified: 11/13/2014 to implement train-on-one, test-on-another
%functionality for betas


conds_orig = [];
for i = 1:length(S.onsets(1,:))
    conds_orig = [conds_orig i*ones(size(S.onsets{i}))];
end
%conds_orig=[1*ones(size(S.onsets{1})) 2*ones(size(S.onsets{2})) 3*ones(size(S.onsets{3})) 4*ones(size(S.onsets{4})) 5*ones(size(S.onsets{5}))];%%%NOTE: this is hardcoded for Circmaze with 5 things to classify. Change to accomodate.
o=[S.onsets{:}];%pull all the onsets (beta numbers) of interest into vector
[~,ixOnsetsSort]=sort(o);
conds_h = conds_orig(ixOnsetsSort);

%% Train on one set of data, test on another (not xvalidation)
if strcmp(S.thisSelector, 'TrainTestOneIterGroup')
    TrainTestOneIter = [S.idxOnsets_train_in_classifier | S.idxOnsets_test_in_classifier];%[1*S.idxOnsets_train_in_classifier 2*S.idxOnsets_test_in_classifier];
    TrainTestOneIter = +TrainTestOneIter; % simple hack converts logical back to double
    TrainTestOneIter(logical(S.idxOnsets_test_in_classifier)) = 2;
    actives =  TrainTestOneIter>0;
    
    %conds = zeros(2,length(TrainTestOneIter));
    
    conds_h2 = zeros(size(TrainTestOneIter));
    conds_h2(actives) = conds_h;
    
    conds = zeros(length(S.condsTrain),length(actives));%Currently set to create n rows of zeros, where n is the number of conds specified for Training. This is likely standard, but there may be weird cases where you want to change this.
    for c = 1:length(S.condsTrain)
        conds(c,conds_h2==c)=1;
    end
    
    %conds(1,conds_h2==1)=1;
    %conds(2,conds_h2==2)=1;
    
    subj = init_object(subj,'selector','TrainTestOneIter');
    subj = set_mat(subj,'selector','TrainTestOneIter', TrainTestOneIter);
    subj = set_objfield(subj, 'selector', 'TrainTestOneIter', 'group_name', 'TrainTestOneIterGroup');
%% NF - random nfold (also used for leave-one-trial-out)    
elseif strcmp(S.thisSelector, 'randomNFold_xval')
    actives =  (S.idxOnsets_train_in_classifier | S.idxOnsets_test_in_classifier);
    train_actives = S.idxOnsets_train_in_classifier;
    
    randomNFold_h = ceil(shuffle(1:sum(actives))/(sum(actives)/S.nFolds));
    
    randomNFold = zeros(size(actives));
    randomNFold(actives) = randomNFold_h;
    
    conds_h2 = zeros(size(actives));
    conds_h2(actives) = conds_h;
    
    % Create conds matrix - contains n boolean rows where 1's
    % correspond to the RELATIVE beta number (that is, of the
    % subset of betas we are analyzing) associated with
    % condition n (e.g. [0 0 0 0 1 1 1 1 0 0 0...])
    conds = zeros(length(S.condsTrain),length(actives));%Currently set to create n rows of zeros, where n is the number of conds specified for Training. This is likely standard, but there may be weird cases where you want to change this.
    for c = 1:length(S.condsTrain)
        conds(c,conds_h2==c)=1;
    end
    
    %                 conds(1,conds_h2==1)=1;
    %                 conds(2,conds_h2==2)=1;
    %                 conds(3,conds_h2==3)=1;
    %                 conds(4,conds_h2==4)=1;
    %                 conds(5,conds_h2==5)=1;
    
    subj = init_object(subj,'selector','trainActives');
    subj = set_mat(subj,'selector','trainActives', train_actives);
    
    subj = init_object(subj,'selector','randomNFold');
    subj = set_mat(subj,'selector','randomNFold', randomNFold);
    subj = set_objfield(subj, 'selector', 'randomNFold', 'group_name', 'randomNFoldGroup');
    subj = TB_create_xvalid_indices(subj,'randomNFold', 'actives_selname', 'trainActives');
%% LOO - leave-one-run-out
elseif strcmp(S.thisSelector, 'leave_one_out_xval')
    %do magic using index from subj.selectors.runs (built from idx.tr.sess)
    actives =  (S.idxOnsets_train_in_classifier | S.idxOnsets_test_in_classifier);
    train_actives = S.idxOnsets_train_in_classifier;
    
    conds_h2 = zeros(size(actives));
    conds_h2(actives) = conds_h;
    
    %Create conds matrix - contains n boolean rows where 1's
    % correspond to the RELATIVE beta number (that is, of the
    % subset of betas we are analyzing) associated with
    % condition n (e.g. [0 0 0 0 1 1 1 1 0 0 0...])
    conds = zeros(length(S.condsTrain),length(actives));%Currently set to create n rows of zeros, where n is the number of conds specified for Training. This is likely standard, but there may be weird cases where you want to change this.
    for c = 1:length(S.condsTrain)
        conds(c,conds_h2==c)=1;
    end
    
    
    
    runinds = subj.selectors{1,1}.mat(1:length(actives));%pull up run numbers for each beta - we'll use this for xvalidation selector
    
    subj = init_object(subj,'selector','trainActives');
    subj = set_mat(subj,'selector','trainActives', train_actives);
    
    subj = init_object(subj,'selector','leave_one_out');
    subj = set_mat(subj,'selector','leave_one_out', runinds);
    subj = set_objfield(subj, 'selector', 'leave_one_out', 'group_name', 'leave_one_outGroup');
    subj = TB_create_xvalid_indices(subj,'leave_one_out', 'actives_selname', 'trainActives');
end




subj = init_object(subj,'regressors','conds');
subj = set_mat(subj,'regressors','conds',conds);
subj = set_objfield(subj,'regressors','conds','condnames',S.condnames);

S.classSelector = S.thisSelector;

% add condensed activation pattern (for betas this is just a
% copy of the existing pattern for housekeeping purposes/to
% keep naming consistent with that of BOLD analysis
subj = duplicate_object(subj,'pattern',S.preprocPatNameFinalMask,S.preprocPatCondensedName);
%subj = set_mat(subj,'pattern',S.preprocPatCondensedName,temporally_condensed_data,'ignore_diff_size',true);
zhist = sprintf('Pattern ''%s'' created by TIB custom code',S.preprocPatCondensedName);
subj = add_history(subj,'pattern',S.preprocPatCondensedName,zhist,true);

end

