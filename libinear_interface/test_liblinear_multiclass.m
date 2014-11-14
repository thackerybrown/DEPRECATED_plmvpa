function [acts scratchpad] = test_liblinear_multiclass(testpats,testtargs,scratchpad)

% **Revision 1.1**
% 11/4/2014
%
% Written by Thackery Brown to handle binary AND multiclass classification
% Generates predictions using a trained logistic regression model
%
% [ACTS SCRATCHPAD] = ll_predict(TESTPATS,TESTTARGS,SCRATCHPAD)
%
%
% Modified from test_liblinear function originally written by Alan Gordon to fit mvpa toolbox conventions
%
% New from 1.0 - added "probconvert" if statement to 1) flexibly switch
% between getting probvals or raw decision vals and 2) handle binary
% classification accordingly
%
% Princeton MVPA License:
%=====================================================================
%
% This is part of the Princeton MVPA toolbox, released under
% the GPL. See http://www.csbmb.princeton.edu/mvpa for more
% information.
%
% The Princeton MVPA toolbox is available free and
% unsupported to those who might find it useful. We do not
% take any responsibility whatsoever for any problems that
% you have related to the use of the MVPA toolbox.
%
% ======================================================================

if scratchpad.constant
    testpats = [testpats; ones(1,cols(testpats))];
end

% output predictions goes into "ACTS"
acts = zeros(size(testtargs));

%[nConds nTimepoints] = size(testtargs);
nConds = rows(testtargs);


thisTest = sparse(double(testpats));
%
% if size(testtargs,1)==2
%     theseTestLabels = 2*(testtargs(1,:)') - 1; %must use -1 and 1 labels
% else
%     for i=1:size(testtargs,2)
%         theseTestLabels(i,1) = find(testtargs(:,i))';
%     end
% end

% re-order to account for liblinear's annoying tendency to order classes by
% the order in which they appear in the training labels (as opposed to a more sensible option like 1,
% 2, 3)

testTargsReoriented = testtargs;%(scratchpad.classOrientation,:);


for i=1:size(testTargsReoriented,2)
    theseTestLabels(i,1) = find(testTargsReoriented(:,i))';
end

if strcmp(scratchpad.classType, 'svm')
    [theseLabels acc probVals]=svmpredict(theseTestLabels, thisTest', scratchpad.model);
    acts = ([probVals -1*probVals]);
    
else
    %*_predict flag '-b 1' is supposed to result in probability values,
    %rather than logits, being returned. However, this does not appear to
    %work with these mac-compiled files.
    probconvert = 0; %output RELATIVE probability values instead of one-vs-rest logits? 1 == yes, 0 == no
    if probconvert == 1
        [theseLabels acc probVals]=ll_predict(theseTestLabels, thisTest', scratchpad.model, '-b 1');%NOTE: changed this from "predict(...)" because I needed to refer to mac-friendly compiled files Alan put together for Valerie
    else
        [theseLabels acc probVals]=ll_predict(theseTestLabels, thisTest', scratchpad.model);
    end
    acts = nan(size(probVals));
    acts = probVals;
    %tries to convert liblinear's stupid shuffling back into correct class
    %values
    % AG code - doesn't work for multiclass
    %     for i = 1:length(scratchpad.classOrientation)
    %          acts(:,i) = probVals(:,scratchpad.classOrientation(i),:);
    %     end
    
    % TIB revision
    % multiclass scenario
    if nConds > 2
        for i = 1:length(scratchpad.classOrientation)
            reor = scratchpad.classOrientation(i);
            acts(:,reor) = probVals(:,(i),:);
        end
        % binary classification scenario
    elseif nConds == 2
        acts_tmp = probVals;
        if probconvert == 0
            acts_tmp(:,2) = 0;%- NOTE, this is required for (and REQUIRES) a situation where you only have 1 logit, not the 2 probability values
        end
        for i = 1:length(scratchpad.classOrientation)
            reor = scratchpad.classOrientation(i);
            acts(:,reor) = acts_tmp(:,(i),:);
        end
    end
    
end



acts = acts';
%probVals = exp(logits) ./ (1 + exp(logits));

%  if size(probVals,2)==2
%      if scratchpad.liblinOrientation==1
%          acts = [probVals'; 1-probVals'];
%      else
%          acts = [1-probVals'; probVals'];
%      end
%  else
%
%     acts = probVals'./repmat(sum(probVals'),size(probVals,2),1);
%  end

% scratchpad.w = [scratchpad.logreg.betas, -1*scratchpad.logreg.betas];



%             % test it
%             thisTest = X_svm(testIdx,:);
%             theseTestLabels = Y_svm(testIdx,:);
%             [theseLabels,~]=predict(theseTestLabels, thisTest, model, S.predictOpts);
%
%             % put the results in the acts file.;
%             guesses(testIdx) = theseLabels;