
function [scratchpad] = train_liblinear_multiclass(trainpats,traintargs,in_args,cv_args)
% **Revision 1.0**
% 9/9/2014
%
% Written by Thackery Brown to handle binary AND multiclass classification
%
% Modified from train_liblinear function originally written by Alan Gordon to fit mvpa toolbox conventions
%


v = (double(trainpats'));

% if size(traintargs,1)==2
%     choice = 2*(traintargs(1,:)') - 1; %must use -1 and 1 labels
% else
for i=1:size(traintargs,2)
    choice(i,1) = find(traintargs(:,i))';
end
% end

voxel_num = size(v,2);
lambda_beta = 0;



%%  set the cost parameter
choice_set = unique(choice);
%note: MUST still MAKE THIS more GENERAL FOR >2 CLASSES
if in_args.chooseOptimalPenalty %pick the optimal cost parameter in unbiased manner
    scratchpad.opt_penalty = TIB_optimal_penalty_beta(choice,choice_set,v,trainpats,traintargs,in_args,cv_args);
    opt_penalty = scratchpad.opt_penalty
else %else manually specify desired penalty
    opt_penalty = in_args.penalty;
    scratchpad.opt_penalty = opt_penalty;
end


%% classify with the selected penalty param
% penalty either provided by XXX_mvpa_params.m or established above in a non-biased fashion by cross-validating the training data

% train classifier, with leave one out cross validation
if strcmp(in_args.classType,'svm')
    trainOpts_orig = in_args.libsvm;
    trainOpts = [trainOpts_orig ' -c ' num2str(opt_penalty) ' -r ' sprintf('%f', opt_r)];
    model = svm_train(choice, sparse(v), trainOpts);
    scratchpad.classOrientation = model.Label';
else
    trainOpts_orig = in_args.libLin ;
    trainOpts = [trainOpts_orig ' -c ' num2str(opt_penalty)];
    model = ll_train(choice, sparse(v), trainOpts);%NOTE: changed this from "train(...)" because I needed to refer to mac-friendly compiled files Alan put together for Valerie
    scratchpad.classOrientation = model.Label';
    
    
    % store the weights for importance maps
    if length(model.Label) > 2 %more than 2 classes, collect together weights
        
        for i = 1:length(scratchpad.classOrientation)
            reor = scratchpad.classOrientation(i);
            scratchpad.logreg.betas = model.w';%(1:(end));
            scratchpad.w(:,reor) = scratchpad.logreg.betas(:,i,:);
        end
        
    elseif length(model.Label) == 2 %binary classification scenario
        scratchpad.logreg.betas = [model.w(end) model.w(1:(end-1))];
        scratchpad.w(:,scratchpad.classOrientation(1),:) = scratchpad.logreg.betas';
        scratchpad.w(:,scratchpad.classOrientation(2),:) = -1*scratchpad.logreg.betas';
    else
        error('It looks like you only have 1 class in your training set... [length(model.Label)]');
    end
    
end
% **Note: liblinear naturally labels the first training class as "class 1".
% This differs from mvpa toolbox, which keeps the '1' labels as
% the first choice regardless of what the first presented training trial label is.
% The choice variable is used to ensure consistency across these toolboxes.

choice_set = unique(choice);
for i = 1:length(choice_set)
    thisChoice = choice_set(i);
    findChoice = find(choice==thisChoice);
    earliestInstanceOfChoice(i) = findChoice(1);
end



% if size(unique(choice),1)==2
%     scratchpad.liblinOrientation = choice(1);
%
%     if scratchpad.liblinOrientation==1
%         scratchpad.logreg.betas(:,1) = [model.w(end) model.w(1:(end-1))];
%     else
%         scratchpad.logreg.betas(:,1) = [-model.w(end) -model.w(1:(end-1))];
%     end
% else
%     scratchpad.logreg.betas = model.w;
% end
scratchpad.classType = in_args.classType;
scratchpad.constant = in_args.constant;
scratchpad.choice = choice;
scratchpad.model = model;


