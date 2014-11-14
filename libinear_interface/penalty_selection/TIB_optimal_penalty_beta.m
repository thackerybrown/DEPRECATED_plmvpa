
function [opt_penalty] = TIB_optimal_penalty_beta(choice,choice_set,v,trainpats,traintargs,in_args,cv_args)
% **Revision 1.1**
% 9/15/2014
%
% Performs unbiased penalty selection for pLR using liblinear
% Function written by Thackery Brown for use with liblinear classification
% Adapts code directly from A Gordon (for binary classification) and J Rissman
% (for totperf n-way classification) for use with liblinear

%MUST MAKE Multiclass scenario more GENERAL - currently only have
%hard-coded 3-way and 5-way

%currently only supports "totperf" and AUC as the performance metrics for identifying the
%best penalty.
%**"totperf" simply calculates the average % correct across the xvalidation, and selects the penalty that yielded the highest.
%**"AUC" calculates area under the curve for n classes, and then aggregates the scores using mean, media, etc.
%**Also computes F score, simply needs to have a flag set up to use this if preferred. Future revisions will include this, and will
%include class-by-class proportion correct (accuracy) as a performance metric as well


doRadialSearch = ~(strcmp(in_args.classType,'libLin') || isempty(in_args.radialBasisSelection));


%% binary classification scenario
if length(choice_set) == 2
    idx.c1 = find(choice==choice_set(1));
    idx.c2 = find(choice==choice_set(2));
    
    guesses = nan(size(choice));
    
    %loop through indices, one at a time
    randomNFold1 = ceil(shuffle(1:length(idx.c1))/(length(idx.c1)/in_args.nFoldsPenaltySelection));
    randomNFold2 = ceil(shuffle(1:length(idx.c2))/(length(idx.c2)/in_args.nFoldsPenaltySelection));
    for s = 1:in_args.nFoldsPenaltySelection %assumes balanced data
        
        thisIdx1 = randomNFold1==s;
        
        omitc1 = idx.c1(randomNFold1==s);
        omitc2 = idx.c2(randomNFold2==s);
        
        theseOmitted = [omitc1; omitc2];
        
        thisChoice = choice;
        thisChoice(theseOmitted) = [];
        
        thisV = v;
        thisV(theseOmitted,:) = [];
        
        
        for i = 1:length(in_args.penaltyRange)
            l = in_args.penaltyRange(i);
            
            if doRadialSearch
                rSet = in_args.radialBasisSelection;
            else
                rSet = 0;
            end
            
            for r = 1:length(rSet)
                thisR = rSet(r);
                if strcmp(in_args.classType,'libLin')
                    trainOpts_orig = in_args.libLin ;
                    trainOptsOptimize = [trainOpts_orig ' -c ' num2str(l)];
                    m = ll_train(thisChoice, sparse(thisV), trainOptsOptimize);%NOTE: changed this from "train(...)" because I needed to refer to mac-friendly compiled files Alan put together for Valerie
                    [theseLabels,~]=ll_predict(choice(theseOmitted), sparse(v(theseOmitted,:)), m);%NOTE: changed this from "predict(...)" because I needed to refer to mac-friendly compiled files Alan put together for Valerie
                elseif strcmp(in_args.classType,'svm')
                    trainOpts_orig = in_args.libsvm ;
                    trainOptsOptimize = [trainOpts_orig ' -c ' num2str(l) ' -r ' num2str(thisR)];
                    m = svm_train(thisChoice, thisV, trainOptsOptimize);
                    testChoice = choice(theseOmitted);
                    [theseLabels,~]=svmpredict(testChoice, v(theseOmitted,:), m);
                    
                end
                
                guesses(theseOmitted,i,r) = theseLabels;
                
            end
        end
    end
    choiceMat = repmat(choice,[1,length(in_args.penaltyRange), length(rSet)]);
    % performance measures
    perf.tp = sum(guesses == choice_set(1) & choiceMat == choice_set(1));   % true pos
    perf.fp = sum(guesses == choice_set(1) & choiceMat == choice_set(2));  % false pos
    perf.fn = sum(guesses == choice_set(2) & choiceMat == choice_set(1));   % false neg
    perf.tn = sum(guesses == choice_set(2) & choiceMat == choice_set(2));  % true neg
    
    perf.Precision = perf.tp./(perf.tp+perf.fp);
    perf.Recall = perf.tp./(perf.tp+perf.fn);
    
    perf.TrueNegRate = perf.tn./(perf.tn+perf.fp);
    %perf.Accuracy = (perf.tp+perf.tn)./ ...
    %    (perf.tp+perf.tn+perf.fp+perf.fn);
    
    perf.AUC = ((perf.tp)./(perf.tp + perf.fn) + (perf.tn)./(perf.tn + perf.fp))*.5;
    
    perf.F_Score = 2.*perf.Precision.*perf.Recall./ ...
        (perf.Precision+perf.Recall);
    
    % What performance measure to use to select the optimal
    % parameter? AUC, F score, Acc, etc
    perfParam = perf.AUC;
    
    optPerf = max(max(perfParam));
    perfParam = squeeze(perfParam);
    
    if doRadialSearch
        [idx.optPerf1 idx.optPerf2] = find(perfParam==optPerf);
        opt_penalty = in_args.penaltyRange(idx.optPerf1(1));
        opt_r = in_args.radialBasisSelection(idx.optPerf2(1));
    else
        idx.optPerf = find(perfParam==optPerf);
        opt_penalty = in_args.penaltyRange(idx.optPerf(1));
    end
    
    %% 3-way classification scenario
elseif length(choice_set) == 3
    
    idx.c1 = find(choice==choice_set(1));
    idx.c2 = find(choice==choice_set(2));
    idx.c3 = find(choice==choice_set(3));
    
    guesses = nan(size(choice));
    
    %loop through indices, one at a time
    randomNFold1 = ceil(shuffle(1:length(idx.c1))/(length(idx.c1)/in_args.nFoldsPenaltySelection));
    randomNFold2 = ceil(shuffle(1:length(idx.c2))/(length(idx.c2)/in_args.nFoldsPenaltySelection));
    randomNFold3 = ceil(shuffle(1:length(idx.c3))/(length(idx.c3)/in_args.nFoldsPenaltySelection));
    
    for s = 1:in_args.nFoldsPenaltySelection %assumes balanced data
        
        thisIdx1 = randomNFold1==s;
        
        omitc1 = idx.c1(randomNFold1==s);
        omitc2 = idx.c2(randomNFold2==s);
        omitc3 = idx.c3(randomNFold3==s);
        
        
        theseOmitted = [omitc1; omitc2; omitc3; ];
        
        thisChoice = choice;
        thisChoice(theseOmitted) = [];
        
        thisV = v;
        thisV(theseOmitted,:) = [];
        
        
        for i = 1:length(in_args.penaltyRange)
            l = in_args.penaltyRange(i);
            
            if doRadialSearch
                rSet = in_args.radialBasisSelection;
            else
                rSet = 0;
            end
            
            for r = 1:length(rSet)
                thisR = rSet(r);
                if strcmp(in_args.classType,'libLin')
                    trainOpts_orig = in_args.libLin ;
                    trainOptsOptimize = [trainOpts_orig ' -c ' num2str(l)];
                    m = ll_train(thisChoice, sparse(thisV), trainOptsOptimize);%NOTE: changed this from "train(...)" because I needed to refer to mac-friendly compiled files Alan put together for Valerie
                    [theseLabels,~]=ll_predict(choice(theseOmitted), sparse(v(theseOmitted,:)), m);%NOTE: changed this from "predict(...)" because I needed to refer to mac-friendly compiled files Alan put together for Valerie
                elseif strcmp(in_args.classType,'svm')
                    trainOpts_orig = in_args.libsvm ;
                    trainOptsOptimize = [trainOpts_orig ' -c ' num2str(l) ' -r ' num2str(thisR)];
                    m = svm_train(thisChoice, thisV, trainOptsOptimize);
                    testChoice = choice(theseOmitted);
                    [theseLabels,~]=svmpredict(testChoice, v(theseOmitted,:), m);
                    
                end
                
                guesses(theseOmitted,i,r) = theseLabels;
                
            end
        end
    end
    choiceMat = repmat(choice,[1,length(in_args.penaltyRange), length(rSet)]);
    % performance measures~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %"all the rest" vectors
    rest1 = [2 3 ];
    rest2 = [1 3 ];
    rest3 = [1 2 ];
    
    
    perf.tp1 = sum(guesses == choice_set(1) & choiceMat == choice_set(1));   % true pos
    perf.fp1 = zeros(size(perf.tp1));
    for x = 1:numel(rest1);
        perf.fp1 = perf.fp1+sum(guesses == choice_set(1) & choiceMat == choice_set(rest1(x)));
    end
    %perf.fp1 = sum(guesses == choice_set(1) & choiceMat == choice_set(2) | choice_set(3) | choice_set(4) | choice_set(5));  % false pos
    perf.fn1 = zeros(size(perf.tp1));
    for x = 1:numel(rest1);
        perf.fn1 = perf.fn1+sum(guesses == choice_set(rest1(x)) & choiceMat == choice_set(1));
    end
    %perf.fn = sum(guesses == choice_set(2) & choiceMat == choice_set(1));   % false neg
    perf.tn1 = zeros(size(perf.tp1));
    for x = 1:numel(rest1);
        perf.tn1 = perf.tn1+sum(guesses == choice_set(rest1(x)) & choiceMat == choice_set(rest1(x)));%collect true postives for the other classes
    end
    perf.tn1 = perf.tn1+sum(guesses == choice_set(3) & choiceMat == choice_set(2));
    perf.tn1 = perf.tn1+sum(guesses == choice_set(2) & choiceMat == choice_set(3));
    
    
    %perf.tn = sum(guesses == choice_set(2) & choiceMat == choice_set(2));  % true neg
    
    perf.tp2 = sum(guesses == choice_set(2) & choiceMat == choice_set(2));   % true pos
    perf.fp2 = zeros(size(perf.tp2));
    for x = 1:numel(rest2);
        perf.fp2 = perf.fp2+sum(guesses == choice_set(2) & choiceMat == choice_set(rest2(x)));
    end
    %perf.fp1 = sum(guesses == choice_set(1) & choiceMat == choice_set(2) | choice_set(3) | choice_set(4) | choice_set(5));  % false pos
    perf.fn2 = zeros(size(perf.tp2));
    for x = 1:numel(rest2);
        perf.fn2 = perf.fn2+sum(guesses == choice_set(rest2(x)) & choiceMat == choice_set(2));
    end
    %perf.fn = sum(guesses == choice_set(2) & choiceMat == choice_set(1));   % false neg
    perf.tn2 = zeros(size(perf.tp2));
    for x = 1:numel(rest2);
        perf.tn2 = perf.tn2+sum(guesses == choice_set(rest2(x)) & choiceMat == choice_set(rest2(x)));%collect true postives for the other classes
    end
    perf.tn2 = perf.tn2+sum(guesses == choice_set(3) & choiceMat == choice_set(1));
    perf.tn2 = perf.tn2+sum(guesses == choice_set(1) & choiceMat == choice_set(3));
    
    
    %perf.tn = sum(guesses == choice_set(2) & choiceMat == choice_set(2));  % true neg
    
    perf.tp3 = sum(guesses == choice_set(3) & choiceMat == choice_set(3));   % true pos
    perf.fp3 = zeros(size(perf.tp3));
    for x = 1:numel(rest3);
        perf.fp3 = perf.fp3+sum(guesses == choice_set(3) & choiceMat == choice_set(rest3(x)));
    end
    %perf.fp3 = sum(guesses == choice_set(1) & choiceMat == choice_set(2) | choice_set(3) | choice_set(4) | choice_set(5));  % false pos
    perf.fn3 = zeros(size(perf.tp3));
    for x = 1:numel(rest3);
        perf.fn3 = perf.fn3+sum(guesses == choice_set(rest3(x)) & choiceMat == choice_set(3));
    end
    %perf.fn = sum(guesses == choice_set(2) & choiceMat == choice_set(1));   % false neg
    perf.tn3 = zeros(size(perf.tp3));
    for x = 1:numel(rest3);
        perf.tn3 = perf.tn3+sum(guesses == choice_set(rest3(x)) & choiceMat == choice_set(rest3(x)));%collect true postives for the other classes
    end
    perf.tn3 = perf.tn3+sum(guesses == choice_set(2) & choiceMat == choice_set(1));
    perf.tn3 = perf.tn3+sum(guesses == choice_set(1) & choiceMat == choice_set(2));
    
    
    %perf.tn = sum(guesses == choice_set(2) & choiceMat == choice_set(2));  % true neg
    
    
    % quality metrics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    perf.Precision1 = perf.tp1./(perf.tp1+perf.fp1);
    perf.Recall1 = perf.tp1./(perf.tp1+perf.fn1);
    
    perf.TrueNegRate1 = perf.tn1./(perf.tn1+perf.fp1);%specificity
    %perf.Accuracy = (perf.tp+perf.tn)./ ...
    %    (perf.tp+perf.tn+perf.fp+perf.fn);
    
    perf.AUC1 = (perf.Recall1+perf.TrueNegRate1)*.5;%compute binary class AUC ((recall+spec)/2)
    
    perf.F_Score1 = 2.*perf.Precision1.*perf.Recall1./ ...
        (perf.Precision1+perf.Recall1);
    
    
    perf.Precision2 = perf.tp2./(perf.tp2+perf.fp2);
    perf.Recall2 = perf.tp2./(perf.tp2+perf.fn2);
    
    perf.TrueNegRate2 = perf.tn2./(perf.tn2+perf.fp2);%specificity
    %perf.Accuracy = (perf.tp+perf.tn)./ ...
    %    (perf.tp+perf.tn+perf.fp+perf.fn);
    
    perf.AUC2 = (perf.Recall2+perf.TrueNegRate2)*.5;%compute binary class AUC ((recall+spec)/2)
    
    perf.F_Score2 = 2.*perf.Precision2.*perf.Recall2./ ...
        (perf.Precision2+perf.Recall2);
    
    
    perf.Precision3 = perf.tp3./(perf.tp3+perf.fp3);
    perf.Recall3 = perf.tp3./(perf.tp3+perf.fn3);
    
    perf.TrueNegRate3 = perf.tn3./(perf.tn3+perf.fp3);%specificity
    %perf.Accuracy = (perf.tp+perf.tn)./ ...
    %    (perf.tp+perf.tn+perf.fp+perf.fn);
    
    perf.AUC3 = (perf.Recall3+perf.TrueNegRate3)*.5;%compute binary class AUC ((recall+spec)/2)
    
    perf.F_Score3 = 2.*perf.Precision3.*perf.Recall3./ ...
        (perf.Precision3+perf.Recall3);
    
    
    % use AUC or F score as the performance measure to select the optimal
    % parameter
    
    perfParam1 = perf.AUC1;
    
    optPerf1 = max(max(perfParam1));
    perfParam1 = squeeze(perfParam1);
    
    perfParam2 = perf.AUC2;
    
    optPerf2 = max(max(perfParam2));
    perfParam2 = squeeze(perfParam2);
    
    
    perfParam3 = perf.AUC3;
    
    optPerf3 = max(max(perfParam3));
    perfParam3 = squeeze(perfParam3);
    
    
    if doRadialSearch
        [idx.optPerf1 idx.optPerf2] = find(perfParam==optPerf);
        opt_penalty = in_args.penaltyRange(idx.optPerf1(1));
        opt_r = in_args.radialBasisSelection(idx.optPerf2(1));
    else
        
        stacked_perf = cat(3,perfParam1, perfParam2, perfParam3);
        
        if strcmp(in_args.multiclassPenaltySummary,'mean')
            tot_perf = mean(stacked_perf,3);
        elseif strcmp(in_args.multiclassPenaltySummary,'median')
            tot_perf = median(stacked_perf,3);
        end
        optPerfs = max(max(tot_perf));
        idx.optPerf = find(tot_perf==optPerfs);
        
        
        opt_penalty = in_args.penaltyRange(idx.optPerf(1));
    end
    
    
    %% 5-way classification scenario
elseif length(choice_set) == 5
    
    idx.c1 = find(choice==choice_set(1));
    idx.c2 = find(choice==choice_set(2));
    idx.c3 = find(choice==choice_set(3));
    idx.c4 = find(choice==choice_set(4));
    idx.c5 = find(choice==choice_set(5));
    
    guesses = nan(size(choice));
    
    %loop through indices, one at a time
    randomNFold1 = ceil(shuffle(1:length(idx.c1))/(length(idx.c1)/in_args.nFoldsPenaltySelection));
    randomNFold2 = ceil(shuffle(1:length(idx.c2))/(length(idx.c2)/in_args.nFoldsPenaltySelection));
    randomNFold3 = ceil(shuffle(1:length(idx.c3))/(length(idx.c3)/in_args.nFoldsPenaltySelection));
    randomNFold4 = ceil(shuffle(1:length(idx.c4))/(length(idx.c4)/in_args.nFoldsPenaltySelection));
    randomNFold5 = ceil(shuffle(1:length(idx.c5))/(length(idx.c5)/in_args.nFoldsPenaltySelection));
    for s = 1:in_args.nFoldsPenaltySelection %assumes balanced data
        
        thisIdx1 = randomNFold1==s;
        
        omitc1 = idx.c1(randomNFold1==s);
        omitc2 = idx.c2(randomNFold2==s);
        omitc3 = idx.c3(randomNFold3==s);
        omitc4 = idx.c4(randomNFold4==s);
        omitc5 = idx.c5(randomNFold5==s);
        
        theseOmitted = [omitc1; omitc2; omitc3; omitc4; omitc5 ];
        
        thisChoice = choice;
        thisChoice(theseOmitted) = [];
        
        thisV = v;
        thisV(theseOmitted,:) = [];
        
        
        for i = 1:length(in_args.penaltyRange)
            l = in_args.penaltyRange(i);
            
            if doRadialSearch
                rSet = in_args.radialBasisSelection;
            else
                rSet = 0;
            end
            
            for r = 1:length(rSet)
                thisR = rSet(r);
                % if we are interested in using total % correct of classifier
                if strcmp(in_args.OptimalPenaltyMetric,'totperf')
                    if strcmp(in_args.classType,'libLin')
                        trainOpts_orig = in_args.libLin ;
                        trainOptsOptimize = [trainOpts_orig ' -c ' num2str(l)];
                        m = ll_train(thisChoice, sparse(thisV), trainOptsOptimize);%NOTE: changed this from "train(...)" because I needed to refer to mac-friendly compiled files Alan put together for Valerie
                        [theseLabels acc ]=ll_predict(choice(theseOmitted), sparse(v(theseOmitted,:)), m);
                        %guesses(theseOmitted,i,r) = theseLabels;
                        accs(s,i,r) = acc;
                        
                    end
                    % if we would rather compute metrics like AUC
                else
                    if strcmp(in_args.classType,'libLin')
                        trainOpts_orig = in_args.libLin ;
                        trainOptsOptimize = [trainOpts_orig ' -c ' num2str(l)];
                        m = ll_train(thisChoice, sparse(thisV), trainOptsOptimize);%NOTE: changed this from "train(...)" because I needed to refer to mac-friendly compiled files Alan put together for Valerie
                        
                        [theseLabels,~]=ll_predict(choice(theseOmitted), sparse(v(theseOmitted,:)), m);%NOTE: changed this from "predict(...)" because I needed to refer to mac-friendly compiled files Alan put together for Valerie
                        
                    elseif strcmp(in_args.classType,'svm')
                        trainOpts_orig = in_args.libsvm ;
                        trainOptsOptimize = [trainOpts_orig ' -c ' num2str(l) ' -r ' num2str(thisR)];
                        m = svm_train(thisChoice, thisV, trainOptsOptimize);
                        testChoice = choice(theseOmitted);
                        [theseLabels,~]=svmpredict(testChoice, v(theseOmitted,:), m);
                        
                    end
                    
                    guesses(theseOmitted,i,r) = theseLabels;
                end
                
            end
        end
    end
    if strcmp(in_args.OptimalPenaltyMetric,'totperf')
        accbypen = mean(accs,1);
        optPerfs = max(max(accbypen));
        idx.optPerf = find(accbypen==optPerfs);
        %fine tune penalty?
        tune = 1;
        if tune == 1
            if idx.optPerf(1) == 1
                lower = in_args.penaltyRange(1);%current choice
                upper = (in_args.penaltyRange(2)-in_args.penaltyRange(1))/2;%halfway to next penalty
                newpenRange = lower:upper/6:upper;
                newpenRange = sort([newpenRange in_args.penaltyRange(idx.optPerf(1))]);%try 7 penalties disributed within this range
                for s = 1:in_args.nFoldsPenaltySelection %assumes balanced data
                    
                    thisIdx1 = randomNFold1==s;
                    
                    omitc1 = idx.c1(randomNFold1==s);
                    omitc2 = idx.c2(randomNFold2==s);
                    omitc3 = idx.c3(randomNFold3==s);
                    omitc4 = idx.c4(randomNFold4==s);
                    omitc5 = idx.c5(randomNFold5==s);
                    
                    theseOmitted = [omitc1; omitc2; omitc3; omitc4; omitc5 ];
                    
                    thisChoice = choice;
                    thisChoice(theseOmitted) = [];
                    
                    thisV = v;
                    thisV(theseOmitted,:) = [];
                    
                    
                    for i = 1:length(newpenRange)
                        l = newpenRange(i);
                        
                        if doRadialSearch
                            rSet = in_args.radialBasisSelection;
                        else
                            rSet = 0;
                        end
                        
                        for r = 1:length(rSet)
                            thisR = rSet(r);
                            
                            if strcmp(in_args.classType,'libLin')
                                trainOpts_orig = in_args.libLin ;
                                trainOptsOptimize = [trainOpts_orig ' -c ' num2str(l)];
                                m = ll_train(thisChoice, sparse(thisV), trainOptsOptimize);%NOTE: changed this from "train(...)" because I needed to refer to mac-friendly compiled files Alan put together for Valerie
                                [theseLabels acc]=ll_predict(choice(theseOmitted), sparse(v(theseOmitted,:)), m);
                                guesses(theseOmitted,i,r) = theseLabels;
                                accs(s,i,r) = acc;
                            end
                        end
                    end
                    
                end
                accbypen = mean(accs,1);
                optPerfs = max(max(accbypen));
                idx.optPerf = find(accbypen==optPerfs);
                opt_penalty = in_args.penaltyRange(idx.optPerf(1));
                
            elseif idx.optPerf(1) == length(in_args.penaltyRange)
                lower = (in_args.penaltyRange(idx.optPerf(1))-in_args.penaltyRange(idx.optPerf(1)-1))/2;%halfway to next penalty
                upper = in_args.penaltyRange(idx.optPerf(1));%current choice
                newpenRange = lower:upper/6:upper;
                newpenRange = sort([newpenRange in_args.penaltyRange(idx.optPerf(1))]);%try 7 penalties disributed within this range
                for s = 1:in_args.nFoldsPenaltySelection %assumes balanced data
                    
                    thisIdx1 = randomNFold1==s;
                    
                    omitc1 = idx.c1(randomNFold1==s);
                    omitc2 = idx.c2(randomNFold2==s);
                    omitc3 = idx.c3(randomNFold3==s);
                    omitc4 = idx.c4(randomNFold4==s);
                    omitc5 = idx.c5(randomNFold5==s);
                    
                    theseOmitted = [omitc1; omitc2; omitc3; omitc4; omitc5 ];
                    
                    thisChoice = choice;
                    thisChoice(theseOmitted) = [];
                    
                    thisV = v;
                    thisV(theseOmitted,:) = [];
                    
                    
                    for i = 1:length(newpenRange)
                        l = newpenRange(i);
                        
                        if doRadialSearch
                            rSet = in_args.radialBasisSelection;
                        else
                            rSet = 0;
                        end
                        
                        for r = 1:length(rSet)
                            thisR = rSet(r);
                            
                            if strcmp(in_args.classType,'libLin')
                                trainOpts_orig = in_args.libLin ;
                                trainOptsOptimize = [trainOpts_orig ' -c ' num2str(l)];
                                m = ll_train(thisChoice, sparse(thisV), trainOptsOptimize);%NOTE: changed this from "train(...)" because I needed to refer to mac-friendly compiled files Alan put together for Valerie
                                [theseLabels acc]=ll_predict(choice(theseOmitted), sparse(v(theseOmitted,:)), m);
                                guesses(theseOmitted,i,r) = theseLabels;
                                accs(s,i,r) = acc;
                            end
                        end
                    end
                    
                end
                accbypen = mean(accs,1);
                optPerfs = max(max(accbypen));
                idx.optPerf = find(accbypen==optPerfs);
                opt_penalty = in_args.penaltyRange(idx.optPerf(1));
                
            else
                lower = (in_args.penaltyRange(idx.optPerf(1))-in_args.penaltyRange(idx.optPerf(1)-1))/2;%halfway to next penalty
                upper = (in_args.penaltyRange(idx.optPerf(1)+1)-in_args.penaltyRange(idx.optPerf(1)))/2;%halfway to next penalty
                newpenRange = lower:upper/6:upper;
                newpenRange = sort([newpenRange in_args.penaltyRange(idx.optPerf(1))]);%try 7 penalties disributed within this range
                for s = 1:in_args.nFoldsPenaltySelection %assumes balanced data
                    
                    thisIdx1 = randomNFold1==s;
                    
                    omitc1 = idx.c1(randomNFold1==s);
                    omitc2 = idx.c2(randomNFold2==s);
                    omitc3 = idx.c3(randomNFold3==s);
                    omitc4 = idx.c4(randomNFold4==s);
                    omitc5 = idx.c5(randomNFold5==s);
                    
                    theseOmitted = [omitc1; omitc2; omitc3; omitc4; omitc5 ];
                    
                    thisChoice = choice;
                    thisChoice(theseOmitted) = [];
                    
                    thisV = v;
                    thisV(theseOmitted,:) = [];
                    
                    
                    for i = 1:length(newpenRange)
                        l = newpenRange(i);
                        
                        if doRadialSearch
                            rSet = in_args.radialBasisSelection;
                        else
                            rSet = 0;
                        end
                        
                        for r = 1:length(rSet)
                            thisR = rSet(r);
                            
                            if strcmp(in_args.classType,'libLin')
                                trainOpts_orig = in_args.libLin ;
                                trainOptsOptimize = [trainOpts_orig ' -c ' num2str(l)];
                                m = ll_train(thisChoice, sparse(thisV), trainOptsOptimize);%NOTE: changed this from "train(...)" because I needed to refer to mac-friendly compiled files Alan put together for Valerie
                                [theseLabels acc]=ll_predict(choice(theseOmitted), sparse(v(theseOmitted,:)), m);
                                guesses(theseOmitted,i,r) = theseLabels;
                                accs(s,i,r) = acc;
                            end
                        end
                    end
                    
                end
                accbypen = mean(accs,1);
                optPerfs = max(max(accbypen));
                idx.optPerf = find(accbypen==optPerfs);
                opt_penalty = in_args.penaltyRange(idx.optPerf(1));
                
            end
            
            
            
            
        else
            opt_penalty = in_args.penaltyRange(idx.optPerf(1));
            
            
        end
        
    else
        choiceMat = repmat(choice,[1,length(in_args.penaltyRange), length(rSet)]);
        % performance measures~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %"all the rest" vectors reflect
        rest1 = [2 3 4 5];
        rest2 = [1 3 4 5];
        rest3 = [1 2 4 5];
        rest4 = [1 2 3 5];
        rest5 = [1 2 3 4];
        
        perf.tp1 = sum(guesses == choice_set(1) & choiceMat == choice_set(1));   % true pos
        perf.fp1 = zeros(size(perf.tp1));
        for x = 1:numel(rest1);
            perf.fp1 = perf.fp1+sum(guesses == choice_set(1) & choiceMat == choice_set(rest1(x)));
        end
        %perf.fp1 = sum(guesses == choice_set(1) & choiceMat == choice_set(2) | choice_set(3) | choice_set(4) | choice_set(5));  % false pos
        perf.fn1 = zeros(size(perf.tp1));
        for x = 1:numel(rest1);
            perf.fn1 = perf.fn1+sum(guesses == choice_set(rest1(x)) & choiceMat == choice_set(1));
        end
        %perf.fn = sum(guesses == choice_set(2) & choiceMat == choice_set(1));   % false neg
        perf.tn1 = zeros(size(perf.tp1));
        for x = 1:numel(rest1);
            perf.tn1 = perf.tn1+sum(guesses == choice_set(rest1(x)) & choiceMat == choice_set(rest1(x)));%collect true postives for the other classes
        end
        perf.tn1 = perf.tn1+sum(guesses == choice_set(3) & choiceMat == choice_set(2));
        perf.tn1 = perf.tn1+sum(guesses == choice_set(4) & choiceMat == choice_set(2));
        perf.tn1 = perf.tn1+sum(guesses == choice_set(5) & choiceMat == choice_set(2));
        perf.tn1 = perf.tn1+sum(guesses == choice_set(2) & choiceMat == choice_set(3));
        perf.tn1 = perf.tn1+sum(guesses == choice_set(4) & choiceMat == choice_set(3));
        perf.tn1 = perf.tn1+sum(guesses == choice_set(5) & choiceMat == choice_set(3));
        perf.tn1 = perf.tn1+sum(guesses == choice_set(2) & choiceMat == choice_set(4));
        perf.tn1 = perf.tn1+sum(guesses == choice_set(3) & choiceMat == choice_set(4));
        perf.tn1 = perf.tn1+sum(guesses == choice_set(5) & choiceMat == choice_set(4));
        perf.tn1 = perf.tn1+sum(guesses == choice_set(2) & choiceMat == choice_set(5));
        perf.tn1 = perf.tn1+sum(guesses == choice_set(3) & choiceMat == choice_set(5));
        perf.tn1 = perf.tn1+sum(guesses == choice_set(4) & choiceMat == choice_set(5));
        %perf.tn = sum(guesses == choice_set(2) & choiceMat == choice_set(2));  % true neg
        
        perf.tp2 = sum(guesses == choice_set(2) & choiceMat == choice_set(2));   % true pos
        perf.fp2 = zeros(size(perf.tp2));
        for x = 1:numel(rest2);
            perf.fp2 = perf.fp2+sum(guesses == choice_set(2) & choiceMat == choice_set(rest2(x)));
        end
        %perf.fp1 = sum(guesses == choice_set(1) & choiceMat == choice_set(2) | choice_set(3) | choice_set(4) | choice_set(5));  % false pos
        perf.fn2 = zeros(size(perf.tp2));
        for x = 1:numel(rest2);
            perf.fn2 = perf.fn2+sum(guesses == choice_set(rest2(x)) & choiceMat == choice_set(2));
        end
        %perf.fn = sum(guesses == choice_set(2) & choiceMat == choice_set(1));   % false neg
        perf.tn2 = zeros(size(perf.tp2));
        for x = 1:numel(rest2);
            perf.tn2 = perf.tn2+sum(guesses == choice_set(rest2(x)) & choiceMat == choice_set(rest2(x)));%collect true postives for the other classes
        end
        perf.tn2 = perf.tn2+sum(guesses == choice_set(3) & choiceMat == choice_set(1));
        perf.tn2 = perf.tn2+sum(guesses == choice_set(4) & choiceMat == choice_set(1));
        perf.tn2 = perf.tn2+sum(guesses == choice_set(5) & choiceMat == choice_set(1));
        perf.tn2 = perf.tn2+sum(guesses == choice_set(1) & choiceMat == choice_set(3));
        perf.tn2 = perf.tn2+sum(guesses == choice_set(4) & choiceMat == choice_set(3));
        perf.tn2 = perf.tn2+sum(guesses == choice_set(5) & choiceMat == choice_set(3));
        perf.tn2 = perf.tn2+sum(guesses == choice_set(1) & choiceMat == choice_set(4));
        perf.tn2 = perf.tn2+sum(guesses == choice_set(3) & choiceMat == choice_set(4));
        perf.tn2 = perf.tn2+sum(guesses == choice_set(5) & choiceMat == choice_set(4));
        perf.tn2 = perf.tn2+sum(guesses == choice_set(1) & choiceMat == choice_set(5));
        perf.tn2 = perf.tn2+sum(guesses == choice_set(3) & choiceMat == choice_set(5));
        perf.tn2 = perf.tn2+sum(guesses == choice_set(4) & choiceMat == choice_set(5));
        %perf.tn = sum(guesses == choice_set(2) & choiceMat == choice_set(2));  % true neg
        
        perf.tp3 = sum(guesses == choice_set(3) & choiceMat == choice_set(3));   % true pos
        perf.fp3 = zeros(size(perf.tp3));
        for x = 1:numel(rest3);
            perf.fp3 = perf.fp3+sum(guesses == choice_set(3) & choiceMat == choice_set(rest3(x)));
        end
        %perf.fp3 = sum(guesses == choice_set(1) & choiceMat == choice_set(2) | choice_set(3) | choice_set(4) | choice_set(5));  % false pos
        perf.fn3 = zeros(size(perf.tp3));
        for x = 1:numel(rest3);
            perf.fn3 = perf.fn3+sum(guesses == choice_set(rest3(x)) & choiceMat == choice_set(3));
        end
        %perf.fn = sum(guesses == choice_set(2) & choiceMat == choice_set(1));   % false neg
        perf.tn3 = zeros(size(perf.tp3));
        for x = 1:numel(rest3);
            perf.tn3 = perf.tn3+sum(guesses == choice_set(rest3(x)) & choiceMat == choice_set(rest3(x)));%collect true postives for the other classes
        end
        perf.tn3 = perf.tn3+sum(guesses == choice_set(2) & choiceMat == choice_set(1));
        perf.tn3 = perf.tn3+sum(guesses == choice_set(4) & choiceMat == choice_set(1));
        perf.tn3 = perf.tn3+sum(guesses == choice_set(5) & choiceMat == choice_set(1));
        perf.tn3 = perf.tn3+sum(guesses == choice_set(1) & choiceMat == choice_set(2));
        perf.tn3 = perf.tn3+sum(guesses == choice_set(4) & choiceMat == choice_set(2));
        perf.tn3 = perf.tn3+sum(guesses == choice_set(5) & choiceMat == choice_set(2));
        perf.tn3 = perf.tn3+sum(guesses == choice_set(2) & choiceMat == choice_set(4));
        perf.tn3 = perf.tn3+sum(guesses == choice_set(1) & choiceMat == choice_set(4));
        perf.tn3 = perf.tn3+sum(guesses == choice_set(5) & choiceMat == choice_set(4));
        perf.tn3 = perf.tn3+sum(guesses == choice_set(2) & choiceMat == choice_set(5));
        perf.tn3 = perf.tn3+sum(guesses == choice_set(1) & choiceMat == choice_set(5));
        perf.tn3 = perf.tn3+sum(guesses == choice_set(4) & choiceMat == choice_set(5));
        %perf.tn = sum(guesses == choice_set(2) & choiceMat == choice_set(2));  % true neg
        
        perf.tp4 = sum(guesses == choice_set(4) & choiceMat == choice_set(4));   % true pos
        perf.fp4 = zeros(size(perf.tp4));
        for x = 1:numel(rest4);
            perf.fp4 = perf.fp4+sum(guesses == choice_set(4) & choiceMat == choice_set(rest4(x)));
        end
        %perf.fp1 = sum(guesses == choice_set(1) & choiceMat == choice_set(2) | choice_set(3) | choice_set(4) | choice_set(5));  % false pos
        perf.fn4 = zeros(size(perf.tp4));
        for x = 1:numel(rest4);
            perf.fn4 = perf.fn4+sum(guesses == choice_set(rest4(x)) & choiceMat == choice_set(4));
        end
        %perf.fn = sum(guesses == choice_set(2) & choiceMat == choice_set(1));   % false neg
        perf.tn4 = zeros(size(perf.tp4));
        for x = 1:numel(rest4);
            perf.tn4 = perf.tn4+sum(guesses == choice_set(rest4(x)) & choiceMat == choice_set(rest4(x)));%collect true postives for the other classes
        end
        perf.tn4 = perf.tn4+sum(guesses == choice_set(3) & choiceMat == choice_set(2));
        perf.tn4 = perf.tn4+sum(guesses == choice_set(1) & choiceMat == choice_set(2));
        perf.tn4 = perf.tn4+sum(guesses == choice_set(5) & choiceMat == choice_set(2));
        perf.tn4 = perf.tn4+sum(guesses == choice_set(2) & choiceMat == choice_set(3));
        perf.tn4 = perf.tn4+sum(guesses == choice_set(1) & choiceMat == choice_set(3));
        perf.tn4 = perf.tn4+sum(guesses == choice_set(5) & choiceMat == choice_set(3));
        perf.tn4 = perf.tn4+sum(guesses == choice_set(2) & choiceMat == choice_set(1));
        perf.tn4 = perf.tn4+sum(guesses == choice_set(3) & choiceMat == choice_set(1));
        perf.tn4 = perf.tn4+sum(guesses == choice_set(5) & choiceMat == choice_set(1));
        perf.tn4 = perf.tn4+sum(guesses == choice_set(2) & choiceMat == choice_set(5));
        perf.tn4 = perf.tn4+sum(guesses == choice_set(3) & choiceMat == choice_set(5));
        perf.tn4 = perf.tn4+sum(guesses == choice_set(1) & choiceMat == choice_set(5));
        %perf.tn = sum(guesses == choice_set(2) & choiceMat == choice_set(2));  % true neg
        
        perf.tp5 = sum(guesses == choice_set(5) & choiceMat == choice_set(5));   % true pos
        perf.fp5 = zeros(size(perf.tp5));
        for x = 1:numel(rest5);
            perf.fp5 = perf.fp5+sum(guesses == choice_set(5) & choiceMat == choice_set(rest5(x)));
        end
        %perf.fp1 = sum(guesses == choice_set(1) & choiceMat == choice_set(2) | choice_set(3) | choice_set(4) | choice_set(5));  % false pos
        perf.fn5 = zeros(size(perf.tp5));
        for x = 1:numel(rest5);
            perf.fn5 = perf.fn5+sum(guesses == choice_set(rest5(x)) & choiceMat == choice_set(5));
        end
        %perf.fn = sum(guesses == choice_set(2) & choiceMat == choice_set(1));   % false neg
        perf.tn5 = zeros(size(perf.tp5));
        for x = 1:numel(rest5);
            perf.tn5 = perf.tn5+sum(guesses == choice_set(rest5(x)) & choiceMat == choice_set(rest5(x)));%collect true postives for the other classes
        end
        perf.tn5 = perf.tn5+sum(guesses == choice_set(3) & choiceMat == choice_set(2));
        perf.tn5 = perf.tn5+sum(guesses == choice_set(4) & choiceMat == choice_set(2));
        perf.tn5 = perf.tn5+sum(guesses == choice_set(1) & choiceMat == choice_set(2));
        perf.tn5 = perf.tn5+sum(guesses == choice_set(2) & choiceMat == choice_set(3));
        perf.tn5 = perf.tn5+sum(guesses == choice_set(4) & choiceMat == choice_set(3));
        perf.tn5 = perf.tn5+sum(guesses == choice_set(1) & choiceMat == choice_set(3));
        perf.tn5 = perf.tn5+sum(guesses == choice_set(2) & choiceMat == choice_set(4));
        perf.tn5 = perf.tn5+sum(guesses == choice_set(3) & choiceMat == choice_set(4));
        perf.tn5 = perf.tn5+sum(guesses == choice_set(1) & choiceMat == choice_set(4));
        perf.tn5 = perf.tn5+sum(guesses == choice_set(2) & choiceMat == choice_set(1));
        perf.tn5 = perf.tn5+sum(guesses == choice_set(3) & choiceMat == choice_set(1));
        perf.tn5 = perf.tn5+sum(guesses == choice_set(4) & choiceMat == choice_set(1));
        %perf.tn = sum(guesses == choice_set(2) & choiceMat == choice_set(2));  % true neg
        % quality metrics
        perf.Precision1 = perf.tp1./(perf.tp1+perf.fp1);
        perf.Recall1 = perf.tp1./(perf.tp1+perf.fn1);
        
        perf.TrueNegRate1 = perf.tn1./(perf.tn1+perf.fp1);%specificity
        %perf.Accuracy = (perf.tp+perf.tn)./ ...
        %    (perf.tp+perf.tn+perf.fp+perf.fn);
        
        perf.AUC1 = (perf.Recall1+perf.TrueNegRate1)*.5;%compute binary class AUC ((recall+spec)/2)
        
        perf.F_Score1 = 2.*perf.Precision1.*perf.Recall1./ ...
            (perf.Precision1+perf.Recall1);
        
        
        perf.Precision2 = perf.tp2./(perf.tp2+perf.fp2);
        perf.Recall2 = perf.tp2./(perf.tp2+perf.fn2);
        
        perf.TrueNegRate2 = perf.tn2./(perf.tn2+perf.fp2);%specificity
        %perf.Accuracy = (perf.tp+perf.tn)./ ...
        %    (perf.tp+perf.tn+perf.fp+perf.fn);
        
        perf.AUC2 = (perf.Recall2+perf.TrueNegRate2)*.5;%compute binary class AUC ((recall+spec)/2)
        
        perf.F_Score2 = 2.*perf.Precision2.*perf.Recall2./ ...
            (perf.Precision2+perf.Recall2);
        
        
        perf.Precision3 = perf.tp3./(perf.tp3+perf.fp3);
        perf.Recall3 = perf.tp3./(perf.tp3+perf.fn3);
        
        perf.TrueNegRate3 = perf.tn3./(perf.tn3+perf.fp3);%specificity
        %perf.Accuracy = (perf.tp+perf.tn)./ ...
        %    (perf.tp+perf.tn+perf.fp+perf.fn);
        
        perf.AUC3 = (perf.Recall3+perf.TrueNegRate3)*.5;%compute binary class AUC ((recall+spec)/2)
        
        perf.F_Score3 = 2.*perf.Precision3.*perf.Recall3./ ...
            (perf.Precision3+perf.Recall3);
        
        
        perf.Precision4 = perf.tp4./(perf.tp4+perf.fp4);
        perf.Recall4 = perf.tp4./(perf.tp4+perf.fn4);
        
        perf.TrueNegRate4 = perf.tn4./(perf.tn4+perf.fp4);%specificity
        %perf.Accuracy = (perf.tp+perf.tn)./ ...
        %    (perf.tp+perf.tn+perf.fp+perf.fn);
        
        perf.AUC4 = (perf.Recall4+perf.TrueNegRate4)*.5;%compute binary class AUC ((recall+spec)/2)
        
        perf.F_Score4 = 2.*perf.Precision4.*perf.Recall4./ ...
            (perf.Precision4+perf.Recall4);
        
        
        perf.Precision5 = perf.tp5./(perf.tp5+perf.fp5);
        perf.Recall5 = perf.tp5./(perf.tp5+perf.fn5);
        
        perf.TrueNegRate5 = perf.tn5./(perf.tn5+perf.fp5);%specificity
        %perf.Accuracy = (perf.tp+perf.tn)./ ...
        %    (perf.tp+perf.tn+perf.fp+perf.fn);
        
        perf.AUC5 = (perf.Recall5+perf.TrueNegRate5)*.5;%compute binary class AUC ((recall+spec)/2)
        
        perf.F_Score5 = 2.*perf.Precision5.*perf.Recall5./ ...
            (perf.Precision5+perf.Recall5);
        
        % use AUC or F score as the performance measure to select the optimal
        % parameter
        
        perfParam1 = perf.AUC1;
        
        optPerf1 = max(max(perfParam1));
        perfParam1 = squeeze(perfParam1);
        
        perfParam2 = perf.AUC2;
        
        optPerf2 = max(max(perfParam2));
        perfParam2 = squeeze(perfParam2);
        
        
        perfParam3 = perf.AUC3;
        
        optPerf3 = max(max(perfParam3));
        perfParam3 = squeeze(perfParam3);
        
        perfParam4 = perf.AUC4;
        
        optPerf4 = max(max(perfParam4));
        perfParam4 = squeeze(perfParam4);
        
        perfParam5 = perf.AUC5;
        
        optPerf5 = max(max(perfParam5));
        perfParam5 = squeeze(perfParam5);
        
        if doRadialSearch
            [idx.optPerf1 idx.optPerf2] = find(perfParam==optPerf);
            opt_penalty = in_args.penaltyRange(idx.optPerf1(1));
            opt_r = in_args.radialBasisSelection(idx.optPerf2(1));
        else
            
            stacked_perf = cat(3,perfParam1, perfParam2, perfParam3, perfParam4, perfParam5);
            %tot_perf = perfParam1 + perfParam2 + perfParam3 + perfParam4 + perfParam5;
            
            if strcmp(in_args.multiclassPenaltySummary,'mean')
                tot_perf = mean(stacked_perf,3);
            elseif strcmp(in_args.multiclassPenaltySummary,'median')
                tot_perf = median(stacked_perf,3);
            end
            optPerfs = max(max(tot_perf));
            idx.optPerf = find(tot_perf==optPerfs);
            
            %         idx.optPerf1 = find(perfParam1==optPerf1);
            %         idx.optPerf2 = find(perfParam2==optPerf2);
            %         idx.optPerf3 = find(perfParam3==optPerf3);
            %         idx.optPerf4 = find(perfParam4==optPerf4);
            %         idx.optPerf5 = find(perfParam5==optPerf5);
            
            opt_penalty = in_args.penaltyRange(idx.optPerf(1));
        end
    end
end
%scratchpad.opt_penalty = opt_penalty;



