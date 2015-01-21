function subj = TB_balanceTrainPats(S, subj)
% modification of PM_balanceTrainPats - includes code to balance based on
% two factors (specifically, designed for balancing number of trials based
% on both the goal location AND the cues for the goal location)

% Rev 1.1 - 1/20/2015 - fixed bug in TrainTestOneIter balancing, added
% option to hardcode max trial-count in the training set

%flags
setmax = 0;%if a non-zero value is specified, override number of cues below to have a training trial-count of 2x this value (e.g. 6 = 12 training trials per cond)

%files and data
cuefile = ['/mnt/wagner/thackery/Circmaze/subs_preprocessed_fermi_nospikerepair/' S.subj_id '/model_plan_allruns_includemarginals/cues_' S.subj_id '_allruns.mat']
load(cuefile)
sels = squeeze(get_group_as_matrix(subj, 'selector',S.thisSelector));%create matrix of bins/iterations*trials (1s = training, 2s = testing trials)
regs = get_mat(subj, 'regressors', 'conds');%get matrix of conditions/classes*trials (1s = trial X = this class, 0s = trial X ~= this class)

if strcmp(S.thisSelector, 'TrainTestOneIterGroup')%for some reason, the selectors need to be translated to work with "train on one set, test on another" style classification
    sels = sels';
    
    %sort the regs by selector to ensure testing set is temporarily at the end (this is
    %because we want to balance the cues are at half the number)
    %[x ind] = sort(sels);
    %y = regs(:,ind);
    %regst = y;
    %sels = x;
    
    % convert cues into reg space (optional - developed for TIB Circmaze study)
    allcues = zeros(size(regs)); %initialize empty cues matrix
    
    
    for r = 1:size(regs,1)
        %testidx
        regstrain = regs(r,:).*sels; %S.trainActives) %zero out test trials temporarily
        regidx = find(regstrain==1);
        %insert each cue value into its corresponding allcues cell
        for c = 1:length(cues{r})
            allcues(r,regidx(c)) = cues{r}(c);
        end
    end
    %regidx = find(sels==1);
    
    %% begin balancing
    for i = 1:S.numBalancedIts%runs the balancing multiple times
        for s = 1:size(sels,1)
            selIt = (sels(s,:)); %row 's' of sels
            trainRegIdx = find(selIt==1); %index training trials for this iteration
            testRegIdx = find(selIt==2); %index testing trials for this iteration
            
            trainregs = regs(:,selIt==1); %filter 'regs' down into training trial subset for this iteration
            % cue balancing..........................
            traincues = allcues(:,selIt==1); %filter 'cues' down into training trial subset for this iteration
            % end cue balancing......................
            nTrainPats = selIt(selIt==1) * trainregs'; %calculate # trials per condition in the training set
            
            newBinSize = min(nTrainPats(nTrainPats~=0)); %the min number of trials in the class bins
            active_trials_train = zeros(1,size(trainregs,2)); %creates vector of zeros of length of training trial subset
            
            for t = 1:length(nTrainPats)
                if (nTrainPats(t)>0)
                    theseTrials_h = shuffle(find((trainregs(t,:))==1)); % index 1s amongst _ training trials, then randomize (shuffle) the order
                    theseTrials = trainRegIdx(theseTrials_h); % now load up the index among ALL (train and test) trials corresponding to these training trials
                    idxInclude{t} = theseTrials(1:newBinSize); % store these trials, thresholded at min trial count, in new array
                    
                    % cue balancing..........................
                    newTraincues = traincues(t,theseTrials_h);
                    newTraincues_short = newTraincues(1:newBinSize);
                    idxbycue{1,t} = idxInclude{t}(newTraincues_short==1);
                    idxbycue{2,t} = idxInclude{t}(newTraincues_short==2);
                    % end cue balancing......................
                    
                    active_trials_train(idxInclude{t}) = 1; % fill in 1s to the new "replacement" active training trials selector
                end
            end
            
            % cue balancing..........................
            fewestcues = min(min(cellfun('size', idxbycue, 2)));%find the cue indices with the fewest instances
            for x = 1:length(idxInclude)
                idxInclude2{x} = [idxbycue{1,x}(1:fewestcues) idxbycue{2,x}(1:fewestcues)];
            end
            % end cue balancing......................
            
            %combine training and testing indices to create new "active_trials"
            %selector
            allIncluded = sort([idxInclude2{:}]);%allIncluded = sort([idxInclude{:}]);
            active_trials = zeros(1,size(selIt,2));
            active_trials(allIncluded) = 1;
            active_trials(testRegIdx) = 1;
            %         NDiff = diff(nTrainPats);
            %         active_trials = sum(regs);
            %
            %         if diff(nTrainPats)<0
            %             theseTrials_h = shuffle(find((trainregs(1,:))==1));
            %             theseTrials = trainRegIdx(theseTrials_h);
            %             idxRemove = theseTrials(1:abs(NDiff));
            %             active_trials(idxRemove) = 0;
            %         elseif diff(nTrainPats)>0
            %             theseTrials_h = shuffle(find((trainregs(2,:))==1));
            %             theseTrials = trainRegIdx(theseTrials_h);
            %             idxRemove = theseTrials(1:abs(NDiff));
            %             active_trials(idxRemove) = 0;
            %         end
            
            selItName = [S.thisSelector num2str(s) '_balanced' num2str(i)];
            selItGroupName = [S.thisSelector 'balanced'];
            
            subj = init_object(subj,'selector', selItName);
            subj = set_mat(subj,'selector',selItName, selIt .* active_trials);
            subj = set_objfield(subj, 'selector', selItName, 'group_name', selItGroupName);
            
        end
    end
    
    
else
    
    
    %% convert cues into reg space (optional - developed for TIB Circmaze study)
    allcues = zeros(size(regs)); %initialize empty cues matrix
    for r = 1:size(regs,1)
        regidx = find((regs(r,:))==1);
        %insert each cue value into its corresponding allcues cell
        for c = 1:length(cues{r})
            allcues(r,regidx(c)) = cues{r}(c);
        end
    end
    
    
    %% begin balancing
    for i = 1:S.numBalancedIts%runs the balancing multiple times
        for s = 1:size(sels,1)
            selIt = (sels(s,:)); %row 's' of sels
            trainRegIdx = find(selIt==1); %index training trials for this iteration
            testRegIdx = find(selIt==2); %index testing trials for this iteration
            
            trainregs = regs(:,selIt==1); %filter 'regs' down into training trial subset for this iteration
            % cue balancing..........................
            traincues = allcues(:,selIt==1); %filter 'cues' down into training trial subset for this iteration
            % end cue balancing......................
            nTrainPats = selIt(selIt==1) * trainregs'; %calculate # trials per condition in the training set
            
            newBinSize = min(nTrainPats(nTrainPats~=0)); %the min number of trials in the class bins
            active_trials_train = zeros(1,size(trainregs,2)); %creates vector of zeros of length of training trial subset
            
            for t = 1:length(nTrainPats)
                if (nTrainPats(t)>0)
                    theseTrials_h = shuffle(find((trainregs(t,:))==1)); % index 1s amongst _ training trials, then randomize (shuffle) the order
                    theseTrials = trainRegIdx(theseTrials_h); % now load up the index among ALL (train and test) trials corresponding to these training trials
                    idxInclude{t} = theseTrials(1:newBinSize); % store these trials, thresholded at min trial count, in new array
                    
                    % cue balancing..........................
                    newTraincues = traincues(t,theseTrials_h);
                    newTraincues_short = newTraincues(1:newBinSize);
                    idxbycue{1,t} = idxInclude{t}(newTraincues_short==1);
                    idxbycue{2,t} = idxInclude{t}(newTraincues_short==2);
                    % end cue balancing......................
                    
                    active_trials_train(idxInclude{t}) = 1; % fill in 1s to the new "replacement" active training trials selector
                end
            end
            
            % cue balancing..........................
            if setmax > 0
                fewestcues = setmax;
            else
                fewestcues = min(min(cellfun('size', idxbycue, 2)));%find the cue indices with the fewest instances
            end
            
            for x = 1:length(idxInclude)
                idxInclude2{x} = [idxbycue{1,x}(1:fewestcues) idxbycue{2,x}(1:fewestcues)];
            end
            % end cue balancing......................
            
            %combine training and testing indices to create new "active_trials"
            %selector
            allIncluded = sort([idxInclude2{:}]);%allIncluded = sort([idxInclude{:}]);
            active_trials = zeros(1,size(selIt,2));
            active_trials(allIncluded) = 1;
            active_trials(testRegIdx) = 1;
            %         NDiff = diff(nTrainPats);
            %         active_trials = sum(regs);
            %
            %         if diff(nTrainPats)<0
            %             theseTrials_h = shuffle(find((trainregs(1,:))==1));
            %             theseTrials = trainRegIdx(theseTrials_h);
            %             idxRemove = theseTrials(1:abs(NDiff));
            %             active_trials(idxRemove) = 0;
            %         elseif diff(nTrainPats)>0
            %             theseTrials_h = shuffle(find((trainregs(2,:))==1));
            %             theseTrials = trainRegIdx(theseTrials_h);
            %             idxRemove = theseTrials(1:abs(NDiff));
            %             active_trials(idxRemove) = 0;
            %         end
            
            selItName = [S.thisSelector num2str(s) '_balanced' num2str(i)];
            selItGroupName = [S.thisSelector 'balanced'];
            
            subj = init_object(subj,'selector', selItName);
            subj = set_mat(subj,'selector',selItName, selIt .* active_trials);
            subj = set_objfield(subj, 'selector', selItName, 'group_name', selItGroupName);
            
        end
    end
    
end