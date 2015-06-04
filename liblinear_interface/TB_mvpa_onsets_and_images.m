function [S] = TB_mvpa_onsets_and_images(S)
%derived from Alan Gordon's PM_mvpa_onsets_and_images. Slight tweaks,
%particularly to flexibly handle both "raw" and beta series classification
%scenarios

if S.xval
    [onsets_train, S.idxOnsets_train_in_classifier] = MakeOnsets(S, 'train');
    [onsets_test, S.idxOnsets_test_in_classifier] = MakeOnsets(S, 'test');
    
    for k=1:length(onsets_train)
        S.onsets{k} = union(onsets_train{k}, onsets_test{k});
    end
else
    [onsets_train, S.idxOnsets_train_in_classifier] = MakeOnsets(S, 'train_set');%MakeOnsets(S, 'train');
    [onsets_test, S.idxOnsets_test_in_classifier] = MakeOnsets(S, 'test_set');%MakeOnsets(S, 'test');
    
    for k=1:length(onsets_train)
        %S.onsets{k} = union(onsets_train{k}, S.durTrain + onsets_test{k});%version where the test onsets "restart" at 0s
        S.onsets{k} = union(onsets_train{k}, onsets_test{k});
    end
end

S.onsets_train_in_classifier = onsets_train;
S.onsets_test_in_classifier = onsets_test;

end

function [onsets idxOnsets_in_classifier] = MakeOnsets(S, part)
%works for x-validation, because currently assigns same onsets to both sets
if strcmp(part, 'train')
    thisOns = load(fullfile(S.onsetsTrainDir, S.onsets_filename));%good, comes from PM_mvpa_params
    if strcmp(S.inputformat, 'betas')
        thisbetidx = load(fullfile(S.onsetsTrainDir, S.betaidx_filename));%load the betas_idx file
    end
    thisConds = S.condsTrain;%good, comes from PM_mvpa_params
    theseFiles = S.filenames_train;%good, comes from PM_mvpa_params
    theseTRWeights = S.TR_weights{1};%good, comes from PM_mvpa_params
    thisIdx = S.idxTr;%-----------------what to doooooo
elseif strcmp(part, 'test')
    thisOns = load(fullfile(S.onsetsTestDir, S.onsets_filename));%good, comes from PM_mvpa_params
    if strcmp(S.inputformat, 'betas')
        thisbetidx = load(fullfile(S.onsetsTestDir, S.betaidx_filename));%load the betas_idx file
    end
    thisConds = S.condsTest;%good, comes from PM_mvpa_params
    theseFiles = S.filenames_test;%good, comes from PM_mvpa_params
    theseTRWeights = S.TR_weights{2};%good, comes from PM_mvpa_params
    thisIdx = S.idxTe;%-----------------what to doooooo
%added for train on one set, test on the other set circumstance -
%6/11/2014, TIB
elseif strcmp(part, 'train_set')
    thisOns = load(fullfile(S.onsetsTrainDir, S.onsets_filename_tr));%good, comes from PM_mvpa_params
    if strcmp(S.inputformat, 'betas')
        thisbetidx = load(fullfile(S.onsetsTrainDir, S.betaidx_filename_tr));%load the betas_idx file
        thisbetidx.betaidx = thisbetidx.betaidx_tr
    end
    thisConds = S.condsTrain;%good, comes from PM_mvpa_params
    theseFiles = S.filenames_train;%good, comes from PM_mvpa_params
    theseTRWeights = S.TR_weights{1};%good, comes from PM_mvpa_params
    thisIdx = S.idxTr;%-----------------what to doooooo
elseif strcmp(part, 'test_set')
    thisOns = load(fullfile(S.onsetsTestDir, S.onsets_filename_tst));%good, comes from PM_mvpa_params
    if strcmp(S.inputformat, 'betas')
        thisbetidx = load(fullfile(S.onsetsTestDir, S.betaidx_filename_te));%load the betas_idx file
        thisbetidx.betaidx = thisbetidx.betaidx_te
    end
    thisConds = S.condsTest;%good, comes from PM_mvpa_params
    theseFiles = S.filenames_test;%good, comes from PM_mvpa_params
    theseTRWeights = S.TR_weights{2};%good, comes from PM_mvpa_params
    thisIdx = S.idxTe;%-----------------what to doooooo    

end





onsets = cell(size(thisConds));
%if ~strcmp(S.patternType, 'betas')%if we are dealing with betas, do the following
for i = 1:length(thisConds)%e.g. 5, for {{one} {two} {three} {four} {five}}
    for j = 1:length(thisConds{i})
        idxThisCond = find(strcmp(thisOns.names, thisConds{i}{j}));
        %betidxnum = find(strcmp(thisbetidx.bnames, thisConds{i}{j}));%finds idx in thisbetidx corresponding to current condition name, added 1/7/2015
        if ~isempty(idxThisCond)
            %enoughTRs{i} = logical(vertcat(enoughTRs{i}, enoughTRs_h));
            %            theseOnsets = asRow(thisOns.onsets{idxThisCond});%spits the
            %onsets into a row vector format (if not already)
            
            if strcmp(S.patternType, 'betas') %we are going to convert the logical indexing for betas into "onsets" (i.e. image numbers)
                
                c = 1
                for x = 1:numel(thisbetidx.betaidx{i})%x = 1:numel(thisbetidx.betaidx{betidxnum})%
                    if thisbetidx.betaidx{i}(x) == 1%if thisbetidx.betaidx{betidxnum}(x) == 1%
                        onsets{i}(c) = x
                        c = c + 1
                    end
                end
                %onsets{i} = sort(horzcat(betaidx{i}, theseOnsets));
            else
                time_idx = floor(thisOns.onsets{idxThisCond}/S.TR) + 1;
                enoughTRs_h = (time_idx + length(theseTRWeights)) <= length(theseFiles);
                
                %onsets{i} = sort(horzcat(onsets{i}, theseOnsets(enoughTRs_h)));
                onsets{i} = sort(horzcat(onsets{i}, thisOns.onsets{idxThisCond}(enoughTRs_h)));%put the onsets for cond{i} into an array,
                %restricting to onsets for which we have enough TRs acquired.
                %NOTE: This will literally REcreate the onsets component of my model files,
                %BUT is still useful IF we haven't prescreened them based
                %on how many TRs we acquired (i.e. we might get 31/32
                %onsets of the 32nd exceeds the number of TRs we acquired,
                %based on how many post-onset TRs we want to average across
                %in our MVPA analysis
            end
        else
            error('SANITY CHECK: the condition names you are referencing do not appear to exist in the file loaded')
        end
    end
end
%end

%idxOnsets_in_classifier = asRow(ismember(thisIdx.alltrials, [onsets{:}]));%ismember is making a vector of 1s and 0s, with 1s indicating onsets that exist in both thisIdx.alltrials and our onsets array.
%We just want a row vector of all 1s for trials for which onsets are listed here.

if strcmp(S.patternType, 'betas')%added by TIB 10/27/2014 to handle binary classification with betas
     idxOnsets_in_classifier = zeros([1,length(thisOns.names)]);
     for z = 1:length(onsets)
         idxOnsets_in_classifier(onsets{z}) = 1;
     end
else
idxOnsets_in_classifier = ones([1,numel([onsets{:}])]);
end

end
