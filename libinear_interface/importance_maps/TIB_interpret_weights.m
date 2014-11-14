function [subj] = TIB_interpret_weights(subj,results,results_IW,S, varargin)

% USAGE : [subj] = interpret_weights(subj,results,varargin)
%
% This algorithm is used to determine the level of influence
% of a given voxel to a given category. The value produced
% by the algorithm is the 'importance value' of the
% voxel. In the context of a neural network classifier, an
% input unit contains signal from a single voxel. The
% activation value of an input unit can be positive or
% negative.
%
%     imp_{i,j} = w_{i,j} * a_{i,j}      [Latex formatting]
%
% This eqn describes the calculation of this importance
% value. The subscript i is used to step across input units
% in the network. The subscript j is used to step across
% output units in the network. Thus, a single input unit has
% an importance value for each category. This importance
% value is defined as the product of two values, the weight
% between input unit i and output unit j, and the average
% activity of input unit i for the training patterns drawn
% from category j. The importance value is meant to reflect
% how effectively this input unit can alter the value of the
% corresponding output unit. The sign of the importance
% value reflects whether activity in this input unit tends
% to drive the activation value of the corresponding output
% unit up or down. If both wij and aij are positive (or both
% are negative), signal in the input unit causes the output
% unit to turn on; it has a positive importance. If the two
% have opposite sign, then the input unit??????s signal acts to
% turn off this output unit, and it has a negative
% importance value. Network classifiers use both positive
% and negative evidence to reduce error on the training set
% (Polyn et al., 2005). For example, if the signal in a
% voxel increases when a face is viewed, then a strong
% positive weight to the face output unit will reduce
% error. A strong negative weight from this input unit to
% the location output unit will also reduce error, by
% turning off the location output when face is the correct
% category.  (Explanation taken from
% http://polyn.com/struct/PolynDissertation.pdf)
%
% NOTE : This algorithm is designed to operate on a cell
% array, RESULTS, containing a separate individual 'results'
% structure. So if you've run the same cross-validation
% multiple times with a backprop classifier, each time
% initialized randomly, then you can average over these
% multiple attempts to get a better estimate of each voxel's
% contribution. If you feed in a single 'results' structure,
% it will be turned into a 1-cell cell array.
%
% In the end we can get two type importance
% maps one with both positive and negative canonicals and
% the other with only positive canonicals, depending on the
% 'type_canon' argument.
%
% VARARGIN :
% ---------
% TYPE_CANON:
% You can give an option to what type of importance maps you
% want. The default is both positive and negative canonical,
% but if you need only the positive then you can give this
% option in varargin as 'TYPE_CANON' as 'pos'.
%
% IMPMAP_NAME:
% An importance map a is cell array with either both the
% importance maps (all and positives) or only positives.This
% cell array will be added to the subj structure, as new
% patterns using the default group name : IMPMAP_NAME. So by
% default this names are going to be IMPMAP_ALL_1,
% IMPMAP_ALL_2.,etc.. based on the no. of iterations of N-1
% you have in your subj structure. But you can your use own
% name by using varargin 'impmap_name'.
%
% WHAT DOES THE STRUCTURE LOOK LIKE:
%
% The importance map will be an average over all the times
% you run your cross validation script. Thus if you have 8
% runs and 3 conditions and you run the cross validation 50
% times in the end you will end up averaging over all the 50
% iterations, and you will get 8 new maps holding a matrix
% of nVoxels (per N-1 iteration since it different for each
% iteration due to the no peeking anova.) by nConds (3 in
% this case).


% setting defaults
defaults.type_canon = 'all';
defaults.impmap_name = 'impmap';
args = propval(varargin,defaults);


% if RESULTS is not a cell array, then convert into a 1-cell
% cell array
if ~iscell(results)
    temp = results;
    results=[];
    results{1} = temp;
end

% number of times the cross_validation function has been run
nTimes = length(results);
% number of iterations of the N-1
nRuns = length(results{1}.iterations);

sanity_check(results,nTimes,nRuns,args);

% now getting the patterns to average
for i=1:nRuns
    % using RESULTS{1} on the assumption that all of the
    % RESULTS cells will be using the same object names
    pats_name = results{1}.iterations(i).created.patname;
    % since i need this for masked by field
    mask_name{i} = results{1}.iterations(i).created.maskname;
    sel_name  = results{1}.iterations(i).created.selname;
    reg_name  = results{1}.iterations(i).created.regsname;
    
    % really useful fuction Greg!
    masked_pats = get_masked_pattern(subj,pats_name,mask_name{i});
    selectors   = get_mat(subj,'selector',sel_name);
    regressors  = get_mat(subj,'regressors',reg_name);
    
    % this picks out only the training patterns.
    %selectors ==2 <-- changed from 1 to select testing patterns only by AG
    %10/13/09
    %NOTE: changed BACK to selectors == 1 by TIB, 9/2/2014 - using the
    %testing patterns needs work for a scenario where the testing set is
    %only 1 trial. Regardless, the benefit of using the testing values is not clear over the training set per Polyn, Rissman, and others
    
    pats_to_avg  = masked_pats(:,(selectors == 1));
    targs_to_avg = regressors(:,(selectors == 1));
    nConds       = size(regressors,1);
    
    % average the patterns per category.
    for j=1:nConds
        curCols  = find(targs_to_avg(j,:)==1);
        curData  = pats_to_avg(:,curCols);
        nVox_mean = mean(curData,2);
        
        %         if strcmp(args.type_canon,'pos')
        %             nVox_mean(find(nVox_mean < 0)) = 0;
        %         end
        iteration{i}.canonical(:,j) = nVox_mean;
        
        
        for it = 1:length(S.impType)
            %for cds = 1:length(S.regNames)
            importance.(S.impType{it}).(S.regNames{j}).map= zeros(size(iteration{1}.canonical(:,j)));
            
            %     importance{i}.pos.map= zeros(size(iteration{i}.canonical));
            %     importance{i}.neg.map= zeros(size(iteration{i}.canonical));
            %     importance{i}.both.map= zeros(size(iteration{i}.canonical));
            %end
        end
    end
end



% this step collapses across every condition per iteration
% per nTimes
for i=1:nTimes
    %for j=1:nRuns
    for cds = 1:length(S.regNames)
        
        if nRuns==1
            for j=1:nRuns
                curWts_h{cds}(:,j) = (results_IW{j}.iterations(1).scratchpad.net.IW{1}(cds,:)');
                curCanonical_h{cds}(:,j) = iteration{j}.canonical(:,cds);
            end
            curWts{cds} = curWts_h{cds};
            curCanonical{cds} = curCanonical_h{cds};
        else
            for j=1:nRuns % these are the weights of the voxels per iterations.
                curWts_h{cds}(:,j) = (results_IW{j}.iterations(1).scratchpad.net.IW{1}(cds,:)');
                curCanonical_h{cds}(:,j) = iteration{j}.canonical(:,cds);
            end
            curWts{cds} = mean(curWts_h{cds}')';
            curCanonical{cds} = mean(curCanonical_h{cds}')';
        end
        
        
        % these are the weights of the voxels per iterations.
        %for k=1:nConds
        % for all the canonical
        % these are the weights per nTimes per nRuns per cat
        %curWts{cds}=results_IW{j}.iterations(1).scratchpad.net.IW{1}(cds,:)'; % JR: get weights for output unit 1
        %curWts{2}=results_IW{i}.iterations(j).scratchpad.net.IW{1}(2,:)'; % JR: get weights for output unit 2
        
        
        
        %curCanonical{cds} = iteration{j}.canonical(:,cds); % JR: get mean activity for output unit 1
        %curCanonical{2} = iteration{j}.canonical(:,2); % JR: get mean activity for output unit 2
        
        
        %%      compile mean activity (act), weights (wts), and sort by sign
        %information
        Vox.pos.(S.regNames{cds}).act = find(curCanonical{cds}>0);
        Vox.pos.(S.regNames{cds}).wts = find(curWts{cds}>0);
        Vox.pos.(S.regNames{cds}).sameSign = intersect(Vox.pos.(S.regNames{cds}).wts , Vox.pos.(S.regNames{cds}).act);
        
        Vox.neg.(S.regNames{cds}).act = find(curCanonical{cds}<=0);
        Vox.neg.(S.regNames{cds}).wts = find(curWts{cds}<=0);
        Vox.neg.(S.regNames{cds}).sameSign = intersect(Vox.neg.(S.regNames{cds}).wts , Vox.neg.(S.regNames{cds}).act);
        
        %Vox.both.(S.regNames{cds}).act = find(curCanonical{cds});
        %Vox.both.(S.regNames{cds}).wts = find(curWts{cds});
        %Vox.both.(S.regNames{cds}).sameSign = intersect(Vox.both.(S.regNames{cds}).wts , Vox.both.(S.regNames{cds}).act);
        
        Vox.both.(S.regNames{cds}).act = 1: length(curCanonical{cds});
        Vox.both.(S.regNames{cds}).wts = 1: length(curWts{cds});
        Vox.both.(S.regNames{cds}).sameSign = intersect(Vox.both.(S.regNames{cds}).wts , Vox.both.(S.regNames{cds}).act);
        
        curImp.pos.(S.regNames{cds}) = zeros(size(curWts{cds}));
        curImp.neg.(S.regNames{cds}) = zeros(size(curWts{cds}));
        curImp.both.(S.regNames{cds}) = zeros(size(curWts{cds}));
        
        if S.generateBetaMaps %output is raw weights maps. We multiply all values by 1000000 to ensure maps are all on same scale for visualization
            curImp.pos.(S.regNames{cds})(Vox.pos.(S.regNames{cds}).sameSign)  = curWts{cds}(Vox.pos.(S.regNames{cds}).sameSign) * 1000000;
            curImp.neg.(S.regNames{cds})(Vox.neg.(S.regNames{cds}).sameSign) = curWts{cds}(Vox.neg.(S.regNames{cds}).sameSign) * 1000000;
            
            %"both" is the weights where where wts & acts
            % are the same sign (so +wt&+act + -wt&-act)
            curImp.both.(S.regNames{cds})(Vox.both.(S.regNames{cds}).sameSign) = curImp.pos.(S.regNames{cds}) + curImp.neg.(S.regNames{cds});
            
            %"raw" is all weights
            curImp.raw.(S.regNames{cds}) = curWts{cds} * 1000000;
        else %output is true "importance maps" as defined above - wt*act yields voxels with positive importance (+*+ and -*-) as well as negative importance (+*- or -*+). **TIB: scaled here by 1mil like AG did for the betas - this makes maps more displayable + overlayable
            curImp.pos.(S.regNames{cds})(Vox.pos.(S.regNames{cds}).sameSign)  = (curWts{cds}(Vox.pos.(S.regNames{cds}).sameSign) .* curCanonical{cds}(Vox.pos.(S.regNames{cds}).sameSign))*1000000;
            curImp.neg.(S.regNames{cds})(Vox.neg.(S.regNames{cds}).sameSign) = (curWts{cds}(Vox.neg.(S.regNames{cds}).sameSign) .* curCanonical{cds}(Vox.neg.(S.regNames{cds}).sameSign))*1000000;
            curImp.both.(S.regNames{cds})(Vox.both.(S.regNames{cds}).sameSign) = (curImp.pos.(S.regNames{cds}) + curImp.neg.(S.regNames{cds}))*1000000;
            
            curImp.raw.(S.regNames{cds}) = (curWts{cds}.*curCanonical{cds})*1000000;
        end
        
        
        
        %                 for it = 1:length(S.impType)
        %                 % identify voxels with positive activity and positive weights
        %
        %                     if strcmp('pos', S.impType{s})
        %                     Vox.(S.impType{it}).(S.regNames{cds}).act = find(curCanonical{cds}>0);
        %                     Vox.(S.impType{it}).(S.regNames{cds}).wts = find(curWts{cds}>0);
        %                     Vox.(S.impType{it}).(S.regNames{cds}).sameSign = intersect( Vox.(S.impType{it}).(S.regNames{cds}).act, Vox.(S.impType{it}).(S.regNames{cds}).wts);
        %                     elseif strcmp('neg', S.impType{s})
        %                     Vox.(S.impType{it}).(S.regNames{cds}).act = find(curCanonical{cds}>0);
        %                     Vox.(S.impType{it}).(S.regNames{cds}).wts = find(curWts{cds}>0);
        %                     Vox.(S.impType{it}).(S.regNames{cds}).sameSign = intersect( Vox.(S.impType{it}).(S.regNames{cds}).act, Vox.(S.impType{it}).(S.regNames{cds}).wts);
        %                     elseif strcmp('both', S.impType{s})
        %                     Vox.(S.impType{it}).(S.regNames{cds}).act = find(curCanonical{cds});
        %                     Vox.(S.impType{it}).(S.regNames{cds}).wts = find(curWts{cds});
        %                     Vox.(S.impType{it}).(S.regNames{cds}).sameSign = intersect( Vox.(S.impType{it}).(S.regNames{cds}).act, Vox.(S.impType{it}).(S.regNames{cds}).wts);
        %                     end
        %
        %                     FinalImp.(S.impType{it}).(S.regNames{cds}) = zeros(size(Vox.(S.impType{it}).(S.regNames{cds}).wts ));
        %
        %                 end
        
        
        
        %                 pos_act_vox{1} = find(curCanonical{1}>0);
        %                 pos_wts_vox{1} = find(curWts{1}>0);
        %                 both_pos{1} = intersect(pos_act_vox{1},pos_wts_vox{1});
        %                 neg_act_vox{1} = find(curCanonical{1}<0);
        %                 neg_wts_vox{1} = find(curWts{1}<0);
        %                 both_neg{1} = intersect(neg_act_vox{1},neg_wts_vox{1});
        %
        %                 pos_act_vox{2} = find(curCanonical{2}>0);
        %                 pos_wts_vox{2} = find(curWts{2}>0);
        %                 both_pos{2} = intersect(pos_act_vox{2},pos_wts_vox{2});
        %                 neg_act_vox{2} = find(curCanonical{2}<0);
        %                 neg_wts_vox{2} = find(curWts{2}<0);
        %                 both_neg{2} = intersect(neg_act_vox{2},neg_wts_vox{2});
        
        %                 curImp{1} = zeros(size(curWts{1}));
        %                 curImp{2} = zeros(size(curWts{2}));
        
        %                 FinalImp.pos{1} = zeros(size(curWts{1}));
        %                 FinalImp.pos{2} = zeros(size(curWts{2}));
        %
        %                 FinalImp.neg{1} = zeros(size(curWts{1}));
        %                 FinalImp.neg{2} = zeros(size(curWts{2}));
        %
        %                 FinalImp.both{1} = zeros(size(curWts{1}));
        %                 FinalImp.both{2} = zeros(size(curWts{2}));
        
        %                 curImp{1}(both_pos{1})= curWts{1}(both_pos{1}) .* curCanonical{1}(both_pos{1});
        %                 curImp{1}(both_neg{1})= -1* curWts{1}(both_neg{1}) .* curCanonical{1}(both_neg{1});
        %
        %                 curImp{2}(both_pos{2})= curWts{2}(both_pos{2}) .* curCanonical{2}(both_pos{2});
        %                 curImp{2}(both_neg{2})= -1* curWts{2}(both_neg{2}) .* curCanonical{2}(both_neg{2});
        %
        %
        %
        %                 FinalImp.pos{1}(both_pos{1}) = curImp{1}(both_pos{1});
        %                 FinalImp.pos{2}(both_pos{2}) = curImp{2}(both_pos{2});
        %
        %                 FinalImp.neg{1}(both_neg{1}) = curImp{1}(both_neg{1});
        %                 FinalImp.neg{2}(both_neg{2}) = curImp{2}(both_neg{2});
        %
        %                 FinalImp.both{1} = curImp{1};
        %                 FinalImp.both{2} = curImp{2};
        
        %
        %         % adding it to the old importance map of the same n-1
        %         % iteration and for that condition. (here we don't
        %         % separate it per nTimes of the cross validation run)
        %         importance{j}.map(:,1) = importance{j}.map( :,1) + FinalImp{1};
        %         importance{j}.map(:,2) = importance{j}.map( :,2) + FinalImp{2};
        
        
        importance.pos.(S.regNames{cds}).map = importance.pos.(S.regNames{cds}).map + curImp.pos.(S.regNames{cds});
        importance.neg.(S.regNames{cds}).map =  importance.neg.(S.regNames{cds}).map + curImp.neg.(S.regNames{cds});
        importance.both.(S.regNames{cds}).map = importance.both.(S.regNames{cds}).map + curImp.both.(S.regNames{cds});
        
        importance.raw.(S.regNames{cds}).map = importance.raw.(S.regNames{cds}).map + curImp.raw.(S.regNames{cds});
        %importance.raw(S.regNames{cds}.map = importance.raa(S.regNames{cds}).map
        
        %
        %                 importance{j}.pos.map(:,1) = importance{j}.pos.map(:,1) + FinalImp.pos{1};
        %                 importance{j}.pos.map(:,2) =  importance{j}.pos.map(:,2) + FinalImp.pos{2};
        %                 importance{j}.neg.map(:,1) = importance{j}.neg.map(:,1) + FinalImp.neg{1};
        %                 importance{j}.neg.map(:,2) = importance{j}.neg.map(:,2) + FinalImp.neg{2};
        %                 importance{j}.both.map(:,1) = importance{j}.both.map(:,1) + FinalImp.both{1};
        %                 importance{j}.both.map(:,2) = importance{j}.both.map(:,2) + FinalImp.both{2};
    end
end
%end

threshold_maps = 0;%0 = save all voxels to map, 1 = save n voxels to map

% we divide by the no. of times we add each block to
% get an average importance map per iteration of n-1 per condition
% and put everything into the subj structure.
for it = 1:length(S.impType)
    for cds = 1:length(S.regNames)
        %for j=1:nRuns
        if threshold_maps == 1
            %custom code to create importance maps thresholded by n most important
            %voxels
            n = 500; %number of voxels desired
            importance_tmp = abs(importance.(S.impType{it}).(S.regNames{cds}).map);%convert importance values to abs values
            [srtd,idx] = sort(importance_tmp(:),'descend');%sort and maintain idx of original position
            idxd = idx(1:n); %idx of desired voxels to keep in importance map
            empty_map = zeros(size(importance.(S.impType{it}).(S.regNames{cds}).map));%create importance vector of zeros
            empty_map(idxd) = importance.(S.impType{it}).(S.regNames{cds}).map(idxd);%assign values from original importance map @ idxd to new vector
            importance.(S.impType{it}).(S.regNames{cds}).map = empty_map; %and then assign those values to the desired structure
            %sanity checks...
            a = sum(importance.(S.impType{it}).(S.regNames{cds}).map(:,1)~=0); %count non-zero voxels in map vector
            if a > n
                error('It looks like the thresholding did not work!');
            elseif a < n
                warning('FYI, you have fewer important voxels than your threshold');
            end
        
        else
            importance.(S.impType{it}).(S.regNames{cds}).map= importance.(S.impType{it}).(S.regNames{cds}).map;%original code, no voxel threshold
        end
        
        patname = strcat(args.impmap_name, '_', S.impType{it}, '_', S.regNames{cds} );
        subj = init_object(subj,'pattern',patname);
        subj = set_mat(subj,'pattern',patname,importance.(S.impType{it}).(S.regNames{cds}).map);
        subj = set_objfield(subj,'pattern', patname ,'masked_by', mask_name{1});
        subj = set_objfield(subj,'pattern', patname ,'group_name',args.impmap_name);
        created.type_canonicals = args.type_canon;
        subj = add_created(subj,'pattern',patname,created);
    end
end
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = sanity_check(results,nTimes,nRuns,args)

% checking if all the names are ok.
pats_name = results{1}.iterations(1).created.patname;
mask_name = results{1}.iterations(1).created.maskname;
sel_name  = results{1}.iterations(1).created.selname;
reg_name  = results{1}.iterations(1).created.regsname;

if nTimes > 1
    for i=2:nTimes
        newpats_name = results{i}.iterations(1).created.patname;
        newmask_name = results{i}.iterations(1).created.maskname;
        newsel_name  = results{i}.iterations(1).created.selname;
        newreg_name  = results{i}.iterations(1).created.regsname;
        
        if ~strcmp( pats_name,newpats_name) || ~strcmp( mask_name,newmask_name)  ||  ~strcmp( sel_name,newsel_name)|| ~strcmp(reg_name,newreg_name)
            error('The names of your patterns, masks , selectors and regressors matrices must be the same for iterations');
        end
    end
end

for i=2:nTimes
    if ~isequal(nRuns,length(results{i}.iterations))
        error('The no. of runs in each of the cross validations results structures should be the same');
    end
end

% this is to check if it's a backpropagation classifier and
% if so then check if there are hidden units;
% for i=1:nTimes
%     if ~strcmp(results{i}.iterations(1).scratchpad.class_args.train_funct_name,'train_bp')
%         error('This only works with the backpropagation classifier');
%     end
%
%     if ~isequal(results{i}.iterations(1).scratchpad.class_args.nHidden,0)
%         error('No. of nHidden units must be 0');
%     end
% end

if ~(strcmp(args.type_canon,'pos') || strcmp(args.type_canon,'all'))
    error('type_canon can be only ''all'' or ''pos''');
end
%end
