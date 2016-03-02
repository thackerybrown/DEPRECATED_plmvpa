function [subj] = TB_scramble_regressors(subj,regsname,selname,activesname, new_regsname,varargin)
% TB modification of JR script. In this variant, we add additional
% constraint that shuffling will repeat until there are no adjacent
% repetitions of a class in time. This only makes sense for a study where
% this was true in the true design and we want to mimic that property in
% our scrambled data.

% Scrambles your regressors for sanity-checking
%
% [SUBJ] = SCRAMBLE_REGRESSORS(SUBJ,REGSNAME,SELNAME,NEW_REGSNAME,...)
%
% SUBJ = the subj structure
%
% REGSNAME = the name of the regressors you want to scramble
%
% SELNAME = the selector that you want to constrain your scrambling
% with. Usually you will want to scramble within runs, so that you
% have the same number of timepoints in each run. Therefore,
% reference your 'runs' variable for the selname
%
% NEW_REGSNAME = the name you want to give your new scrambled
% regressors matrix
%
%   xxx - shouldn't this be optional???
%
% IGNORE_1OFN (optional, default = false). If your regressors
% are continuous-valued, contain rest or contain multiple active
% conditions in a timepoint, then you might want to scramble them
% in a more sophisticated way that ensures that their
% characteristics are preserved. By default, you'll get warned if
% your regressors aren't in basic 1-of-n form, unless you set this
% to true.

% License:
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

%add no adjacency rule
defaults.ignore_adjacency = 1; %by default, it's ok if trials of the same class wind up adjacent in time

defaults.ignore_1ofn = false;
args = propval(varargin,defaults);

regs = get_mat(subj,'regressors',regsname);
runs = get_mat(subj,'selector',selname);
actives = get_mat(subj,'selector',activesname);
nonactive_trials = find(actives == 0);
runs(nonactive_trials)=0;


% See comment about IGNORE_1OFN above
[isbool isrest isoveractive] = check_1ofn_regressors(regs);
if ~isbool || isrest || isoveractive
    if ~args.ignore_1ofn
        warn = ['Your regressors aren''t in basic 1-of-n form - ' ...
            'you may want to consider more sophisticated shuffling techniques'];
        warning(warn);
    end
end

% These next lines will shuffle your regressors within each run
if args.ignore_adjacency == 1
    for i = nonzeros(unique(runs))'
        thisrun = find(runs == i);
        regs(:,thisrun) = shuffle(regs(:,thisrun),2);
    end
elseif args.ignore_adjacency == 0 %we want to apply additional constraint where scrambling can't put trials of the same class next to each other in time
    for i = nonzeros(unique(runs))'
        thisrun = find(runs == i);
        
        adj = 1; %initialize state 1, where it is assumed there are adjacent class labels of same class (which we don't want)
        while adj == 1
            regs(:,thisrun) = shuffle(regs(:,thisrun),2);
            %b = shuffle(a,2)
            
            adjpresent = []; %initialize vector of adjacency hits
            %for i = 1:length(b(:,1))
            for i = 1:length(regs(:,1))
                %x = find(b(i,:)==1)
                x = find(regs(i,thisrun)==1)
                y = diff(x(1,:))==1
                if sum(y) > 0 %if there was at least one adjacency, the sum will be > 0, and we flag this
                    adjpresent = [adjpresent 1];
                end
            end
            
            if sum(adjpresent) > 0
                adj = 1
            else
                adj = 0
            end
        end
        
        
        %regs(:,thisrun) = shuffle(regs(:,thisrun),2);
    end
end

subj = duplicate_object(subj,'regressors',regsname,new_regsname);
subj = set_mat(subj,'regressors',new_regsname,regs);

hist = sprintf('Regressors ''%s'' created by scramble_regressors',new_regsname);
subj = add_history(subj,'regressors',new_regsname,hist,true);

created.function = 'scramble_regressors';
created.regsname = regsname;
created.selname = selname;
subj = add_created(subj,'regressors',new_regsname,created);

