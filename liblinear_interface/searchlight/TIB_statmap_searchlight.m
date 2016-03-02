function [subj] = TIB_statmap_searchlight(subj,data_patname,regsname,selname,new_map_patname,extra_arg)

%modified by TIB to include additional memory-saving features. These remove
%potentially useful data (scratch.class_scratch) and should really only be
%employed in memory-extreme circumstances (e.g. wholebrain high resolution
%data)

% [SUBJ] = STATMAP_SEARCHLIGHT( SUBJ, ...
%   DATA_PATNAME,REGSNAME,SELNAME,NEW_MAP_PATNAME,EXTRA_ARG)
%
% OBJ_FUNCT (required, default = ''),
% e.g. 'statmap_classify'
%
% ADJ_LIST (required) - as produced by CREATE_ADJ_LIST.M. If
% you want to run this univariately, feed in an ADJ_LIST
% where every cell is empty.
%
% IGNORE_EMPTY_ADJ_LIST (optional, default = false). By
% default, will not allow the ADJ_list to be empty (since
% this is designed for multivariate algorithms, you want the
% searchlights to have more than one voxel in them). Set
% this to false to allow single-voxel searchlights.
%
% ADD_CENTER_VOXEL (optional, default = false). By default,
% don't add the center voxel to its own sphere, since the
% default for ADJ_SPHERE.M is to include the center
% voxel.
%
% SCRATCH (optional, default = []). Gets passed into the
% OBJ_FUNCT as a third argument.

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


defaults.cur_iteration = [];
defaults.obj_funct = '';
defaults.adj_list = [];
defaults.add_center_voxel = true;
defaults.scratch = [];
defaults.ignore_empty_adj_list = false;
defaults.memsave = 0;%by default we allow scratch.class_scratch to grow across voxels
defaults.parallel = 0;%by default we run searchlights serially

args = propval(extra_arg, defaults);

scratch = args.scratch;

% objective function is required
if isempty(args.obj_funct)
    error('Objective function required');
end
[funct_h funct_n] = arg2funct(args.obj_funct);

% adj list is required
if ~isfield(args,'adj_list')
    error('ADJ_LIST required');
end

if size(args.adj_list,2)<2
    % because this is a multivariate algorithm, require there to
    % be some spheres that have more than one voxel in them
    warning('Feeding in empty adj cells');
end

pat  = get_mat(subj,'pattern',data_patname);
regs = get_mat(subj,'regressors',regsname);
sel  = get_mat(subj,'selector',selname);

sanity_check(pat,regs,sel,args);

xval1 = find(sel==1);
pat1  = pat(:,xval1);
regs1 = regs(:,xval1);

xval3_idx = find(sel==3);
% if the user has any 3s in their xval selector, create the
% pat3 and regs3 matrices, e.g. for use in
% statmap_classify.m.
if length(xval3_idx)
    
    % the actual subset of pat3 that gets sent in will be
    % different for every sphere, so we're not going to put
    % this in the scratch yet
    pat3 = pat(:,xval3_idx);
    % this can go in the scratch, since it'll be the same
    % for each sphere
    scratch.regs3 = regs(:,xval3_idx);
end

[nVox nTRs] = size(pat1);
% this is useful for the objective function to know, in case it wants
% to preallocate a matrix for storing things for all spheres
scratch.nVox = nVox;

map = nan(nVox,1);



% loop over the spheres (or whatever adjacency neighborhood you're
% using)
if args.parallel == 1
    if args.memsave < 3
        multilog = zeros([nVox,size(scratch.regs3)]); %initialize matrix to store logits from the parfor scenario
    end
    
    cluster = parcluster('local');
    tmpdirforpool = '/hsgs/projects/awagner/thackery/parfortest';
    mkdir(tmpdirforpool);
    cluster.JobStorageLocation = tmpdirforpool;
    cluster.NumWorkers = str2num(getenv('NSLOTS'));
    
    matlabpool(cluster);
    
    cluster;
    
    parfor v = 1:nVox
        
        cur_adj_list = args.adj_list(v,:);
        % this looks like madness, but it's just saying 'throw
        % away all non-zero values from me', i.e. cur_adj_list =
        % cur_adj_list(find(cur_adj_list));
        cur_adj_list = cur_adj_list(cur_adj_list~=0);
        
        progress(v,nVox);
        
        if args.add_center_voxel
            % add the current (center) voxel to its own sphere
            vox_idx = [v cur_adj_list];
        else
            vox_idx = cur_adj_list;
        end
        
        % get the spherical subset of voxel data
        sphere = pat1(vox_idx,:);
        
        scratch1 = scratch;
        
        scratch1.pat3 = pat3(vox_idx,:);
        % store the current voxel # being operated on
        scratch1.v_counter = v;
        
        % calculate value of this sphere, using the objective
        % function, e.g. STATMAP_SIMST_LOGIC.M or
        % STATMAP_CLASSIFY.M
        %
        % it's important to catch the SCRATCH, since some
        % functions may want to cache things for themselves in
        % there, to avoid recomputing them each time
        
        %[map(v) subj.acts scratch] = funct_h(sphere, regs1, scratch);
        [map(v) scratch1] = funct_h(sphere, regs1, scratch1);
        
        if args.memsave < 3
            multilog(v,:,:) = scratch1.multi_acts(v,:,:);
        end
        
        %added by TIB to optionally clear class_scratch. This info is potentially
        %useful for additional data exploration and/or debugging, but NOT needed to simply generate searchligh perf maps
        if args.memsave == 3
            scratch1 = rmfield(scratch1, 'class_scratch');
        end
        
    end
    
    if args.memsave < 3
        scratch.multi_acts = multilog;%append the multilog to scratch now to write out prior to clearing everything
    end
    
    matlabpool close;
    
else
    % loop over the spheres (or whatever adjacency neighborhood you're
    % using)
    for v = 1:nVox
        
        cur_adj_list = args.adj_list(v,:);
        % this looks like madness, but it's just saying 'throw
        % away all non-zero values from me', i.e. cur_adj_list =
        % cur_adj_list(find(cur_adj_list));
        cur_adj_list = cur_adj_list(cur_adj_list~=0);
        
        progress(v,nVox);
        
        if args.add_center_voxel
            % add the current (center) voxel to its own sphere
            vox_idx = [v cur_adj_list];
        else
            vox_idx = cur_adj_list;
        end
        
        % get the spherical subset of voxel data
        sphere = pat1(vox_idx,:);
        
        scratch.pat3 = pat3(vox_idx,:);
        % store the current voxel # being operated on
        scratch.v_counter = v;
        
        % calculate value of this sphere, using the objective
        % function, e.g. STATMAP_SIMST_LOGIC.M or
        % STATMAP_CLASSIFY.M
        %
        % it's important to catch the SCRATCH, since some
        % functions may want to cache things for themselves in
        % there, to avoid recomputing them each time
        
        %[map(v) subj.acts scratch] = funct_h(sphere, regs1, scratch);
        [map(v) scratch] = funct_h(sphere, regs1, scratch);
        
        %added by TIB to optionally clear class_scratch. This info is potentially
        %useful for additional data exploration and/or debugging, but NOT needed to simply generate searchligh perf maps
        %if args.memsave == 2
        %    scratch = rmfield(scratch, 'class_scratch');
        %end
        
    end
    
end



if length(find(isnan(map)))
    warning('There are NaNs in the map');
end

masked_by = get_objfield(subj,'pattern',data_patname,'masked_by');
subj = initset_object(subj,'pattern',new_map_patname,map, ...
    'masked_by',masked_by);

hist = sprintf('Created by %s',mfilename());
subj = add_history(subj,'pattern',new_map_patname,hist);

created.function = mfilename();
created.obj_funct_name = funct_n;
created.obj_funct_handle = funct_h;
created.data_patname = data_patname;
created.regsname = regsname;
created.selname = selname;
created.new_map_patname = new_map_patname;
created.extra_arg = extra_arg;
created.extra_arg.adj_list = NaN; % to save memory
created.args = args;
created.args.adj_list = NaN; % to save memory
%if args.memsave == 3
%    created.scratch = scratch1;
%else
    created.scratch = scratch;
%end
created.unused = NaN; % to save memory
subj = add_created(subj,'pattern',new_map_patname,created);

%added by TIB to optionally move pattern information to HDD and out of RAM.
if args.memsave >= 1
    %patpath = ['/hsgs/projects/awagner/thackery/Circmaze_mvpa_files/betas/' subj.header.id]
    patpath = ['/Volumes/group/awagner/wagner/thackery/testo/' subj.header.id];
    %subj.logits = scratch.multi_acts; %add the class*trial logits matrix from this classifier iteration (assumes we are moving patterns to HD on a "run-by-run" iteration basis
    [subj] = TIB_move_pattern_to_hd(subj, new_map_patname, 'subdir', patpath, 'memsave', args.memsave);
    
    %added by TIB to optionally clear class_scratch. This info is potentially
    %useful for additional data exploration and/or debugging, but NOT needed to simply generate searchligh perf maps
    if args.memsave == 2
        %         if args.parallel == 1
        %             scratch1 = rmfield(scratch1, 'class_scratch');
        %         else
        scratch
        %scratch = rmfield(scratch, 'class_scratch');
        %        end
        %    elseif args.memsave == 2
        
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = sanity_check(pat,regs,sel,args)

if size(pat,2) ~= size(regs,2)
    error('Different number of timepoints in pattern and regressors');
end
if size(pat,2) ~= size(sel,2)
    error('Different number of timepoints in pattern and selector');
end



