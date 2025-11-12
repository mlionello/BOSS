function [masks, params, downsampleRoi2] = loadCheckpoint(output_folder, roi2_pattern, opts)
    global deriv_path atlas_path cort_labels subcort_labels ...
        cort_atlas_path subcort_atlas_path;
    downsampleRoi2 = false;

    output_folder = dir([output_folder '*/']).folder;
    
    % Load masks and params from checkpoint
    masks_file = fullfile(output_folder, 'masks.mat');
    params_file = fullfile(output_folder, 'params.mat');

    if exist(masks_file, 'file') && exist(params_file, 'file')
        load(masks_file, 'masks');
        load(params_file, 'params');
        fprintf("Loaded masks and params from checkpoint.\n");
    else
        error("ERROR: Checkpoint files not found in %s", output_folder);
    end

    if ~isempty(roi2_pattern) && ~strcmp(params.roi2_pattern, "allbrain")
        error('ERROR: Cannot subsample results if original analysis is not allbrain');
    elseif ~isempty(roi2_pattern) && isempty(opts.suffix)
        error('ERROR: when performing subsampling add a suffix');
    elseif ~isempty(roi2_pattern)
        opts.searchlight = 1;
        params.roi2_pattern = roi2_pattern;
        downsampleRoi2 = true;
        masks4pattern2 = get_masks(atlas_path, deriv_path, [], ...
            cort_labels, roi2_pattern, subcort_labels, ...
            cort_atlas_path, subcort_atlas_path, opts);
        masks.roi2_mask = masks4pattern2.roi2_mask;
    end
end
