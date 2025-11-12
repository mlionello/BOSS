function [masks, params, downsampleRoi2, opts] = initializeMasksAndParams(roi1_pattern, roi2_pattern, opts, output_folder)
    global atlas_path deriv_path cort_labels subcort_labels HCP_labels ...
        cort_atlas_path subcort_atlas_path HCP_atlas_path feature_list;

    channel_changed = false;
    radius_changed = false;
    padding_changed = false;

    downsampleRoi2 = false; % Default value
    if ~isempty(opts.checkpoint)
        if ~isempty(roi1_pattern)
            error("ERROR: Checkpoint provided with additional ROI patterns.");
        end
        % Load checkpoint data
        % TODO: if opts.suffix the folder must loaded from the original dir
        [masks, params, downsampleRoi2] = loadCheckpoint(output_folder, roi2_pattern, opts);
        params.changed_features = false;

        channel_changed = ~isequal(params.channel, opts.channel);
        radius_changed = ~isequal(params.radius, opts.radius);
        padding_changed = ~isequal(params.padding, opts.padding);
        if opts.update_params
            if all([isempty(opts.feature_list), ~radius_changed, ~padding_changed])
                error("opts.feature_list is true but htere are no params to update");
            end
            if ~isempty(opts.feature_list)
                params = get_music_features(params, opts.feature_list, 24);
                params.changed_features = 1;
                fprintf("updating feature_list: if the model has already been fit, " + ...
                    "this should be done only when the training does not depend on the features!\n");
            end
            if length(opts.radius) ~= length(opts.padding)
                error("radius and padding have different number of elements");
            end
            if radius_changed || padding_changed
                params.radius = opts.radius;
                params.padding = opts.padding;
            end
            if channel_changed
                params.channel = opts.channel;
            end
        end
        if any([~isempty(opts.feature_list), radius_changed, padding_changed]) && ~opts.update_params
            error("you give params to be updated but opts.feature_list is false");
        end
    else
        % Initialize parameters and compute masks
        if any([opts.update_params, ~isempty(opts.feature_list), radius_changed, padding_changed, channel_changed])
            error("if no checkpoint is given, model cannot be updated for the anlaysis.");
        end

        params = initializeParams(opts, roi1_pattern, roi2_pattern, output_folder);
        params.changed_features = false;
        masks = get_masks(atlas_path, deriv_path, roi1_pattern, ...
            cort_labels, roi2_pattern, subcort_labels, HCP_labels, ...
            cort_atlas_path, subcort_atlas_path, HCP_atlas_path, opts);

        fprintf("Computed masks using get_masks.\n");
        saveMasksAndParams(output_folder, masks, params);
    end

    if ~isempty(roi2_pattern) && matches(roi2_pattern, "allbrain")
        opts.searchlight = 1;
    end
end
