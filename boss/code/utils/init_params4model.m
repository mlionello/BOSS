function init_model = init_params4model(output_folder, params, masks, opts)
    global HCP_atlas_path commonmask_path;

    load(commonmask_path, 'common_mask')
    
    init_model.output_folder = output_folder;
    init_model.feature_regressors = params.feature_regressors;
    init_model.feature_names = params.feature_names;
    init_model.sourceRegionMask = cellfun(@(x) x .* common_mask, masks.roi1_mask, 'UniformOutput', false);
    init_model.targetRegionMask = masks.roi2_mask .* common_mask;
    init_model.targetRadius = params.radius;
    init_model.targetPadding = params.padding;
    init_model.targetPropInVx = params.pVxIn;
    if matches(params.target_model, "IterateAllAtlas_CanonCorr")
        init_model.atlas_img = single(niftiread(HCP_atlas_path));
        init_model.source_label_id = single(masks.roi1_label_id);
    end
    init_model.changed_features = params.changed_features;
end
