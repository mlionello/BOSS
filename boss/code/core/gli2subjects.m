function [pval_full, corr_full, pval_res, corr_res, ...
    mean_activity, corr_res_perm] = gli2subjects( ...
    sourceRoi, targetRoi, infeat, opts)
arguments
    sourceRoi;
    targetRoi;
    infeat;
    opts.coeff = [];
    opts.vx_shuffling = 0;
    opts.save_activity = 0;
    opts.nb_feat_shuffling = 0;
    opts.analyser = [];
end

    % apply trasnformation to the two nifti regions
    if isempty(opts.analyser)
        source_stats = sourceRoi;
        target_stats = targetRoi;
    else
        [source_stats, target_stats] = opts.analyser(sourceRoi, targetRoi, opts.coeff);
    end

    mean_activity = struct();
    if opts.save_activity
        mean_activity.source = squeeze(mean(source_stats, 1));
        mean_activity.target = squeeze(mean(target_stats, 1));
    end

    % initialiaze permutation schema to be the same for resiudal and
    % fullmodel studies
    permSchema = initPermutationTests(opts.vx_shuffling, size(infeat{1}, 1));

    % calculate correlation and pvalues for the activity between the two
    % regions for full model and when subtracted it (residual) 

    [pval_full, corr_full] = corr_and_maxstat_on_perm(source_stats, target_stats, permSchema);
    corr_res = cell(1, numel(infeat));
    pval_res = cell(1, numel(infeat));
    corr_res_perm = cell(1, numel(infeat));

    for feat_set_i = 1 : numel(infeat)
        [pval_res{feat_set_i}, corr_res{feat_set_i}] = corr_and_maxstat_on_perm(source_stats, target_stats, permSchema, infeat{feat_set_i});
        corr_res_perm{feat_set_i} = permute_feat(source_stats, target_stats, infeat{feat_set_i}, opts.nb_feat_shuffling);
    end

end


function corr_res_perm = permute_feat(x, y, infeat, nb_rep)
    % Perform residual analysis for different subsets of features
    % permuting single features
    num_features = 1; %size(infeat, 2); TODO: FIX THIS!!!
    % corr_res_perm = cell(num_features, nb_rep);
    corr_res_perm = nan(num_features, nb_rep, size(x, 2));
    if nb_rep == 0; return; end

    for i = 1:num_features
        for k = 1:nb_rep
            feat_subset = infeat;
            feat_subset = feat_subset(randperm(size(feat_subset, 1)), :);
            corr_res_perm(i, k, :) = fast_corr( ...
                get_residual(x, feat_subset), ...
                get_residual(y, feat_subset) ...
                );
        end
    end
end