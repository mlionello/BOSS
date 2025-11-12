function draw_ttest_results_on_subj_AVGcanocorr(AVG_targetsource, ...
    targetRegionMask, sourceRegionMask, template, info_nifti_nifti, ...
    output_folder, target_canvars_mean_activity, source_canvars_mean_activity, pval)

    targetMask_indices = find(targetRegionMask);
    sourceMask_indices = find(sourceRegionMask);
    
    values2nifti = cat(1, zscore(AVG_targetsource.A_AVG, 0, 2), zscore(AVG_targetsource.B_AVG, 0, 2));
    values2nifti = values2nifti(:, 1: min(size(values2nifti, 2), 100));
    indices2nifti = cat(1, targetMask_indices, sourceMask_indices);
    draw_and_save_nifti(template, fullfile(output_folder, 'avg_subj_comp'), ...
        indices2nifti, info_nifti_nifti, values2nifti);
    
    target_canvars_mean_activity_mat = cat(1, target_canvars_mean_activity{:});
    source_canvars_mean_activity_mat = cat(1, source_canvars_mean_activity{:});
    
    p_passing_indices = find(pval(1,:)<0.05);
    values2nifti = cell(1,2);
    for p_i = 1 : numel(p_passing_indices)
        p_index = p_passing_indices(p_i);
        values2nifti{1} = cat(1, ...
            AVG_targetsource.A_AVG(:, p_index), ...
            AVG_targetsource.B_AVG(:, p_index));
    
        values2nifti{2} = squeeze(cat(2, ...
            zscore(mean( target_canvars_mean_activity_mat(:, p_index)*AVG_targetsource.A_AVG(:,p_index)', 1), 0, 2), ...
            zscore(mean( source_canvars_mean_activity_mat(:, p_index)*AVG_targetsource.B_AVG(:,p_index)', 1), 0, 2)));
    
        draw_and_save_nifti(template, fullfile(output_folder, sprintf('all_subj_comp%d_n%d', p_index, p_i)), ...
            indices2nifti, info_nifti_nifti, values2nifti);
    end
end
