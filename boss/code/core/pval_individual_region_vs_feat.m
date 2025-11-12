function pval_individual_region_vs_feat(numPermutations, feature_set, raw_corr, resCorr, info_nifti, template)
    msg = compose("computing pval...");
    fprintf(cr0 + msg);
    cr0 = repmat('\b', 1, strlength(msg));
    for feat_type = feature_set
        feat_name = feat_type{1};
        [res_p, res_punc, ~] = compute_signflip_test(resCorr.(feat_name), raw_corr, numPermutations);
        pval.(feat_name) = res_p;
        pval_unc.(feat_name) = res_punc;
    end
    passing_comp = find(pval.feat_all(1,:)<0.05);
    indices2nifti = cat(1, targetMask_indices, sourceMask_indices);
    for subj = 1:length(sublist)
        source_volume_cell = mat2cell(B{subj}(:, passing_comp), size(B{subj}, 1), repelem(1, numel(passing_comp)) );
        target_volume_cell = mat2cell(A{subj}(:, passing_comp), size(A{subj}, 1), repelem(1, numel(passing_comp)) );
        concat_volume_cell = cellfun(@(x,y) cat(1, x, y), target_volume_cell, source_volume_cell, 'UniformOutput', false);
        draw_and_save_nifti(template, ...
             fullfile(output_folder, "subj-" + num2str(subj) + "_cca_heatmap_concat"), ...
             indices2nifti, ...
             info_nifti, ...
             concat_volume_cell);
    end
    
    save(fullfile(output_folder, 'pval.mat'),  'pval', 'pval_unc');
    save(fullfile(output_folder, 'info_nifti.mat'),  'info_nifti');
end