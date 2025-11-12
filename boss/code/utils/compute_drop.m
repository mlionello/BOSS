function [corr_res_perm, corr_full_tmp] = compute_drop(paths_nifti, suffix, model, output_folder, center_indices_passing, params)

    if ~exist(fullfile(output_folder, "corr_res_perm" + suffix + ".mat"), 'file')
        [~, corr_full_tmp, ~, corr_res_tmp, ~, corr_res_perm] = infer(paths_nifti, model, output_folder, ...
            0, 10, center_indices_passing, params, "drop");

        save(fullfile(output_folder, "corr_res_perm" + suffix), 'corr_res_perm', 'corr_res_tmp', 'corr_full_tmp', '-v7.3');
    else
        load(fullfile(output_folder, "corr_res_perm" + suffix + ".mat"), 'corr_res_perm', 'corr_full_tmp');
    end

end
