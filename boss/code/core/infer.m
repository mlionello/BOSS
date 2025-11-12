function [pval_full, corr_full, pval_res, corr_res, ...
            mean_activity, corr_res_perm] = infer( ...
            paths_nifti, model, output_folder, nb_perm_vx, nb_perm_feat, indices, params, suffix, opts)
    pval_full = {};
    corr_full = {};
    pval_res = {};
    corr_res = {};
    corr_res_perm = {};
    mean_activity = nan;

    tic;
    cr1 = ""; cr0 = "";
    for sub_i = 1 : length(paths_nifti)
        if exist(fullfile(output_folder, "corr_pval_subj_" + string(sub_i)) + suffix + '.mat', "file")
            continue
        end
        msg = compose("%s subject n. %d/%d ", string(duration(0, 0, toc)), sub_i, length(paths_nifti));
        fprintf("" + cr1 + cr0 + msg);
        cr0 = repmat('\b', 1, strlength(msg));

        pathnifti = paths_nifti{sub_i};
        msg = compose(" ...loading image");
        fprintf( msg);
        cr1 = repmat('\b', 1, strlength(msg));
        nifti_image = niftiread(pathnifti);

        msg = compose(" ...applying cca coeff");
        fprintf(cr1 + msg);
        cr1 = repmat('\b', 1, strlength(msg));
        [pval_full2save, ...
            corr_full2save, ...
            pval_res2save, ...
            corr_res2save, ...
            corr_res_perm2save] = model.predict(...
                nifti_image, nb_perm_vx, nb_perm_feat, ...
                'out_dir', output_folder, ...
                'radius', params.radius, ...
                'channel', params.channel, ...
                'sub_indices', indices);

        pval_full(sub_i, :, :, :) = pval_full2save;
        corr_full(sub_i, :, :, :) = corr_full2save;
        pval_res(sub_i, :, :, :) = pval_res2save;
        corr_res(sub_i, :, :, :) = corr_res2save;
        corr_res_perm(sub_i, :, :, :) = corr_res_perm2save;
    
        save(fullfile(output_folder, "corr_pval_subj_" + string(sub_i)) + suffix, 'corr_res2save', ...
            'corr_full2save', 'pval_res2save', 'pval_full2save', 'corr_res_perm2save')
    end

    cr0 = "";
    for sub_i = 1: length(paths_nifti)
        msg = compose("loading individual inferences %d/%d ...", sub_i, length(paths_nifti));
        fprintf("" + cr0 + msg);
        cr0 = repmat('\b', 1, strlength(msg));
        %if isempty(corr_full) || isempty(corr_full(sub_i))
            load(fullfile(output_folder, "corr_pval_subj_" + string(sub_i)) + suffix + ".mat");
            pval_full(sub_i, :, :, :) = pval_full2save;
            corr_full(sub_i, :, :, :) = corr_full2save;
            pval_res(sub_i, :, :, :) = pval_res2save;
            corr_res(sub_i, :, :, :) = corr_res2save;
            if ismember('corr_res_perm2save', who)
                corr_res_perm(sub_i, :, :, :) = corr_res_perm2save;
            end
        %end
    end
end
