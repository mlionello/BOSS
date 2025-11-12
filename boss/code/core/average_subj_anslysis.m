function average_subj_anslysis(info_nifti, paths_nifti, ...
            masks, feature_names, sublist, ...
            feature_regressors, output_folder, meansubject_path, commonmask_path)

    meansubj_exists = exist(meansubject_path + ".mat", 'file');
    common_mask_exists = exist(commonmask_path + ".mat", 'file');
    
    all_brain_mask = masks.cortAndSubMaks;
    perform_individual_mapping = false; 

    if ~meansubj_exists || ~common_mask_exists
        [nifti_image_avg, r2, common_mask] = get_avgsubj_and_individualR2(paths_nifti, ...
            all_brain_mask, feature_names, feature_regressors, sublist, perform_individual_mapping);

        if ~meansubj_exists
            niftiwrite(single(nifti_image_avg), meansubject_path, 'Compressed', true)
            save(meansubject_path, 'nifti_image_avg', 'r2', '-v7.3')
        end
        if ~common_mask_exists
            save(commonmask_path, 'common_mask', '-v7.3')
        end
    end

    % concatenate mean r2 to individual ones
    if ~perform_individual_mapping
        return
    end

    if ~exist(fullfile(output_folder, 'r2_sub_and_avg.mat'), 'file')
        load(fullfile(meansubject_path + ".mat"), 'nifti_image_avg', 'r2');
    else
        return
    end
    
    msg = compose("avg r2 vxwise...");
    fprintf(msg);
    cr0 = repmat('\b', 1, strlength(msg));
    for feat_type = feature_names
        feat_name = feat_type{1};

        tmp_cat = cat(2, r2.(feat_name){:,1});
        r2_meanch1 = mean(tmp_cat, 2);
        tmp_cat = cat(2, r2.(feat_name){:,2});
        r2_meanch2 = mean(tmp_cat, 2);

        r2.(feat_name) = [{r2_meanch1, r2_meanch2}; r2.(feat_name)];
    end

    msg = compose("computing average subject...");
    fprintf(cr0 + msg);
    cr0 = repmat('\b', 1, strlength(msg));

    for ch = 1:2
        r2AVG{ch} = r2_features_img(nifti_image_avg, masks.cortAndSubMaks, ...
            feature_names, cellfun(@(x) x(:,:,ch), feature_regressors, 'UniformOutput', false));
    end
    r2AVG = cellArrayToStruct(r2AVG);

    msg = compose("writing to nifti files...");
    fprintf(cr0 + msg);
    cr0 = repmat('\b', 1, strlength(msg));

    template = create_template(size(masks.roi1_mask{1}), 1, 1);

    fileID = fopen(fullfile(output_folder, 'r2subj_list.txt'),'w'); fprintf(fileID,'%s\n',string(paths_nifti)); fclose(fileID);
    for feat_type = feature_names
        feat_name = feat_type{1};
        draw_and_save_nifti(template, fullfile(output_folder, "r2_allsubj_" + feat_name + "VsAllbrain"), find(masks.cortAndSubMaks), info_nifti, r2.(feat_name));
        draw_and_save_nifti(template, fullfile(output_folder, "r2_AVGsubj_" + feat_name + "VsAllbrain"), find(masks.cortAndSubMaks), info_nifti, r2AVG.(feat_name));
    end
    save(fullfile(output_folder, 'r2_sub_and_avg.mat'), 'r2', 'r2AVG');

    fprintf(cr0);

end
