function [nifti_image_avg, r2, common_mask] = get_avgsubj_and_individualR2(paths_nifti, ...
    mask, feature_names, feature_regressors, sublist, do_feature_mapping)

    r2 = struct();
    
    cr0="";
    tic;
    for sub_i=1:length(paths_nifti)
        msg = compose("%s subject n. %d/%d ", string(duration(0,0,toc)), sub_i, length(paths_nifti));
        fprintf("" + cr0 + msg);
        cr0 = repmat('\b', 1, strlength(msg));
    
        pathnifti = paths_nifti{sub_i};
        msg = compose(" ...loading image");
        fprintf(msg);
        cr1 = repmat('\b', 1, strlength(msg));
        nifti_image = niftiread(pathnifti);
        nifti_image = zscore(nifti_image, 0, 4);
        fprintf( cr1);
    
        if sub_i == 1
            nifti_image_avg = nifti_image/length(sublist);
            common_mask = any(nifti_image, 4);
        else
            nifti_image_avg = nifti_image_avg + nifti_image/length(sublist);
            common_mask = common_mask .* any(nifti_image, 4);
        end
        if do_feature_mapping
            for feat_type = feature_names
            feat_name = feat_type{1};
                for ch = 1:2
                    r2.(feat_name){sub_i, ch} = r2_features_img(nifti_image, mask, ...
                        feature_names, cellfun(@(x) x(:,:,ch), feature_regressors, 'UniformOutput', false));
                end
            end
        end
    end
    
    if do_feature_mapping
        % Convert cell array to struct of cell arrays
        r2 = cellArrayToStruct(r2);
    end
    
    cr0="\n";
    msg = compose("saving to average epi...");
    fprintf(cr0 + msg);
end