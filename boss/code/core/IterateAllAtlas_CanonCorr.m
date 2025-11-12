classdef IterateAllAtlas_CanonCorr < CanonCorr_Base
    properties
        analyser
        mdata
        saveactivity
        atlas_img
        source_label_id
        feat_pca
        featOnTarget

        center_indices
    end
    methods
        function obj = IterateAllAtlas_CanonCorr(opts)
            arguments
                opts (1,1) struct
            end
            opts.targetRegionMask = nan;
            obj = obj@CanonCorr_Base(opts);
            fields = { 'atlas_img', ...
                'analyser', 'saveactivity', 'source_label_id', 'feat_pca'};
            init = struct('analyser', @perform_cca_on_searchlight, ...
                'saveactivity', 0);
            obj = parse_args(obj, fields, init, opts);

            obj.mdata = struct('sourceRegionMask', obj.sourceRegionMask, ...
                'feat_all', obj.feat_all, ...
                'winkler', obj.winkler, ...
                'nb_perm', obj.nb_perm, ...
                'extra', obj.extra);

            % the function given the 3d volume atlas_img, returns a cell
            % array with all the indices for each unique value
            obj.center_indices = IterateAllAtlas_CanonCorr.get_indices_per_roi_from_atlas(obj.atlas_img, obj.source_label_id);
            obj.coeff = cell(numel(obj.center_indices), 1);
        end

        function obj = fit(obj, nifti_image_avg)
            sourceRegion2D = extract2DRois(find(obj.sourceRegionMask), nifti_image_avg);
            
            for i = 1:numel(obj.center_indices)
                if isempty(obj.center_indices{i})
                    continue
                end
                targetRoi2D = extract2DRois( ...
                    obj.center_indices{i}, ...
                    nifti_image_avg);

                featvsTarget = obj.analyser(obj.feat_all, targetRoi2D);
                obj.featOnTarget{i} = featvsTarget.r;

                obj.coeff{i} = obj.analyser(sourceRegion2D, targetRoi2D);
            end
        end

        function [pval_full, corr_full, ...
                pval_res, corr_res, ...
                mean_activity, corr_res_perm] = predict(obj, nifti_image, vx_shuffling, feat_shuffling, indices)
            if isempty(indices)
                indices = 1 : numel(obj.center_indices);
            end
            center_indices_run = obj.center_indices;

            pval_full = cell(numel(indices), 1);
            corr_full = cell(numel(indices), 1);
            pval_res = cell(numel(indices), 1);
            corr_res = cell(numel(indices), 1);
            mean_activity = cell(numel(indices), 1);

            % extract from the nifti the timepoints by voxels images for the given a mask
            sourceRoi = extract2DRois(find(obj.sourceRegionMask), nifti_image);

            feat_all_run = obj.feat_all;
            feat_pca_run = obj.feat_pca;
            saveactivity_run = obj.saveactivity;
            coeff_run = obj.coeff;

            parfor_progress(numel(indices));
            init_time = tic;
            cr = "";
            for i = 1 : numel(indices)
                current_indx = indices(i);
                if isempty(center_indices_run{current_indx})
                    continue
                end
                if isempty(getCurrentTask())
                    msg = compose("%d out of %d; elapsed %s; eta: %s", i, numel(indices), ...
                    duration(0,0,toc(init_time)), ...
                    duration(0,0,toc(init_time)/i*(numel(indices)-i)));
                    fprintf(cr + msg);
                    cr = repmat('\b', 1, strlength(msg));
                else
                    parfor_progress;
                end

                targetRoi = extract2DRois(center_indices_run{current_indx}, nifti_image);

                [pval_full{i}, corr_full{i}, pval_res{i}, corr_res{i}, ...
                    mean_activity{i}, corr_res_perm{i}] = gli2subjects( ...
                    sourceRoi, targetRoi, coeff_run{current_indx}, ...
                    vx_shuffling, struct('feat', feat_all_run, 'pca', feat_pca_run), saveactivity_run, feat_shuffling);
            end
            fprintf(cr);
            parfor_progress(0);
        end

    end
    methods (Static)
        function center_indices = get_indices_per_roi_from_atlas(atlas_img, exclude_inds)
            % Initialize a cell array to store indices for each unique value
            labels = unique(atlas_img);
            nb_rois = numel(labels);
            center_indices = cell(nb_rois, 1);

            % Iterate through each unique value in the atlas image
            for i = 1:nb_rois
                % Find indices where the current value occurs in the atlas image
                if ismember(labels(i), [0, exclude_inds'])
                    continue;
                end
                center_indices{i} = find(atlas_img == labels(i)); % or ismember()

            end
        end


        function drop2nifti(info_nifti, percent_drop, center_nonempty_passing,output_file, other)
            outnifti = zeros([info_nifti.ImageSize(1 : 3), size(percent_drop, 2) + 1]);
            outnifti = drawnifti(outnifti, center_nonempty_passing, ...
                max(percent_drop, [], 2), 1);
            info_nifti.ImageSize = size(outnifti);
            cr = "";
            for i = 1:size(percent_drop, 2)
                msg = compose("brick n %d of %d", i, size(percent_drop, 2));
                fprintf("" + cr + msg);
                cr = repmat('\b', 1, strlength(msg));
                outnifti = drawnifti(outnifti, center_nonempty_passing, ...
                    percent_drop(:, i), i + 1);
            end
            niftiwrite(single(outnifti), output_file, info_nifti, 'Compressed', true);
        end


        function passing_regions2nifti(info_nifti, center_indices_passing, center_indices, output_file, other)
            outnifti = zeros([info_nifti.ImageSize(1 : 3), numel(center_indices_passing)]);
            info_nifti.ImageSize = size(outnifti);
            cr = "";
            for i = 1:numel(center_indices_passing)
                msg = compose("mask n %d of %d", i, numel(center_indices_passing));
                fprintf("" + cr + msg);
                cr = repmat('\b', 1, strlength(msg));

                mask_indices = center_indices{center_indices_passing(i)};
                outnifti = drawnifti(outnifti, {mask_indices}, ...
                    1, i);
            end
            niftiwrite(single(outnifti), output_file, info_nifti, 'Compressed', true);
        end


        function fullcorr2nifti(info_nifti, center_indices, corr_full, corr_res, output_file, other)
            outnifti = zeros([info_nifti.ImageSize(1 : 3), 3]);
            info_nifti.ImageSize = size(outnifti);
            outnifti = drawnifti(outnifti, ...
                center_indices, ...
                squeeze(mean(corr_full(:, :, 1))), ...
                1);

            outnifti = drawnifti(outnifti, ...
                center_indices, ...
                squeeze(mean(corr_res.feat(:, :, 1, 1))), ...
                2);

            outnifti = drawnifti(outnifti, ...
                center_indices, ...
                squeeze(mean(corr_full(:, :, 1)))' - ...
                squeeze(mean(corr_res.feat(:, :, 1, 1)))', ...
                3);
            niftiwrite(single(outnifti), output_file, info_nifti, 'Compressed', true);
        end

    end
end