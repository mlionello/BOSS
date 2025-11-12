classdef Searchlight_CanonCorr_pca2sphere < CanonCorr_Base
    properties
        targetPropInVx
        targetRoiTemplate
        analyserapply
        analyserfit
        mdata
        targetRadius
        targetPadding
        saveactivity
        featOnTarget
        img_size

        center_indices
        pca_source
    end
    methods
        function obj = Searchlight_CanonCorr_pca2sphere(opts)
            arguments
                opts (1, 1) struct
            end
            obj = obj@CanonCorr_Base(opts);
            fields = {'targetRadius', 'targetPadding', 'targetPropInVx', ...
                'analyserapply', 'analyserfit', 'saveactivity'};
            init = struct('targetRadius', [1 2 3 4], 'targetPadding', [1 1 2 3], ...
                'targetPropInVx', 1, 'analyserfit', @perform_cca_on_searchlight.fit, ...
                'analyserapply', @perform_cca_on_searchlight.apply, ...
                'saveactivity', 0, 'featOnTarget', []);
            obj = parse_args(obj, fields, init, opts);

            %assert(isscalar(unique(cellfun(@(x) size(x), obj.targetRegionMask))));
            %obj.img_size = size(obj.targetRegionMask{1});
            obj.img_size = size(obj.targetRegionMask);

            for k = 1:numel(obj.targetRadius)
                obj.targetRoiTemplate{k} = create_template(obj.img_size, ...
                    obj.targetRadius(k), obj.targetPadding(k));
            end

            obj.mdata = struct('targetRegionMask', obj.targetRegionMask, ...
                'winkler', obj.winkler, ...
                'nb_perm', obj.nb_perm, ...
                'extra', obj.extra);
            obj.mdata.feature_regressors = obj.feature_regressors;
            obj.mdata.sourceRegionMask = obj.sourceRegionMask;
            obj.coeff = cell(numel(obj.targetRadius), 2);
            obj.center_indices = cell(numel(obj.targetRadius), 2);
            obj.featOnTarget = cell(numel(obj.targetRadius), 2);
            obj.pca_source = cell(numel(obj.targetRadius), 2);
        end

        function obj = update_params(obj, params)
            if params.changed_features
                obj.feature_regressors = params.feature_regressors;
                obj.mdata.feature_regressors = params.feature_regressors;
                obj.feature_names = params.feature_names;
            end
            if ~isequal(params.targetRadius, obj.targetRadius) || ~isequal(params.targetPadding, obj.targetPadding)
                for i = 1:length(params.targetRadius)
                    rad_i = params.targetRadius(i);
                    padd_i = params.targetPadding(i);

                    indx = find(obj.targetRadius == rad_i);
                    if ~any(padd_i == obj.targetPadding(indx))
                        error("radius/padding pair subsampling not found in training parameters");
                    end
                end
            end
        end

        function obj = fit(obj, nifti_image_avg, opts)
            arguments
                obj;
                nifti_image_avg;
                opts.channel = ["left", "right"];
                opts.radius = obj.targetRadius;
                opts.out_dir = "./";
            end

            for channel = opts.channel
                switch channel
                    case "right"
                        region_index = 1; % source on left emisphere
                        ch_indx = 2; % right channel
                    case "left"
                        region_index = 2; % source on right emisphere
                        ch_indx = 1; % left channel
                end
    
                targetRoi = nifti_image_avg .* obj.targetRegionMask; % ALL BRAIN
                targetRoi2D = reshape(targetRoi, [], size(targetRoi, 4));

                indices_source = find(obj.sourceRegionMask{region_index});
                sourceRegion2D = extract2DRois(indices_source, nifti_image_avg);

                target_mask = obj.targetRegionMask;
                target_mask(indices_source) = 0;

                for radius = opts.radius
                    template_index = find(obj.targetRadius == radius, 1);

                    [pca_on_source{template_index, ch_indx}, source_pca_score] = pca(sourceRegion2D, "NumComponents", min(obj.targetRoiTemplate{template_index}.maxnbvx, size(sourceRegion2D,2)));

                    % Construct file paths using file_out_id and index
                    coeff_file = fullfile(opts.out_dir, sprintf('coeff_ch_%s_rad_%d.mat', channel, radius));
                    center_indices_file = fullfile(opts.out_dir, sprintf('center_indices_ch_%s_rad_%d.mat', channel, radius));
                    pca_source_file = fullfile(opts.out_dir, sprintf('pca_on_source_ch_%s_rad_%d.mat', channel, radius));

                    if ~exist(coeff_file, 'file')
                        mdata_only_one_channel = obj.mdata;
                        mdata_only_one_channel.feature_regressors = cellfun(@(x) x(:, :, ch_indx), mdata_only_one_channel.feature_regressors, 'UniformOutput', false);

                        [coeff, center_indices] = searchlight(source_pca_score, obj.targetPropInVx, ...
                            target_mask, targetRoi2D, ...
                            obj.targetRoiTemplate{template_index}, obj.analyserfit, mdata_only_one_channel); % for template in templates (radius range)

                        % Save results to file
                        save(coeff_file, 'coeff');
                        save(center_indices_file, 'center_indices');
                        save(pca_source_file, 'pca_on_source')
                    end

                    % Store file paths
                    if isempty(obj.coeff{template_index, ch_indx})
                        obj.coeff{template_index, ch_indx} = coeff_file; end
                    if isempty(obj.center_indices{template_index, ch_indx})
                        obj.center_indices{template_index, ch_indx} = center_indices_file; end
                    if isempty(obj.pca_source{template_index, ch_indx})
                        obj.pca_source{template_index, ch_indx} = pca_source_file; end
                end
            end
        end

        function [pval_full, corr_full, ...
                pval_res, corr_res, corr_res_perm] = predict(obj, nifti_image, vx_shuffling, feat_shuffling, opts)
            arguments
                obj;
                nifti_image;
                vx_shuffling;
                feat_shuffling;
                opts.channel = ["left", "right"];
                opts.radius = obj.targetRadius;
                opts.out_dir = "./";
                opts.sub_indices = [];
            end
            warning('off');
            cr0 = ""; cr = "";

            for channel = opts.channel
                switch channel
                    case "right"
                        region_index = 1; % source on left emisphere
                        ch_indx = 2; % right channel
                    case "left"
                        region_index = 2; % source on right emisphere
                        ch_indx = 1; % left channel
                end

                sourceRoi = single(extract2DRois(find(obj.sourceRegionMask{region_index}), nifti_image));

                for radius = opts.radius

                    template_index = find(obj.targetRadius == radius, 1);

                    coeff_file = obj.coeff{template_index, ch_indx};
                    center_indices_file = obj.center_indices{template_index, ch_indx};
                    pca_source_file = obj.pca_source{template_index, ch_indx};

                    % Load the data from saved files
                    if ~isempty(coeff_file) && isfile(coeff_file)
                        load(coeff_file, 'coeff');
                    else
                        error('Coefficient file not found for channel %s and radius %d.', channel, radius);
                    end

                    if ~isempty(center_indices_file) && isfile(center_indices_file)
                        load( center_indices_file, 'center_indices');
                    else
                        error('Center indices file not found for channel %s and radius %d.', channel, radius);
                    end

                    if ~isempty(pca_source_file) && isfile(pca_source_file)
                        load( pca_source_file, 'pca_on_source');
                    else
                        error('Center indices file not found for channel %s and radius %d.', channel, radius);
                    end                    

                    if ~isempty(opts.sub_indices)
                        union_passing_centers = vertcat(opts.sub_indices{template_index, ch_indx, :}); % Concatenation
                        
                        % Compute the union of all elements in the array
                        union_passing_centers = unique(union_passing_centers);
                        for feat_i=1:size(opts.sub_indices, 3)
                            assert(all(ismember(opts.sub_indices{template_index, ch_indx, feat_i}, union_passing_centers)))
                            assert(numel(ismember( ...
                                opts.sub_indices{template_index, ch_indx, feat_i}, ...
                                union_passing_centers)) ...
                                == numel(opts.sub_indices{template_index, ch_indx, feat_i}) ...
                                )
                        end

                        logical_indices = ismember(center_indices, union_passing_centers);
                        coeff = coeff(logical_indices);
                        center_indices = center_indices(logical_indices);

                        assert(all(center_indices == union_passing_centers))
                    end

                    feature_regressors = cellfun(@(x) single(x(:,:,ch_indx)), obj.feature_regressors, 'UniformOutput',false);
                    saveactivity_run = obj.saveactivity;
                    targetRoiTemplate_run = obj.targetRoiTemplate{template_index};

                    source_pca_score = sourceRoi*pca_on_source{template_index, ch_indx};

                    % parfor_progress(numel(indices));
                    init_time = tic;
                    last_toc = toc(init_time);
                    msg0 = compose("channel: %s radius: %d ", channel, radius);
                    fprintf("" + cr0 + cr + msg0);
                    cr0 = repmat('\b', 1, strlength(msg0));

                    cr = "";
                    for i = 1:numel(center_indices)
                        current_indx = center_indices(i);
                        %if isempty(getCurrentTask()) && toc(init_time) > last_toc + 1
                        if toc(init_time) > last_toc + 1
                            msg = compose("%d out of %d; elapsed %s; eta: %s", i, numel(center_indices), ...
                                duration(0,0,toc(init_time)), ...
                                duration(0,0,toc(init_time)/i*(numel(center_indices)-i)));
                            fprintf(cr + msg);
                            cr = repmat('\b', 1, strlength(msg));
                            last_toc = toc(init_time);
                        %elseif toc(init_time) > last_toc + 1
                        %    parfor_progress;
                        end

                        [~, keepIndices] = extract_roi_from_template(nan, current_indx, targetRoiTemplate_run);    

                        % extract from the nifti the timepoints by voxels images for the given a mask
                        targetRoi = single(extract2DRois(keepIndices, nifti_image));

                        [pval_full{i, template_index, ch_indx}, ...
                            corr_full{i, template_index, ch_indx}, ...
                            pval_res{i, template_index, ch_indx}, ...
                            corr_res{i, template_index, ch_indx}, ...
                            ~, ...
                            corr_res_perm{i, template_index, ch_indx}] = gli2subjects( ...
                                source_pca_score, targetRoi, feature_regressors, 'coeff', coeff{i}, ...
                                'vx_shuffling', vx_shuffling, 'save_activity', saveactivity_run, ...
                                'nb_feat_shuffling', feat_shuffling, 'analyser', obj.analyserapply); % features must be passed as cell array
                    end

                end
            end
            fprintf("" + cr + cr0);
            % parfor_progress(0);

        end

    end

    methods (Static)

        function drop2nifti(info_nifti, percent_drop, center_nonempty_passing,output_file, tempalte)
            if all(size(percent_drop,[1, 2]) > [1,1])
                percent_drop = percent_drop(:,1);
            elseif numel(center_nonempty_passing)==1
                percent_drop = percent_drop(1); % in the case there is only one index
            else
                fprintf("error between indices and percdrop size: %d vs %s\n", numel(center_nonempty_passing), string(size(percent_drop)).join(" x "))
            end
            outnifti = zeros([info_nifti.ImageSize(1 : 3), size(percent_drop, 2) + 1]);
            outnifti = drawnifti_from_template(outnifti, center_nonempty_passing, ...
                tempalte, max(percent_drop, [], 2), 1);
            info_nifti.ImageSize = size(outnifti);
            cr = "";
            for i = 1:min(10, size(percent_drop, 2))
                msg = compose("brick n %d of %d", i, size(percent_drop, 2));
                fprintf("" + cr + msg);
                cr = repmat('\b', 1, strlength(msg));
                outnifti = drawnifti_from_template(outnifti, center_nonempty_passing, ...
                    tempalte, percent_drop(:, i), i + 1);
            end
            niftiwrite(single(outnifti), output_file, info_nifti, 'Compressed', true);
        end


        function passing_regions2nifti(info_nifti, center_indices_passing, center_indices, output_file, tempalte)
            outnifti = zeros([info_nifti.ImageSize(1 : 3), numel(center_indices_passing)]);
            info_nifti.ImageSize = size(outnifti);
            for i = 1:numel(center_indices_passing)

                mask_indices = center_indices(center_indices_passing(i));
                outnifti = drawnifti_from_template(outnifti, mask_indices, ...
                    tempalte, 1, i);
            end
            niftiwrite(single(outnifti), output_file, info_nifti, 'Compressed', true);
        end


        function fullcorr2nifti(info_nifti, center_indices, corr_full, corr_res, output_file, template)
            outnifti = zeros([info_nifti.ImageSize(1 : 3), 3]);
            info_nifti.ImageSize = size(outnifti);
            outnifti = drawnifti_from_template(outnifti, ...
                center_indices, ...
                template, ...
                squeeze(mean(corr_full(:, :, 1))), ...
                1);

            outnifti = drawnifti_from_template(outnifti, ...
                center_indices, ...
                template, ...
                squeeze(mean(corr_res(:, :, 1, 1))), ...
                2);

            outnifti = drawnifti_from_template(outnifti, ...
                center_indices, ...
                template, ...
                squeeze(mean(corr_full(:,:,1)))' - ...
                squeeze(mean(corr_res(:, :, 1, 1)))', ...
                3);
            niftiwrite(single(outnifti), output_file, info_nifti, 'Compressed', true);
        end
    end

end