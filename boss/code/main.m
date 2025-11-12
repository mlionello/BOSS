function main(roi1_pattern, roi2_pattern, opts)

arguments
    roi1_pattern;
    roi2_pattern;
    opts.channel = ["left", "right"];
    opts.nb_perm_step1 = 0;
    opts.smooth_regions_padding = 0;
    opts.smooth_cortex_padding = 0;
    opts.numPermutations_individual = 0;
    opts.nb_perm = 0;
    opts.nb_signflip = 10000;
    opts.nb_permAVG = 0;
    opts.extra = 0;
    opts.winkler = 0;
    opts.checkpoint = [];
    opts.suffix_folder = [];
    opts.suffix_results = [];
    opts.cpus = 4;
    opts.target_model = "CanonCorr";
    opts.radius = [1 2 3 4];
    opts.padding = [1 1 2 3];
    opts.pVxIn = 1;
    opts.feature_list = [];
    opts.update_params = 0;
end

    addpath('../');
    constants();
    % Set up default options and validate inputs
    setupEnvironment(opts);
    opts.suffix_folder = getSuffix(opts.suffix_folder);
    opts.suffix_results = getSuffix(opts.suffix_results);

    % Initialize constants and check for existing checkpoint
    global meansubject_path commonmask_path;

    output_folder = initializeOutputFolder(opts);

    % Load or create masks and parameters
    [masks, params, downsampleRoi2, opts] = initializeMasksAndParams(roi1_pattern, roi2_pattern, opts, output_folder);

    % Load NIFTI info
    [info_nifti, paths_nifti, sublist] = loadNiftiInfo(output_folder);

    writeMasksNifti(output_folder, masks, info_nifti);

    % Execute analysis steps
    average_subj_anslysis(info_nifti, paths_nifti, ...
            masks, params.feature_names, sublist, ...
            params.feature_regressors, output_folder, meansubject_path, commonmask_path);

    init_model = init_params4model(output_folder, params, masks, opts);
    model = fit_model(meansubject_path, output_folder, opts.target_model, init_model, opts);


    if exist(fullfile(output_folder, "passing_centers.mat"))
        load(fullfile(output_folder, "passing_centers.mat"));
        load(fullfile(output_folder, "passing_indices.mat"));
        load(fullfile(output_folder, "nb1.mat"));
    else
        %%% Compute Correlation Subject-wise all spheres
        [~, corr_full, ~, corr_res, ...
                    ~] = predict_individuals(output_folder, opts.suffix_results, paths_nifti, model, params, opts);
    
        [corr_full, corr_res, nb1] = reshape_corr_res(corr_full, corr_res);
    
        %%% Compute Sign flipping
        [~, passing_centers, passing_indices, ~] = iterate4signFlipping(output_folder, opts.suffix_results, corr_full, corr_res, nb1, opts.nb_signflip, model.center_indices);
        save(fullfile(output_folder, "passing_centers.mat"), 'passing_centers');
        save(fullfile(output_folder, "passing_indices.mat"), 'passing_indices');
        save(fullfile(output_folder, "nb1.mat"), 'nb1');
    end

    %%% Compute Drop Subject-wise on passing spheres
    [corr_res_perm, corr_full_passing] = compute_drop(paths_nifti, ...
        opts.suffix_results, model, output_folder, ...
        passing_centers.fivecent, params);
    [corr_full_passing_reshaped, corr_res_perm, ~] = reshape_corr_res(corr_full_passing, corr_res_perm);

    for rad = 1:size(corr_res_perm, 1)
        for ch = 1:size(corr_res_perm, 2)
            full_ref = corr_full_passing_reshaped{rad, ch};

            % these are the union of indices across features as seen during drop calc
            union_passing_indices = [passing_indices.fivecent{rad, ch, :}];
            union_passing_indices = any(union_passing_indices, 2);

            union_passing_centers = vertcat(passing_centers.fivecent{rad, ch, :});
            union_passing_centers = unique(union_passing_centers);

            for feat_i = 1 : size(corr_res_perm, 3)
                if isempty(corr_res_perm{rad, ch, feat_i})
                    continue
                end
                % corr_res was calculated on the overall target region
                corr_res_reshaped_tmp = corr_res{rad,ch,feat_i};
                corr_res_reshaped_tmp = corr_res_reshaped_tmp(:,passing_indices.fivecent{rad,ch,feat_i},:);

                % corr_res_perm was calculated only on the union of the passing indices across the features
                % check which common index across features is presente in the current passing indices set
                indices_from_union_to_feat = ismember(find(union_passing_indices), find(passing_indices.fivecent{rad,ch,feat_i}));
                row = corr_res_perm{rad, ch, feat_i};
                row = row(:,indices_from_union_to_feat,:);

%                 percent_drop = (1 - corr_res_reshaped_tmp ./ row) * 100;
%                 mean_perc_drop = squeeze(mean(percent_drop, 1));
                if contains(params.target_model, "mean2mean")
                    mean_perc_drop = squeeze(1 - mean(corr_res_reshaped_tmp, 1) ./ mean(row, [1,3])) * 100;
                else
                    mean_perc_drop = squeeze(1 - mean(corr_res_reshaped_tmp, 1) ./ mean(row, 1)) * 100;
                end
                % mean_perc_drop = squeeze(1 - mean(corr_res_reshaped_tmp, 1) ./ mean(corr_full_reshaped_tmp, 1)) * 100;
                template = model.targetRoiTemplate{rad};

                % Get the corresponding parameters
                current_radius = params.radius(rad); % Assuming params.radius is an array
                current_padding = params.padding(rad); % Assuming params.padding is an array
                feature_name = params.feature_names{feat_i}; % Assuming feature_names is a cell array

                % Create an informative file name
                file_name = fullfile(output_folder, sprintf("Result_rad%d_pad%d_feat_%s_ch_%d%s.nii", ...
                    current_radius, current_padding, feature_name, ch, opts.suffix_results));

                %file_name = fullfile(output_folder, "OMG" + opts.suffix_results); % add rad, ch, feat
                model.drop2nifti(info_nifti, mean_perc_drop, passing_centers.fivecent{rad, ch, feat_i}, file_name, template)
                % model.fullcorr2nifti(info_nifti, passing_centers.fivecent{rad, ch}, squeeze(mean(full_ref)), squeeze(mean(data_extractor)), fullfile(output_folder, "full_corr" + opts.suffix_results), template)
            end

        end
    end
    %[data_extractor(:,1,1)'; full_ref(:,1,1)']'

%     model.passing_regions2nifti(info_nifti, center_indices_passing, model.center_indices, fullfile(output_folder, "OMG_sorted_mask" + opts.suffix_results), template)
%     model.fullcorr2nifti(info_nifti, model.center_indices, corr_full, corr_res, fullfile(output_folder, "full_corr" + opts.suffix_results), opts.suffix_results)
% 
%     %%% Results to NIFTI
%     fprintf("drawing the results2nifti\n")
%     draw_ttest_results_on_subj_AVGcanocorr(model.AVG_results, ...
%         masks.roi2_mask, masks.roi1_mask, ...
%         model.targetRoiTemplate, info_nifti, output_folder, ...
%         target_canvars_mean_activity, source_canvars_mean_activity, pval)
% 
%     % Print completion message
%     fprintf("Happily ever after! You are done!\n");
end

% --- Helper function to set up the environment ---
function setupEnvironment(opts)
    global pkg_folder;

    addpath(fullfile(pkg_folder, 'code'));
    addpath(fullfile(pkg_folder, 'code', 'utils'), ...
        fullfile(pkg_folder, 'code', 'core/'), ...
        fullfile(pkg_folder, 'code', 'extra/'), ...
        fullfile(pkg_folder, 'code', 'searchlight/'));
    addpath('~/Matlab/spm12/');
    addpath('~/Matlab/MultiSolverFolder/');
    addpath('~/Matlab/MultipleQR/');
    maxNumCompThreads(opts.cpus);
    %setupParpool(opts.cpus);
end

% --- Helper function for setting up parallel pool ---
function setupParpool(cpus)
    try
        p = gcp('nocreate');
        if isempty(p) || p.NumWorkers ~= cpus
            delete(gcp('nocreate'));
            parpool(cpus);
        end
    catch
        warning('Could not set up parallel pool.');
    end
end

