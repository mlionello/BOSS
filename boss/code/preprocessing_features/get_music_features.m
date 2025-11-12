function params = get_music_features(params, feature_list, number_features_each_model)
    global features_path;

    % Load and concatenate all soundFeat_* files in data_path
    feature_files = dir(fullfile(features_path, 'soundFeat_*.mat'));
    all_features = [];
    feature_names = [];

    for i = 1:numel(feature_files)
        if ~contains(fullfile(features_path, feature_files(i).name), 'pca')
            continue
        end
        loaded_data = load(fullfile(features_path, feature_files(i).name));
        loaded_fields = fieldnames(loaded_data);
        feature2load = find(ismember(loaded_fields, feature_list), 1);
        if isempty(feature2load)
            continue; end
        cell_data = loaded_data.(loaded_fields{feature2load});
        feat_name = loaded_data.featureNames;
        all_features = [all_features, cellfun(@(x) x(:,1:number_features_each_model,:), cell_data, 'UniformOutput', false)];
        feature_names = [feature_names, feat_name];
    end

    %all_features = cellfun(@(x) x(1:91,:,:), all_features, 'UniformOutput', false);

    %all_features = add_missing_timepoints(all_features, onset_list);
    %keep_only_firstNfeatures
    %all_features = cellfun(@(x) x(:, 1:min(number_features_each_model, size(x, 2)), :), all_features, 'UniformOutput', false);

    %[musicRegrPCA_conv, ~, musicRegr_init_conv, ~] = processMusicRegressors(all_features, feature_names, output_folder);
    TR = 1;
    musicRegr_init_conv = convolveMusicRegressors(all_features, spm_hrf(TR));
    musicRegr_init_conv = perform_zsccore(musicRegr_init_conv);
    %musicRegrPCA_conv = perform_zsccore(musicRegrPCA_conv);
    
    % when new features are inserted, add them in constraints feature_set cell
    % array!!!
    % when new features are inserted, add them in constraints feature_set cell
    % array!!!
    
    % assert(all(cellfun(@exist, feature_set)), 'features list and variables do not match: missing corrisponding variable')
    
    %params.musicRegrPCA_conv = musicRegrPCA_conv;
    slices = 5030;

    params.feature_regressors =  cellfun(@(x) x(1:slices, :, :), musicRegr_init_conv, 'UniformOutput', false);
    params.feature_names = feature_names;

end
