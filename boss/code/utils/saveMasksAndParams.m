function saveMasksAndParams(output_folder, masks, params)
    masks_file = fullfile(output_folder, 'masks.mat');
    params_file = fullfile(output_folder, 'params.mat');
    save(masks_file, 'masks');
    save(params_file, 'params');
    fprintf("Saved masks and params to %s and %s\n", masks_file, params_file);
end
