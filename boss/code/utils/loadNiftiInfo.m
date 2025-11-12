function [info_nifti, paths_nifti, sublist] = loadNiftiInfo(output_folder)
    global data_path deriv_path file_pattern;

    [paths_nifti, sublist] = get_nifti_paths(deriv_path, file_pattern);

    info_file = fullfile(data_path, 'nifti_info.mat');
    if exist(info_file, 'file')
        load(info_file, 'info_nifti');
    else
        info_nifti = niftiinfo(paths_nifti{1});
        save(info_file, 'info_nifti');
    end
end
