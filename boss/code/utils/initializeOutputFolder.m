function output_folder = initializeOutputFolder(opts)
    global res_path;
    
    if ~isempty(opts.checkpoint)
        output_folder = fullfile(res_path, opts.checkpoint);
        output_folder = dir([output_folder '*/']).folder;
        fprintf("Using checkpoint: %s\n", output_folder);
    else
        output_folder = fullfile(res_path, ...
            string(datetime('now', 'Format', 'yyyyMMdd_HHmmss')) + opts.suffix_folder);
        mkdir(output_folder);
    end
end
