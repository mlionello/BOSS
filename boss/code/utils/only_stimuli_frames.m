function only_stimuli_frames(path_files, out_root)
    % Ensure output directory exists
    if ~exist(out_root, 'dir')
        mkdir(out_root);
    end

    % Loop through each NIfTI file
    for i = 1:length(path_files)
        % Load the NIfTI file
        nifti_file = path_files{i};
        nii_data = niftiread(nifti_file);
        nii_info = niftiinfo(nifti_file);
        
        % Determine the corresponding CSV file path
        [file_folder, file_name, ~] = fileparts(nifti_file);
        csv_file = fullfile(file_folder, [file_name '.csv']);
        
        % Load the CSV file
        stim_data = readtable(csv_file, 'ReadVariableNames', false); % Read without variable names
        
        % Extract onset and duration columns by their index
        onsets = stim_data{:, 1};  % First column for onset
        durations = stim_data{:, 2};  % Second column for duration
        
        % Calculate the onset and offset frames
        TR = nii_info.PixelDimensions(4); % Assuming the 4th dimension is time (TR)
        onset_frames = round(onsets / TR) + 1; % Add 1 to convert to 1-based index
        offset_frames = round((onsets + durations) / TR) + 1;
        
        % Determine the frames to keep
        frames_to_keep = [];
        for j = 1:length(onset_frames)
            frames_to_keep = [frames_to_keep, onset_frames(j):offset_frames(j)];
        end
        
        % Remove duplicates and sort frames_to_keep
        frames_to_keep = unique(frames_to_keep);
        % remove first 26 frames (length of kernel filter in afni
        % 3ddeconvolve)
        frames_to_keep = frames_to_keep(frames_to_keep>26);
        
        % Keep only the required bricks
        nii_data_modified = nii_data(:, :, :, frames_to_keep);
        
        % Update the nii_info dimensions
        nii_info.ImageSize(4) = length(frames_to_keep);
        
        % Save the modified NIfTI file
        output_file = fullfile(out_root, [file_name '_modified.nii']);
        niftiwrite(single(nii_data_modified), output_file, nii_info);
        
        fprintf('Processed and saved: %s\n', output_file);
    end
end
