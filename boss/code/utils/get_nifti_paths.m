function [full_paths, main_subdir_names] = get_nifti_paths(main_directory, file_pattern)
% Function: get_nifti_paths
%
% Description:
%   Retrieves full paths of files matching a specified pattern within a main directory,
%   along with the names of the main subdirectories containing these files.
%
% Inputs:
%   main_directory - The main directory containing the files.
%   file_pattern   - The pattern of the files to search for.
%
% Outputs:
%   full_paths         - A cell array containing the full paths of the matched files.
%   main_subdir_names  - A cell array containing the names of the main subdirectories
%                        containing the matched files.
%
% Example:
%   main_directory = 'C:\Data\derivatives\';
%   file_pattern = '*.nii';
%   [full_paths, main_subdir_names] = get_nifti_paths(main_directory, file_pattern);

file_list = dir(fullfile(main_directory, '**', file_pattern));

% Extract full paths and corresponding main subdirectory names
full_paths = cell(numel(file_list), 1);
main_subdir_names = cell(numel(file_list), 1);

for i = 1:numel(file_list)
    full_paths{i} = fullfile(file_list(i).folder, file_list(i).name);

    % Extract main subdirectory name
    subdirectory_parts = strsplit(file_list(i).folder, filesep);
    main_subdir_names{i} = subdirectory_parts{numel(strsplit(main_directory, filesep)) + 1};
end
end
