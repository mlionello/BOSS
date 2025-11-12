function [roi_id, roi_label] = extractROIFromLabels(atlas_path, labels_file, pattern)
    labels_file_path = fullfile(atlas_path, labels_file);
    [roi_id, roi_label] = extractROIIds(labels_file_path, pattern);
end
