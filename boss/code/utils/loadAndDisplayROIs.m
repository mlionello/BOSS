function [roi1_id, roi1_label, roi2_id, roi2_label] = loadAndDisplayROIs(atlas_path, roi1_pattern, cort_labels, roi2_pattern, subcort_labels, HCP_labels)
    [roi1_id, roi1_label] = extractROIIds( roi1_pattern, cort_labels, subcort_labels, HCP_labels);
    displayROIIds(roi1_pattern, roi1_id)
    if matches(roi2_pattern, "allbrain")
        roi2_id = {nan};
        roi2_label = {"allbrain"};
        fprintf("vs all brain\n")
    else
        [roi2_id_cell, roi2_label_cell] = extractROIIds(roi2_pattern, cort_labels, subcort_labels, HCP_labels);
        
        roi2_id = [roi2_id_cell{:}];
        roi2_label = [roi2_label_cell{:}];

        displayROIIds(roi2_pattern, roi2_id)
    end

end

function [roi_id_out, roi_label_out] = extractROIIds(pattern, CORT_LABELS, SUBCORT_LABELS, HCP_LABELS)
    
    roi_label_out = {[], []};  % Cell array to hold L and R results separately
    roi_id_out = {[], []};
    
    for i = 1:length(pattern)
        split_pattern = strsplit(pattern{i}, ':');
        extracted_pattern = split_pattern{2};
        
        % Check if the extracted pattern contains 'L_' or 'R_'
        if ~contains(extracted_pattern, 'L_') && ~contains(extracted_pattern, 'R_')
            % Create patterns for both 'L_' and 'R_'
            extracted_patterns = {['L_' extracted_pattern], ['R_' extracted_pattern]};
        else
            % If it already contains 'L_' or 'R_', use the extracted pattern as is
            extracted_patterns = {extracted_pattern};
        end
        
        for j = 1:length(extracted_patterns)
            switch split_pattern{1}
                case 'cort'
                    [roi_id, roi_label] = extractROIIds_interate_atlas(CORT_LABELS, extracted_patterns{j});
                case 'sub'
                    [roi_id, roi_label] = extractROIIds_interate_atlas(SUBCORT_LABELS, extracted_patterns{j});
                case 'hcp'
                    [roi_id, roi_label] = extractROIIds_interate_atlas(HCP_LABELS, extracted_patterns{j});
            end
            
            % Determine if it's L or R case
            if contains(extracted_patterns{j}, 'L_')
                idx = 1;  % Index for L
            else
                idx = 2;  % Index for R
            end
            
            roi_label_out{idx} = [roi_label_out{idx}, strcat(split_pattern{1}, ":", roi_label)];
            roi_id_out{idx} = [roi_id_out{idx}, roi_id];
        end
    end
end


function [roi_id, roi_label] = extractROIIds_interate_atlas(file_path, pattern)
    roi_id = [];
    roi_label = {};
    fid = fopen(file_path, 'r');
    tline = fgetl(fid);
    while ischar(tline)
        current_id = sscanf(tline, '%f');
        last_string = strsplit(tline);
        last_string = last_string{end};
        if any(contains(lower(last_string), lower(pattern)))
            roi_id = [roi_id, int32(current_id(1))];
            roi_label{end+1} = last_string;
        end
        tline = fgetl(fid);
    end
    fclose(fid);
end
