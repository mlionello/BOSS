function [cort_roi_id, cort_roi_label] = extractROIIds(filePath, pattern)
cort_roi_id = [];
cort_roi_label = [];

fid = fopen(filePath, 'r');
if fid == -1
    error('Error opening the file.');
end

while ~feof(fid)
    currentLine = fgetl(fid);
    roi_id = textscan(currentLine, '%f');

    lastString = split(currentLine, ' ');
    containsSubstring = arrayfun(@(x) contains(lastString{end}, x), pattern);
    if any(containsSubstring)
        cort_roi_id = [cort_roi_id; roi_id{1}(1)];
        cort_roi_label = [cort_roi_label; string(lastString{end})];
    end
end
fclose(fid);
end