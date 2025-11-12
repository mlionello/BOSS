function displayROIIds(pattern, roi_id)
    fprintf('ROI IDs corresponding to lines with %s:', mat2str(pattern));
    fprintf("%d\n", cell2mat(roi_id'));
end
