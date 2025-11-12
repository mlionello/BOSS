function [roiOut, keepIndices] = extract_roi_from_template(sourceRegion2D, centerIndex, template)
    % extract the correct indices in the roi from the template
    roiIndices = centerIndex + template.indices;

    % for each index in the roi, get the coordinates and check they are
    % one radius away from the given center
    [centerCoord(1), centerCoord(2), centerCoord(3)] = ind2sub(template.size, centerIndex);

    [roiCoord(1, :), roiCoord(2, :), roiCoord(3, :)] = ind2sub(template.size, roiIndices);
    rm_coord = sqrt(sum((roiCoord-centerCoord').^2, 1)) > template.radius;

    % remove the voxels not within the roi
    roiCoord(:, rm_coord) = [];

    % get back the indices from the coordinates
    keepIndices = sub2ind(template.size, roiCoord(1,:), roiCoord(2,:), roiCoord(3,:));
    % keepIndices_new = arrayfun(@(x) find(x==template.mask_indices), keepIndices)

    if numel(sourceRegion2D)==1 && isnan(sourceRegion2D)
        roiOut= nan;
    else
        roiOut = sourceRegion2D(keepIndices, :);
    
        % Remove the indices for those voxels outside the mask region
        roiOut(abs(roiOut(:,1))<1e-6, :) = [];
        roiOut = roiOut';
    end
end
