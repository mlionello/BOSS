
function availableCenters = update_availableCenters(availableCenters, centerIndex, template)
% extract the correct indices in the roi from the template
roiIndices = centerIndex + template.indices;

[centerCoord(1), centerCoord(2), centerCoord(3)] = ind2sub(template.size, centerIndex);
[roiCoord(1, :), roiCoord(2, :), roiCoord(3, :)] = ind2sub(template.size, roiIndices);
nonRoiCoord = sqrt(sum((roiCoord-centerCoord').^2, 1))>=template.padding;

% remove the voxels not within the roi
roiCoord(:, nonRoiCoord) = [];

% get back the indices from the coordinates
rmIndices = sub2ind(template.size, roiCoord(1,:), roiCoord(2,:), roiCoord(3,:));

availableCenters(rmIndices) = 0;
end