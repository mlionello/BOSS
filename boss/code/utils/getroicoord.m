function roiCoords = getroicoord(centerIndex, template)
    roiIndices = centerIndex + template.indices;

    [centerCoord(1), centerCoord(2), centerCoord(3)] = ind2sub(template.size, centerIndex);
    [roiCoords(1, :), roiCoords(2, :), roiCoords(3, :)] = ind2sub(template.size, roiIndices);
    keepIndices = sqrt(sum((roiCoords-centerCoord').^2, 1)) <= template.radius;
    roiCoords = roiCoords(:, keepIndices);
end
