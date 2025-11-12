function values_out = get_roi(region, center, radius)
% Extracts a spherical region of interest (ROI) from a 3D binary region
% centered at the specified coordinates with the given radius.

% Ensure the input coordinates are valid
if ~isnumeric(center) || numel(center) ~= 3
    error('Invalid center coordinates. Must be a 3-element numeric array.');
end

% Ensure the input region is a 4D array
if ndims(region) ~= 4
    error('Invalid region. Must be a 3D binary array.');
end

% Get the size of the region
[cols, rows, slices, t_points] = size(region);

% Ensure the center coordinates are within the valid range
if any(center < 1) || any(center > [rows, cols, slices])
    error('Center coordinates are outside the valid range.');
end

% Create a grid of coordinates
[X, Y, Z] = ndgrid(1:cols, 1:rows, 1:slices);

% Compute the distance from each point to the center
distances = sqrt((X - center(1)).^2 + (Y - center(2)).^2 + (Z - center(3)).^2);

% Create a binary mask for the spherical ROI
roi_mask = distances <= radius;

% Extract the values within the ROI
nonzeros_indices = find(roi_mask);
numel_mask = numel(roi_mask);
values_out = zeros(size(region, 4), numel(nonzeros_indices));
for n = 1:numel(nonzeros_indices)
    values_out(:, n) = squeeze(region(nonzeros_indices(n):numel_mask:end));
end
end