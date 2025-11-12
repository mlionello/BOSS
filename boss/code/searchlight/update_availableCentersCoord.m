function updatedCenters = update_availableCentersCoord(availableCenters, currentCenter, padding)
% Get the size of the availableCenters matrix
[rows, cols, slices] = size(availableCenters);

% Create a grid of coordinates
[X, Y, Z] = ndgrid(1:rows, 1:cols, 1:slices);

% Compute the distance from each point to the center
distances = sqrt((X - currentCenter(1)).^2 + (Y - currentCenter(2)).^2 + (Z - currentCenter(3)).^2);

% Define a spherical mask for the ROI
roi_mask = distances < padding;

% Update the availability status of neighboring centers
updatedCenters = availableCenters & ~roi_mask;
end

