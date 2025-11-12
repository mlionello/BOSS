function roi2D = extract2DRois(indices, nifti_image)
    roi2D = nan(size(nifti_image, 4), length(indices));
    nb_vx = prod(size(nifti_image, 1: 3));
    for indx = 1: length(indices)
        roi2D(:, indx) = nifti_image( ...
            indices(indx) : ...
            nb_vx : ...
            end);
    end
end
