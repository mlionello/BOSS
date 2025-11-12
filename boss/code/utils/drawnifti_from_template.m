function outnifti = drawnifti_from_template(outnifti, indices, template, values, brick)
    values_counter = zeros(size(outnifti));
    for i = 1:length(values)
        if isnan(indices(i))
            continue
        end
        roiCoords = getroicoord(indices(i), template);
        for j = 1: size(roiCoords, 2)
            x = roiCoords(1, j); y = roiCoords(2, j); z = roiCoords(3, j);

            outnifti(x, y, z, brick) = outnifti(x, y, z, brick) + values(i);
            values_counter(x, y, z, brick) = values_counter(x, y, z, brick) + 1;
        end
    end
    nonzeroes_counter = values_counter > 0;
    outnifti(nonzeroes_counter) = outnifti(nonzeroes_counter) ./ values_counter(nonzeroes_counter);
end