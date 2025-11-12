function r2 = voxelwiseRegr(roi, infeat, permSchema)
    r2 = zeros(size(roi, 2), size(permSchema, 2));
    cr = "";
    for i = 1:size(roi, 2)
        if isempty(getCurrentJob)
            msg = compose(" permuting vx %06d/%d", i, size(roi, 2));
            fprintf(cr + msg);
            cr = repmat('\b', 1, strlength(msg));
        end
        for k = 1:size(permSchema, 2)
            perm_infeat = infeat(permSchema(:, k), :);
            b = perm_infeat\roi(:, i);
            y = perm_infeat*b;
            ssres = sum((y - roi(:, i)).^2);
            sstot = sum((roi(:, i) - mean(roi(:, i))).^2);
            r2(i, k) = 1 - ssres/sstot;
        end
    end
    if isempty(getCurrentJob); fprintf(cr); end
end
