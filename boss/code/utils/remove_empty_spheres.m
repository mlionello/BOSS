function rmind = remove_empty_spheres(corr_full)

    rmind = [];
    for i = 1: size(corr_full, 2)
        % if at least one subject empty sphere (0s), remove the sphere
        if any(all(corr_full(:, i, :) == 0, 3), 1)
            rmind = [rmind, i];
        end
    end

end
