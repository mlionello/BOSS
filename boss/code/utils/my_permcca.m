function [pval, r, A, B, U, V] = my_permcca(targetRoiAVG, sourceRoiAVG, permschema)
    r = zeros(size(permschema, 2), min(size(targetRoiAVG,2), size(sourceRoiAVG, 2)) );
    [A, B, r(1,:), U, V, ~] = canoncorr(targetRoiAVG, sourceRoiAVG(permschema(:, 1),:));

    parfor k = 2:size(permschema, 2)
        [~, ~, r(k,:), ~, ~, ~] = canoncorr(targetRoiAVG, sourceRoiAVG(permschema(:, k), :));
    end

    r_null = max(r(2: end, :), [], 2);
    pval = tiedrank(-cat(1, r(1, :), repmat(r_null, 1, size(r, 2))))/size(r, 1);
    pval = pval(1, :);
end