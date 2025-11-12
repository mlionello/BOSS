function [pval, corr] = corr_and_maxstat_on_perm(x, y, permSchema, feat)
    % source_canvars and target_canvars: are timpeoints by voxels
    % perm schema: is the schema to permute timepoints
    % feat: is eventually the model from which extract the residual from the
    % regions
    corr = zeros(size(permSchema, 2), size(x, 2));
    for k = 1: size(permSchema, 2)
        % if feat is given, correlation between residual x-feat*(feat\x) 
        % from both regions, otherwise, simple correlation between regions
        if nargin > 3 && (~isempty(feat) || ~any(isnan(feat)))
            corr(k, :) = fast_corr( ...
                get_residual(x, feat(permSchema(:, k), :)), ...
                get_residual(y, feat(permSchema(:, k), :)) ...
                );
        else
            corr(k, :) = fast_corr(x, y(permSchema(:, k), :));
        end
    end
    
    %calculate p value via maxstatistic across voxels
    pval = nan;
    if size(corr, 1) > 1
        pval = get_pval_from_maxstat(corr);
    end
    corr = corr(1, :);
end
