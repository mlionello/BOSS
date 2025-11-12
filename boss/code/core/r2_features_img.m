function r2 = r2_features_img(nifti_image, mask, feature_names, feature_regressors)
% Perform various analyses including canonical correlation between regions and features
%
% Inputs:
%   - nifti_image: Nifti image data
%   - average_subj_path: Path for average subject data
%   - all_cortex_analysis: Flag indicating whether to perform cortex analysis
%   - sourceRegionMask: Mask for the source region
%   - targetRegionMask: Mask for the target region
%   - feature_names: Set of features to analyze
%   - extra: Flag indicating whether to perform extra analysis
%
% Outputs:
%   - A: Canonical correlation coefficients for target region
%   - B: Canonical correlation coefficients for source region
%   - raw_corr: Raw canonical correlations
%   - rMusSource: Canonical correlation between features and source region
%   - rMusTarget: Canonical correlation between features and target region
%   - resCorr: Residual canonical correlations between regions and features
%
% Example usage:
%   [A, B, raw_corr, rMusSource, rMusTarget, resCorr] = myFunction(nifti_image, sublist, sub_i, average_subj_path, all_cortex_analysis, sourceRegionMask, targetRegionMask, feature_names, extra);

r2 = struct();
msg = compose(" ...extracting 2D Masked Roi");
fprintf(msg);
cr1 = repmat('\b', 1, strlength(msg));
% Perform voxelwise regressions; store 1D-volumes of R2 for each subject
allCortex = extract2DRois(find(mask), nifti_image);
for feat_i = 1:numel(feature_names)
    feat_name = feature_names{feat_i};
    msg = compose(" feat %s: %d/%d", feat_name, feat_i, numel(feature_names));
    fprintf(cr1 + msg);
    cr1 = repmat('\b', 1, strlength(msg));
    r2.(feat_name) = voxelwiseRegr_alt(allCortex, feature_regressors{feat_i});
end
fprintf(cr1);

end