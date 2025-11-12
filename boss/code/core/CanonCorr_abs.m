classdef (Abstract) CanonCorr_abs
    properties
        output_folder
        coeff
        sourceRegionMask
        targetRegionMask
        feature_regressors
        feature_names
        winkler
        nb_perm
        extra
    end
    
    methods
        obj = fit(obj, nifti_image_avg)
        results = predict(obj, input_image)
    end
end