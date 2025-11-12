classdef CanonCorr < CanonCorr_Base
    properties
        AVG_resPval;
        AVG_resCorr;
    end
    methods
        function obj = CanonCorr(opts)
            arguments
                opts (1,1) struct
            end
            obj = obj@CanonCorr_Base(opts);
        end

        function obj = load(obj, AVG_targetsource)
            obj.coeff = AVG_targetsource;
        end

        function obj = fit(obj, nifti_image_avg)

            coeff = analyze_average_subject(nifti_image_avg, ...
                obj.sourceRegionMask, obj.targetRegionMask, obj.feature_regressors, ...
                obj.winkler, obj.nb_perm, obj.extra);

            obj.coeff = coeff;

            permSchema = initPermutationTests(obj.nb_perm, size(obj.coeff.U, 1));
            [obj.AVG_resPval, obj.AVG_resCorr] = corr_and_maxstat_on_perm(obj.coeff.U, ...
                obj.coeff.V, permSchema, obj.feature_regressors);

        end

        function [pval_full, corr_full, ...
                pval_res, corr_res, ...
                source_canvars_mean_activity, ...
                target_canvars_mean_activity] = predict(obj, nifti_image)
        
            % extract from the nifti the timepoints by voxels images for the given a mask
            sourceRoi = extract2DRois(find(obj.sourceRegionMask), nifti_image);
            targetRoi = extract2DRois(find(obj.targetRegionMask), nifti_image);

            [pval_full, corr_full, pval_res, corr_res, ...
                source_canvars_mean_activity, ...
                target_canvars_mean_activity] = gli2subjects( ...
                sourceRoi, targetRoi, obj.coeff, obj.nb_perm, obj.feature_regressors);            
        end
        
    end

end