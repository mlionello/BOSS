classdef perform_cca_on_searchlight

    methods (Static)

        function results = fit(targetRoi, sourceRoi, mdata)
            % perform canonical correlation between the two ROIs and
            warning('off')
            [A, B, r, ~, ~, ~] = canoncorr(targetRoi, sourceRoi);
            results.A = A;
            results.B = B;
            results.r = r;
            warning('on')
        end
        
        function [source_stats, target_stats] = apply(sourceRoi, targetRoi, coeff)
            source_stats = sourceRoi * coeff.A;
            target_stats = targetRoi * coeff.B;
        end
    end
end