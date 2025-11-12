classdef perform_full_and_res_corr_on_searchlight

    methods (Static)

        function results = corr(x, y, mdata)
            % perform canonical between the two ROIs and
            results.full_corr = corr(mean(x, 2), mean(y, 2));
            for feat_i = 1:numel(mdata.feature_regressors)
                feat = mdata.feature_regressors{feat_i};
                results.res_corr(feat_i) = corr( ...
                    get_residual(mean(x, 2), feat), ...
                    get_residual(mean(y, 2), feat));
    
            end

        end
        
    end
end