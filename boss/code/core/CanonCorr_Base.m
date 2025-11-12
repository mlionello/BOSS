classdef CanonCorr_Base < CanonCorr_abs
    methods
        function obj = CanonCorr_Base(opts)
            arguments
                opts (1,1) struct
            end

            fields = {'output_folder', 'sourceRegionMask', 'targetRegionMask', ...
                      'feature_regressors', 'feature_names', 'winkler', 'nb_perm', 'nb_perm'};

            init = struct('winkler', 0, 'nb_perm', 0, 'extra', 0);

            obj = parse_args(obj, fields, init, opts);

        end

        function obj = parse_args(obj, fields, init, opts)
            for i = 1:numel(fields)
                field = fields{i};
                if isfield(opts, field)
                    obj.(field) = opts.(field);
                else
                    if isfield(init, field)
                        obj.(field) = init.(field);
                    else
                        error("input field %s needs a value" ,field);
                    end
                end
            end
        end

    end
end
