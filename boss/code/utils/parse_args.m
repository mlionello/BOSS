% function obj = parse_args(obj, fields, init, opts)
%     for i = 1:numel(fields)
%         field = fields{i};
%         if isfield(opts, field)
%             obj.(field) = opts.(field);
%         else
%             if isfield(init, field)
%                 obj.(field) = init.(field);
%             else
%                 error("input field %s needs a value" ,field);
%             end
%         end
%     end
% end
