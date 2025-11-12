function out_struct = cellstruct2structcell(in_cell)
out_struct = struct();
feat_names = fields(in_cell{1,1});
for i = 1: numel(feat_names)
    out_struct.(feat_names{i}) = cellfun(@(x) x.(feat_names{i}), ...
        in_cell, 'UniformOutput', false);
end
end
