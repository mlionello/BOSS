function flatten = do_flat_cells(in_cells)

cls_id = {};
if any(cellfun(@(x) isstruct(x), in_cells{1}))
    itisstruct = cellfun(@(x) isstruct(x), in_cells{1});
    cls_id = fieldnames(in_cells{1}{find(itisstruct,1)});
    flatten = struct();
    for idf = 1:numel(cls_id)
        tmp_cell = in_cells{1}(:);
        max_comps = nan;
        for region = 1 : numel(tmp_cell)
            if ~isempty(tmp_cell{region})
                if any(isnan(max_comps))
                    max_comps = size(tmp_cell{region}.(cls_id{idf}));
                else
                    max_comps = max(max_comps, size(tmp_cell{region}.(cls_id{idf})));
                end
            end
        end
        flatten.(cls_id{idf}) = nan([numel(in_cells), numel(in_cells{1}), ...
            max_comps]);
    end
else
    max_comps = max(cellfun(@(x) length(x), in_cells{1}));
    flatten = nan([numel(in_cells), numel(in_cells{1}), max_comps]);
end

for i = 1:numel(in_cells)
    for j = 1:numel(in_cells{i})
        if isempty(cls_id)
            flatten(i, j, 1:length(in_cells{i}{j}), :, :) = in_cells{i}{j};
        else
            if isempty(in_cells{i}{j})
                continue
            end
            for idf = 1:numel(cls_id)
                values = in_cells{i}{j}.(cls_id{idf});
                flatten.(cls_id{idf})(i, j, 1:size(values, 1), 1:size(values, 2), 1:size(values, 3)) = values;
            end
        end
    end
end
end
