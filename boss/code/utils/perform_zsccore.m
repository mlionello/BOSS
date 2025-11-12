function data = perform_zsccore(data)
    for j = 1:numel(data)
        for i = 1:size(data{j}, 2)
            for k = 1:size(data{j}, 3)
                data{j}(:, i, k) = zscore(data{j}(:, i, k));
            end
        end
    end
end
