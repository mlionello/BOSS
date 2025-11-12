function resultStruct = cellArrayToStruct(cellArray)
    % Initialize result struct
    resultStruct = struct();

    % Extract field names
    fieldNames = fieldnames(cellArray{1});

    % Loop through each field
    for i = 1:numel(fieldNames)
        % Initialize cell array for values
        values = cell(size(cellArray));

        % Loop through each structure in the cell array
        for j = 1:size(cellArray, 1)
            for k = 1:size(cellArray, 2)
                % Extract field value from each structure
                values{j, k} = cellArray{j, k}.(fieldNames{i});
            end
        end

        % Assign the cell array to the result struct
        resultStruct.(fieldNames{i}) = values;
    end
end
