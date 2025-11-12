function results = resCell2Mat(cca_values, winkler)
    cr = "";
    for j = 1:size(cca_values.r, 1)
        msg = compose("transforming to struct array row %d out of %d ...", j, size(cca_values.r, 1));
        fprintf(cr + msg);
        cr = repmat('\b', 1, strlength(msg));
        pause(0.0001)

        results.rMusVsSource(j, :) = cca_values.rMusSource{j};
        results.sourceRoiCoord(j) = cca_values.sourceRoiCenter{j};
        for i = 1:size(cca_values.r, 2)
            results.rMusVsTarget(i, :) = cca_values.rMusTarget{i};
            results.targetRoiCoord(i, :) = cca_values.targetRoiCenter{1, i};
            results.rSourceVsTarget(j, i) = cca_values.r{j, i}(1);
            results.pSourceVsTarget(j, i) = cca_values.p{j, i}(1);
            results.resCorrPerm(j, i, :, :) = cca_values.resCorrPerm{j, i};
            results.R2U(j, i, :, :) = cca_values.R2U{j, i};
            results.R2V(j, i, :, :) = cca_values.R2V{j, i};
        end
        if ~winkler
            for i = 1:size(cca_values.r, 2)
                results.wilks(j, i) = cca_values.wilks{j, i}(1);
                results.pF(j, i) = cca_values.pF{j, i}(1);
                results.F(j, i) = cca_values.F{j, i}(1);
            end
        end
    end
    fprintf(' Done\n');
end
