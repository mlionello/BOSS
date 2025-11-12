function r2 = voxelwiseRegr_alt(roi, infeat)
    nb_vx = size(roi,2);
    batch_size = 100;
    nb_batches = nb_vx/batch_size;
    nb_batches = nb_batches + int32(mod(nb_vx, batch_size)>0);
    r2 = zeros(nb_vx, 1); % or r2 = NaN(1, nb_vx);
    cr = "";
    for batch_i = 1 : nb_batches
        if true || isempty(getCurrentJob)
            msg = compose(" OLS on batch %d/%d", batch_i, nb_batches);
            fprintf(cr + msg);
            cr = repmat('\b', 1, strlength(msg));
        end
        start_ind = (batch_i-1)*batch_size+1;
        end_ind = min(start_ind+batch_size-1, nb_vx);
        batch_indices = start_ind:end_ind;
        B = SliceMultiSolver(infeat(:, :), roi(:, batch_indices));
        Y = infeat * B;
        ssres = sum((Y - roi(:, batch_indices)).^2);
        sstot = sum((roi(:, batch_indices) - mean(roi(:, batch_indices))).^2);
        r2(batch_indices) = squeeze(1 - ssres./(sstot+eps))';
        r2(batch_indices(sstot==0)) = 0;
    end
    fprintf(cr);
end
