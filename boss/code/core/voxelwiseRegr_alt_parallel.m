function r2 = voxelwiseRegr_alt_parallel(roi, infeat)
    nb_vx = size(roi,2);
    batch_size = 40000;
    nb_batches = nb_vx/batch_size;
    nb_batches = nb_batches + int32(mod(nb_vx,batch_size)>0);
    if isempty(gcp('nocreate'))
        parpool(min(maxNumCompThreads(), nb_batches));
    end
    r2 = cell(1, nb_batches); % or r2 = NaN(1, nb_vx);
    cr = "";
    parfor batch_i = 1:nb_batches
        if isempty(getCurrentJob)
            msg = compose(" OLS on batch %d/%d", batch_i, nb_batches);
            fprintf(cr + msg);
            cr = repmat('\b', 1, strlength(msg));
        end
        start_ind = (batch_i-1)*batch_size+1;
        end_ind = min(start_ind+batch_size-1, nb_vx);
        batch_indices = start_ind:end_ind;
        B = MultiProd(infeat', roi(:, batch_indices));
        Y = infeat*B;
        ssres = sum((Y - roi(:, batch_indices)).^2);
        sstot = sum((roi(:, batch_indices) - mean(roi(:, batch_indices))).^2);
        r2{batch_i} = 1 - ssres/sstot;
    end
    r2 = cat(2, r2{:});
    fprintf(cr);
end