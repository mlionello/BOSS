function [AVG_targetsource, AVG_featsource, AVG_feattarget] = fit_on_average_subj(output_folder, average_subj_path, sourceRegionMask, targetRegionMask, feat_all, winkler, extra, nb_perm, fitfun)

    if ~gcp().Connected
        parpool(maxNumCompThreads());
    end
    fprintf("Analysing average subject: ");
    msg = compose(" loading image...");
    cr1 = repmat('\b', 1, strlength(msg));
    fprintf(msg);
    nifti_image_avg = niftiread(average_subj_path);
    model.fit(nifti_image_avg);
    fprintf(cr1);


end
