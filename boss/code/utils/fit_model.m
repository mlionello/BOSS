function model = fit_model(average_subj_path, output_folder, target_model, init_model, opts)

    average_subj_path = average_subj_path + ".mat";
    fprintf(" using average subject from:\n\t%s\n", average_subj_path);

    %%% Model fitting on average subject
    if ~exist(fullfile(output_folder, 'model.mat'), 'file')
        if isfield(opts, 'update_features') && opts.update_features
            error("there is no model, but an update of feature list is required\n");
        end
        model = eval([target_model + "(init_model)"]);

        fprintf("Analysing (fitting) average subject:\n");
        msg = compose("\tloading avg image...");
        fprintf(msg);
        load(average_subj_path, 'nifti_image_avg');
        fprintf('Done\n');
        model = model.fit(nifti_image_avg, 'out_dir', output_folder);

        fprintf('saving model...')
        save(fullfile(output_folder, 'model'), "model", '-v7.3');
        fprintf('Done\n');
    else
        fprintf("\tLoading fitted model... ");
        load(fullfile(output_folder, 'model.mat'), 'model');
        fprintf("Done\n");
        if opts.update_params
            model = model.update_params(init_model);
        end
    end
end
