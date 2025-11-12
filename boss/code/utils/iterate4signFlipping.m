function [pval_reshaped, passing_centers, passing_indices, passing_matrix] = iterate4signFlipping(output_folder, ...
    suffix, corr_full, corr_res, nb, nb_signflip, centers_list_files)

    pval_group_file = fullfile(output_folder, "pval_group" + suffix);
    if exist(pval_group_file + '.mat', 'file')
        load(pval_group_file + '.mat', 'pval', 'maxstat', 'nb_signflip');
    else
        pval = {};
        maxstat = {};
    end
    for radius_indx = 1:nb.rad
        for ch_indx = 1:nb.channels
            nb_spheres = nb.spheres{radius_indx, ch_indx};
            x1 = single(corr_full{radius_indx, ch_indx});
            if isempty(x1)
                continue
            end

            for feat_set = 1:size(corr_res, 3)
                if all([radius_indx, ch_indx, feat_set] <= size(pval, 1:3)) && ~isempty(pval{radius_indx, ch_indx, feat_set})
                    continue
                end
                x2 = single(corr_res{radius_indx, ch_indx, feat_set});

                assert(size(x2, 1) == nb.subj && size(x1, 1) == nb.subj);
                assert(size(x2, 2) == nb_spheres && size(x1, 2) == nb_spheres);

                [pval{radius_indx, ch_indx, feat_set}, ...
                    maxstat{radius_indx, ch_indx, feat_set} ...
                    ] = compute_signflip_test( ...
                    reshape(x2, [size(x2, 1), prod(size(x2, [2,3]))]), ...
                    reshape(x1, [size(x1, 1), prod(size(x1, [2,3]))]), ...
                    nb_signflip ...
                    );

                save(pval_group_file, 'pval', 'maxstat', 'nb_signflip');
            end
        end
    end


    for radius_indx = 1:nb.rad
        for ch_indx = 1:nb.channels
            for feat_set = 1:nb.feat
    
                current_pval = pval{radius_indx, ch_indx, feat_set};
                if isempty(current_pval)
                    continue
                end
                pval_field_names = fieldnames(current_pval);
                corr_full_cell = corr_full{radius_indx, ch_indx};
    
                % Check if a sphere has any subjects or any component nan
                itisempty = any(isnan(corr_full_cell), [1, 3])';
    
                % Update center indices for non-empty cells
                centers_file = centers_list_files{radius_indx, ch_indx};
                load(centers_file, 'center_indices');
    
                % Reshape the 'fivecent' field of the struct
                for field = pval_field_names'
                    pval_reshaped = reshape(current_pval.(field{1}), ...
                        size(corr_full_cell, 2), ...
                        size(corr_full_cell, 3));
    
                    passing_indices_tmp = any(pval_reshaped, 2);
                    valid_passing_indices = passing_indices_tmp & ~itisempty;
                    passing_indices.(field{1}){radius_indx, ch_indx, feat_set} = valid_passing_indices;
                    passing_matrix.(field{1}){radius_indx, ch_indx, feat_set} = pval_reshaped;
                    passing_centers.(field{1}){radius_indx, ch_indx, feat_set} = center_indices(valid_passing_indices);
                end
    
            end
        end
    end
    
end
