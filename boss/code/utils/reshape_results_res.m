function [a, nb_spheres] = reshape_results_res(x)
    nb_subj = size(x, 1);
    nb_rad = size(x, 3);
    nb_channels = size(x, 4);

    nb_feat = unique(cellfun(@(y) max(1, numel(y)), x(1,1,:,:)));
    
    a = cell(nb_rad, nb_channels, nb_feat);
    for i = 1 : nb_rad
        for j = 1 : nb_channels
            nb_spheres{i, j} = sum(cellfun(@(y) ~isempty(y), x(1,:,i,j)));

            tmp = cat(1, x{:, :, i, j});
            for k = 1 : nb_feat
                if ~isempty(tmp)
                    tmp2=cat(1, tmp{:, k});
                    if ndims(tmp2) == 3
                        tmp2 = squeeze(mean(tmp2, 2)); %mean across repetitions for permutations
                    end
                    nb_comps = size(tmp2, 2); %in the case no components dim the second dim is for rept
                    a{i, j, k} = reshape(tmp2, [nb_subj, nb_spheres{i, j}, nb_comps]);
                    if ndims(x{3,5,i,j}{k}) == 3
                        assert(all(squeeze(a{i, j, k}(3,5,:)) == squeeze(mean(x{3,5,i,j}{k}, 2))))
                    else
                        assert(all(squeeze(a{i, j, k}(3,5,:)) == squeeze(x{3,5,i,j}{k}(1,:)'), 'all'))
                    end
                    
                else
                    a{i, j, k} = zeros(nb_subj, nb_spheres{i, j}, 0);
                end
            end
        end
    end
end
