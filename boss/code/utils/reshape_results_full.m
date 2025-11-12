function [a, nb_spheres] = reshape_results_full(x)
    nb_subj = size(x, 1);
    nb_rad = size(x, 3);
    nb_channels = size(x, 4);
    
    a = cell(nb_rad, nb_channels);
    for i = 1 : nb_rad
        for j = 1 : nb_channels
            nb_spheres{i, j} = sum(cellfun(@(y) ~isempty(y), x(1,:,i,j)));

            tmp = cat(1, x{:, :, i, j});
            nb_comps = size(tmp, 2);
            a{i, j} = reshape(tmp, [nb_subj, nb_spheres{i, j}, nb_comps]);
            if ~isempty(a{i, j})
                assert(all(squeeze(a{i,j}(3,5,:)) == squeeze(x{3,5,i,j}(1,:)'), 'all'));
            end
        end
    end
end
