function [corr_full, corr_res, nb] = reshape_corr_res(corr_full, corr_res)

    nb.subj = size(corr_full, 1);
    nb.rad = size(corr_full, 3);
    nb.channels = size(corr_full, 4);
    nb.feat = unique(cellfun(@(y) max(1, numel(y)), corr_res(1,1,:,:)));
    
    [corr_full, nb_spheres0] = reshape_results_full(corr_full);
    [corr_res, nb_spheres] = reshape_results_res(corr_res);
    %assert(all(cellfun(@(x,y) x==y == nb_spheres))
    
    nb.spheres = nb_spheres;

    assert(ndims(corr_full) == 2 && ndims(corr_res) == (2 + int32(nb.feat>1)));
    assert(all(size(corr_full) == [nb.rad, nb.channels]));
    assert(all(size(corr_res, [1, 2]) == [nb.rad, nb.channels]));
    assert(size(corr_res, 3) ==  nb.feat);
    
    %cellfun(@(x) assert(all(size(x, [1, 2])==[nb.subj, nb.spheres])), corr_full);
    %cellfun(@(x) isempty(x) || assert(all(size(x, [1, 2])==[nb.subj, nb.spheres])), corr_res);
end
