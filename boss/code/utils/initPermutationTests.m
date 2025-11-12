function permSchema = initPermutationTests(nb_perm, nb_tpoints)
    permSchema(:, 1) = 1:nb_tpoints;
    for i = 1:nb_perm
        permSchema(:, i + 1) = randperm(nb_tpoints);
    end
end
