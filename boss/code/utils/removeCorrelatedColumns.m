function [musicRegr, labels] = removeCorrelatedColumns(musicRegr_init_conv, labels_init, corr_thr)
    [maxvalue, rm_feats] = max(sum((corrcoef(musicRegr_init_conv) - eye(size(musicRegr_init_conv, 2))) > corr_thr, 1));
    fprintf("Removing columns: ");
    labels = labels_init;
    musicRegr = musicRegr_init_conv;
    while maxvalue > 0
        fprintf("%s; ", labels(rm_feats));
        musicRegr(:, rm_feats) = [];
        labels(rm_feats) = [];
        [maxvalue, rm_feats] = max(sum((corrcoef(musicRegr) - eye(size(musicRegr, 2))) > corr_thr, 1));
    end

end
