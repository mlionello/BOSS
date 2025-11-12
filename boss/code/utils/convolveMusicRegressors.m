function musicRegr_init_conv = convolveMusicRegressors(data, hrf)
    musicRegr_init_conv = {};
    for feat_set = 1:numel(data)
        for ch = [1, 2]
            for inner_feat = 1:size(data{feat_set}, 2)
                tmp_data = data{feat_set}(:, inner_feat, ch);
                conv_data = conv(tmp_data, hrf);
                if isempty(musicRegr_init_conv)
                    musicRegr_init_conv{feat_set} = nan(size(conv_data, 1), size(data{feat_set}, 2), 2);
                end
                musicRegr_init_conv{feat_set}(:, inner_feat, ch) = conv(tmp_data, hrf);
            end
        end
    end
end