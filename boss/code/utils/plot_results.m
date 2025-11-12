function plot_results(a, b, percent_drop, corr_res_fs_feat)
    N = size(a, 2); C = repmat(linspace(0.3, 0.7, N).', 1, 3);
    C(1 : 3, :) = [[0 1 1]; [1 0 1]; [1 1 0]];
    figure; imagesc(percent_drop)
    figure; axes('ColorOrder', C,'NextPlot','replacechildren'); plot(percent_drop, 'LineWidth', 0.8);
    hold on; plot(percent_drop(:, 1), 'LineWidth', 0.8, Color='red')
    plot(percent_drop(:, 2), 'LineWidth', 0.8, Color='yellow')
    plot(percent_drop(:, 3), 'LineWidth', 0.8, Color='blue')
    plot((1 - corr_res_fs_feat' ./ corr_full_passing') * 100, 'LineWidth', 1, Color='black', LineStyle=':')
    title("PC plots")
    xlabel("Vx region passing singificance thr")
    ylabel("correlation drop %")
    
    figure; plot((1 - b' ./ a') * 100);
    title("vx region curves")
    xlabel("Principal Components")
    ylabel("correlation drop %")
    
    figure; plot(a - b);
    figure; plot(a' - b');

end