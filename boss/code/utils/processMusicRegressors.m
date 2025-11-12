function [musicRegrPCA_conv, coeffPCA, musicRegr_init_conv, labels_init] = processMusicRegressors(musicRegr_init, labels_init, outfolder)
% Function: processMusicRegressors
%
% Description:
%   Processes the input music regressors by convolving with the hemodynamic response function (HRF),
%   removing correlated columns, and performing principal component analysis (PCA) with z-score normalization.
%
% Inputs:
%   musicRegr_init  - Initial music regressors matrix.
%   labels_init     - Initial labels corresponding to the regressors.
%   corr_thr        - Correlation threshold for removing correlated columns.
%
% Outputs:
%   musicRegr_crosscorr  - Music regressors matrix after removing correlated columns.
%   labels               - Updated labels after removing correlated columns.
%   musicRegrPCA_conv    - Convolved music regressors matrix after PCA.
%   coeffPCA             - Principal component coefficients.
%   musicRegr_init_conv  - Convolved initial music regressors matrix.
%   labels_init          - Updated labels after replacing underscores with '\_'.
%
% Example:
%   corr_thr = 0.7;
%   [musicRegr_crosscorr, labels, musicRegrPCA_conv, coeffPCA, musicRegr_init_conv, labels_init] = ...
%       processMusicRegressors(musicRegr_init, labels_init, corr_thr);
labels_init = arrayfun(@(x) replace(x, '_', "\_"), labels_init);

TR = 1;
hrf = spm_hrf(TR);
musicRegr_init_conv = convolveMusicRegressors(musicRegr_init, hrf);

%[musicRegr_crosscorr, labels] = removeCorrelatedColumns(musicRegr_init_conv, labels_init, corr_thr);

[coeffPCA_init, musicRegrPCA_init, latent_init] = pca(musicRegr_init);

% Apply Varimax rotation to the PCA coefficients
[coeffPCA, T] = rotatefactors(coeffPCA_init(:,1:10), 'method', 'varimax');

% Transform the PCA scores using the rotation matrix
musicRegrPCA = musicRegrPCA_init(:,1:10) * T;

% Recalculate the explained variance for the rotated components
latent = var(coeffPCA, 1);

musicRegrPCA_conv = convolveMusicRegressors(musicRegrPCA, hrf);

outfolder = fullfile(outfolder, "figures");
if ~exist(outfolder, 'dir'); mkdir(outfolder); end
laodings2table(coeffPCA, labels_init, outfolder);
laodings2graphics(coeffPCA, labels_init, outfolder);
loading2plot(coeffPCA, labels_init, outfolder);
loadings2hist(coeffPCA, labels_init, outfolder);
savelatents(latent, outfolder)
end

function savelatents(latent, outfolder)
    figure('Visible', 'off');  % Create the figure without displaying it
    bar(latent/sum(latent)*100);
    xlabel('Principal Components');
    ylabel('% variance explained');
    title('Latents Bar Plot');
    grid on;
    saveas(gcf, fullfile( outfolder, 'variance_explained.png'));
end


function laodings2table(coeff, labels, outfolder)
    % Number of components
    num_components = size(coeff, 2);
    % Create a formatted table
    formatted_coeff = nan(size(coeff));
    % Loop through each element in the coeff matrix
    for i = 1:size(coeff, 1)
        for j = 1:size(coeff, 2)
            if abs(coeff(i, j)) < 0.1
                formatted_coeff(i, j) = nan;
            else
                formatted_coeff(i, j) = coeff(i, j);
            end
        end
    end
    % Convert formatted_coeff to a table and add feature names as the first column
    T = array2table(round(formatted_coeff,2), 'VariableNames', strcat('PC', string(1:num_components)));
    T.Features = labels';
    % Move the Features column to the front
    T = [T(:, end) T(:, 1:end-1)];
    
    % Save the table to a CSV file
    writetable(T, fullfile( outfolder, 'PCA_Loadings.csv'));
end

function laodings2graphics(coeff, labels, outfolder)
    % Number of components and features
    [num_features, num_components] = size(coeff);
    
    % Create a grouped bar plot
    figure('Visible', 'off');  % Create the figure without displaying it
    bar(coeff);
    set(gca, 'XTickLabel', labels, 'XTick', 1:num_features);
    xlabel('Features');
    ylabel('Loadings');
    legend(arrayfun(@(x) ['PC ', num2str(x)], 1:num_components, 'UniformOutput', false));
    title('Loadings of Features on Principal Components');
    xtickangle(45);
    grid on;
    
    % Save the figure as a PNG file
    saveas(gcf, fullfile( outfolder, 'pca_loadings.png'));
end

function loadings2hist(coeff, labels, outfolder)
    % Number of components
    num_components = size(coeff, 2);
    
    % Create a bar plot for each principal component
    for i = 1:num_components
        figure('Visible', 'off'); % Create a figure without displaying it
        bar(coeff(:, i));
        set(gca, 'XTickLabel', labels, 'XTick', 1:numel(labels));
        xlabel('Features');
        ylabel('Loadings');
        title(['Principal Component ', num2str(i)]);
        xtickangle(45);
        grid on;
        
        % Save the figure as a PNG file
        saveas(gcf,  fullfile( outfolder, ['Principal_Component_' num2str(i) '.png']));
        close(gcf); % Close the figure after saving
    end
end

function loading2plot(coeff, labels, outfolder)
    compcomb = [1 2 3 1];
    for k = 1 : numel(compcomb) - 1
        figure('Visible', 'off'); % Create a figure without displaying it
        hold on;
        
        % Plot the unit circle
        theta = linspace(0, 2*pi, 100);
        x = cos(theta);
        y = sin(theta);
        plot(x, y, 'k--'); % Unit circle
        
        % Plot the loadings
        quiver(zeros(size(coeff, 1), 1), zeros(size(coeff, 1), 1), ...
            coeff(:, compcomb(k)), coeff(:,  compcomb(k + 1)), 0, 'r');
        
        % Annotate the plot with feature names
        for i = 1:length(labels)
            text(coeff(i, compcomb(k)), coeff(i, compcomb(k + 1)), labels{i}, 'FontSize', 8, 'FontWeight', 'bold');
        end
        
        % Set plot limits and labels
        xlim([-0.5, 0.5]);
        ylim([-0.5, 0.5]);
        xlabel("Principal Component " + compcomb(k));
        ylabel("Principal Component " + compcomb(k + 1));
        title('Loading Plot (Circle of Loadings)');
        grid on;
        axis equal;
        
        % Save the figure as a PNG file
        saveas(gcf,  fullfile( outfolder, ['Circle_of_Loadings_' num2str(k) '.png']));
        close(gcf); % Close the figure after saving
    end
end