classdef AudioFeatureExtractor
    properties
        audioPath
        stimuliOrderFiles
        listAudioFiles
        stimuliOrderAllRun
    end
    
    properties (Abstract)
        featureNames
    end
    
    methods
        function obj = AudioFeatureExtractor(audioPath, stimuliOrderFiles)
            obj.audioPath = audioPath;
            obj.stimuliOrderFiles = stimuliOrderFiles;
            obj.listAudioFiles = obj.get_audio_files();
            obj.stimuliOrderAllRun = obj.load_stimuli_orders();
        end
        
        function audioFiles = get_audio_files(obj)
            files = dir(fullfile(obj.audioPath, '*.wav'));
            audioFiles = fullfile({files.folder}, {files.name});
            audioFiles = sort(audioFiles);
        end
        
        function stimuliOrderAllRun = load_stimuli_orders(obj)
            stimuliOrderAllRun = [];
            for i = 1:length(obj.stimuliOrderFiles)
                fid = fopen(obj.stimuliOrderFiles{i}, 'r');
                row_tmp = textscan(fid, '%n ', 'delimiter', ' ');
                stimuliOrderAllRun = [stimuliOrderAllRun, row_tmp{1}'];
            end
            stimuliOrderAllRun = stimuliOrderAllRun + 1;
        end

        function [out_raw, out_avg, out_std, out_std_long, durations] = extract_features(obj)
            [out_raw, out_avg, out_std, out_std_long, durations] = obj.iterate_audio_files();
        end
    end

    methods(Hidden)    
        function [out_raw, out_avg, out_std, out_std_long, durations] = iterate_audio_files(obj)
            out_raw = struct();
            out_avg = struct();
            out_std = struct();
            out_std_long = struct();
            durations = nan(1, length(obj.stimuliOrderAllRun));
            for j = 1 : length(obj.stimuliOrderAllRun)
                stimuli_i = obj.stimuliOrderAllRun(j);
        
                if j > 1
                    prev_run_ind = find(obj.stimuliOrderAllRun(1 : j-1) == stimuli_i);
                else
                    prev_run_ind = [];
                end
        
                if isempty(prev_run_ind)
                    [out_features, num_sec] = obj.FPU(obj.listAudioFiles{stimuli_i});
                    durations(j) = num_sec;
        
                    for k = 1:numel(obj.featureNames)
                        fieldName = obj.featureNames{k};
                        if ~isfield(out_raw, fieldName)
                            out_raw.(fieldName) = {};
                            out_avg.(fieldName) = {};
                            out_std.(fieldName) = {};
                            out_std_long.(fieldName) = {};
                        end
                        out_raw.(fieldName){end+1} = out_features{k};
                        [out_avg.(fieldName){end+1}, ...
                            out_std.(fieldName){end+1}, ...
                            out_std_long.(fieldName){end+1}] = obj.get_avg_per_sec(out_features{k}, num_sec);
                    end
            
                else
                    durations(j) = durations(prev_run_ind(1));
                    for k = 1:numel(obj.featureNames)
                        fieldName = obj.featureNames{k};
                        out_raw.(fieldName){end+1} = out_raw.(fieldName){prev_run_ind(1)};
                        out_avg.(fieldName){end+1} = out_avg.(fieldName){prev_run_ind(1)};
                        out_std.(fieldName){end+1} = out_std.(fieldName){prev_run_ind(1)};
                        out_std_long.(fieldName){end+1} = out_std_long.(fieldName){prev_run_ind(1)};
                    end
                end
            end
        end
    end
    
    methods (Abstract)
        [out_features, num_sec] = FPU(obj, audioFile);
    end
    
    methods (Static)
        function [avgFeature, stdevFeature, stdevFeature_long] = get_avg_per_sec(features, numSeconds)
            frameOnset = 0 : numSeconds / size(features, 1) : numSeconds;
            if size(features, 1) < numel(frameOnset)
                numRepeat = numel(frameOnset) - size(features, 1);
                features = cat(1, features, repmat(features(end, :, :), numRepeat, 1, 1));
            end
            features = features(1 : numel(frameOnset), :, :);
            avgFeature = zeros(numSeconds, size(features, 2), size(features, 3));
            stdevFeature = zeros(numSeconds, size(features, 2), size(features, 3));
            stdevFeature_long = zeros(numSeconds, size(features, 2), size(features, 3));

            for sec = 0 : numSeconds
                featIndices = floor(frameOnset) == sec;
                featIndices = featIndices(1 : min(numel(featIndices), size(features, 1)));
                avgFeature(sec + 1, :, :) = mean(features(featIndices, : ,:), 1);
                if sum(featIndices) == 1
                    stdevFeature(sec + 1, :, :) = 0;
                else
                    stdevFeature(sec + 1, :, :) = std(features(featIndices, : ,:), 1);
                end

                featIndices_long = floor(frameOnset) >= sec-2 & floor(frameOnset) <= sec;
                featIndices_long = featIndices_long(1 : min(numel(featIndices_long), size(features, 1)));
                if sum(featIndices_long) == 1
                    stdevFeature_long(sec + 1, :, :) = 0;
                else
                    stdevFeature_long(sec + 1, :, :) = std(features(featIndices_long, : ,:), 1);
                end
            end
        end
        
        function [coeffCh, scoreCh, latentCh, explainedCh] = get_pca(data, global_onsets, opts)
            arguments
                data cell;
                global_onsets;
                opts.compute_individual = 0;
            end
            assert(numel(data) == size(global_onsets, 1))

            explainedCh = [];
            scoreCh = [];
            coeffCh = [];
            latentCh = [];
            coeffSong = {};
            nbChannels = size(data{1}, 3);
            if opts.compute_individual
                coeffSong = cell(size(data));
            end
            concatData = nan([round(global_onsets(32, 2)), size(data{1}, [2, 3])]);
            for stimulus = 1 : numel(data)
                start_stim = round(global_onsets(stimulus, 1)) + 1;
                end_stim = round(global_onsets(stimulus, 2));
                dur_stim = end_stim-start_stim+1;
                concatData(start_stim:end_stim, :, :) = data{stimulus}(1 : round(dur_stim), :, :);
                slice_tmp = data{stimulus}(1 : round(dur_stim), :, :);
                z = find(any(isnan(slice_tmp), [2,3]));
            end
            missing_indices = find(any(isnan(concatData), [2,3]));
            for ind = 1: numel(missing_indices)
                concatData(missing_indices, :, :)= concatData(missing_indices-1, :, :);
            end

            for c = 1 : nbChannels
                Xconcat = concatData(:, :, c);
                [coeff, scoreConcat, latent, ~, explained] = pca(Xconcat(:, ~(sum(isnan(Xconcat), 1) > 0.7 * size(Xconcat, 1))));
                if isempty(explainedCh)
                    explainedCh = explained;
                    scoreCh = scoreConcat;
                    coeffCh = coeff;
                    latentCh = latent;
                else
                    explainedCh(:, 2) = explained;
                    scoreCh(:, :, 2) = scoreConcat;
                    coeffCh(:, :, 2) = coeff;
                    latentCh(:, :, 2) = latent;
                end
                if opts.compute_individual
                    for j = 1:numel(data)
                        X = data{j}(:, :, c);
                        coeffSong{j}(:, :, c) = X * coeff;
                    end
                end
            end
        end
        
        function [coeffCh, scoreCh, latentCh, explainedCh] = get_pca_between_fields(data, global_onsets, opts)
            arguments
                data struct;
                global_onsets;
                opts.compute_individual = 0;
            end
            fields = fieldnames(data);
            egCell = data.(fields{1});
            for i = 1:numel(egCell)
                temp = cell(1, numel(fields));
                for j = 1:numel(fields)
                    row = data.(fields{j}){i};
                    row = row(:, 1: min(16,size(row, 2)), :);
                    if size(row, 3) == 1
                        row = repmat(row, 1, 1, 2);
                    end
                    temp{j} = row;
                end
                concatenatedResults{i, 1} = cat(2, temp{:, :, 1});
            end
            [coeffCh, scoreCh, latentCh, explainedCh] = AudioFeatureExtractor.get_pca(concatenatedResults, global_onsets, 'compute_individual', opts.compute_individual);
        end
    end
end
