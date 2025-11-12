classdef VggishExtractor < AudioFeatureExtractor
    properties
        % featureNames = ["pool4", "fc1_1", "relu5_1",  "fc1_2", "relu5_2", "EmbeddingBatch"];
        featureNames = ["fc1_1",  "fc1_2", "EmbeddingBatch"];
        vggishNet;
    end

    methods
        function obj = VggishExtractor(audioPath, stimuliOrderFiles)
            obj@AudioFeatureExtractor(audioPath, stimuliOrderFiles);

            VGGishLocation = '../../data/stimuli/';
            if ~exist(fullfile(VGGishLocation,"vggish"), 'dir')
                zip_file = fullfile(VGGishLocation,"vggish.zip");
                if ~exist(zip_file, 'file')
                    websave(zip_file, "https://ssd.mathworks.com/supportfiles/audio/vggish.zip");
                end
                unzip(zip_file, VGGishLocation)
            end
            addpath(fullfile(VGGishLocation, "vggish"))
            obj.vggishNet = vggish; %audioPretrainedNetwork("vggish");
        end
        
        function [embedings, num_sec] = FPU(obj, audioFile)
            [audio, sr] = audioread(audioFile);
            num_sec = ceil(size(audio, 1) / sr);

            features = vggishPreprocess(audio, sr, OverlapPercentage = 66);
            embedings = obj.predict_vggish(features, size(audio, 2));
        end

        function embedings = predict_vggish(obj, features, nb_channels)
            embedings = cell(numel(obj.featureNames), 1);
            for l_i = 1:numel(obj.featureNames)
                %obj.vggishNet.OutputNames = obj.featureNames(l_i);
                outfeat = activations(obj.vggishNet,features,obj.featureNames(l_i));
                %outfeat = predict(obj.vggishNet, features);
                if ndims(outfeat) > 2
                    outfeat = reshape(outfeat, [], size(outfeat, ndims(outfeat)));
                else
                    outfeat = permute(outfeat, [2 1]);
                end
                outfeat = reshape(outfeat, size(outfeat, 1), size(outfeat, 2) / nb_channels, nb_channels);
        
                embedings{l_i} = permute(outfeat, [2 1 3]);
            end
        
        end

 
    end
end
