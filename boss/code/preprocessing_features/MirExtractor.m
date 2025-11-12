classdef MirExtractor < AudioFeatureExtractor
    properties
        featureNames = {'zero_crossing_rate', 'rms', 'centroid', 'spread', 'rolloff', 'entropy', ...
                        'flatness', 'brightness', 'roughness', 'flux', 'sub_band_flux', ...
                        'irregularity', 'inharmonicity', 'fluctuation', 'fluctuation_peak', ...
                        'fluctuation_centroid_entropy'};
        funcsNames = {@mirzerocross, @mirrms, @mircentroid, @mirspread, @mirrolloff, @mirentropy, ...
                       @mirflatness, @mirbrightness, @mirroughness, @mirflux, @mirflux, ...
                       @mirregularity, @mirinharmonicity, @mirfluctuation, @mirpeaks, @mircentroid};
    end

    methods
        function obj = MirExtractor(audioPath, stimuliOrderFiles)
            obj@AudioFeatureExtractor(audioPath, stimuliOrderFiles);
        end
        
        function [out_features, num_sec] = FPU(obj, audioFile)
            audio = miraudio(audioFile, 'Mono', 0);
            num_sec = ceil(mirgetdata(mirlength(audio)));
            out_features = obj.mir_tool(audio, obj.funcsNames, obj.featureNames);
        end
        
        function rawFeatures = mir_tool(~, audio, funcs, featureNames)
             % Initialize feature extraction parameters
            ShortFrameLength = 0.025;
            Hop_short_perc = 50;
            LongFrameLength = 3.0;
            Hop_Long_seconds = 1;

            frames_short = mirframe(audio, 'Length', ShortFrameLength, 's', 'Hop', Hop_short_perc, '%');
            audio_spectrum_short = mirspectrum(frames_short);
            peaks = mirpeaks(audio_spectrum_short, 'Contrast', 0.01, 'Order', 'Abscissa');
            mono_audio = miraudio(audio, 'Mono', 1);
            frames_short_mono = mirframe(mono_audio, 'Length', ShortFrameLength, 's', 'Hop', Hop_short_perc, '%');
            audio_spectrum_short_mono = mirspectrum(frames_short_mono);

            rawFeatures = cell(length(featureNames), 1);
            fluctuationData = [];
            for i = 1:length(featureNames)
                featureName = featureNames{i};
                func = funcs{i};
                switch featureName
                    case 'sub_band_flux'
                        featureData = func(audio_spectrum_short, 'SubBand');
                    case 'rms'
                        featureData = func(frames_short);
                    case 'fluctuation'
                        fluctuationData = func(mono_audio, 'Frame', 'Summary');
                        featureData = fluctuationData;
                    case {'fluctuation_peak'}
                        if isempty(fluctuationData)
                            error('Fluctuation data is required for %s', featureName);
                        end
                        featureData = func(fluctuationData, 'Total', 1);
                    case {'fluctuation_centroid_entropy'}
                        if isempty(fluctuationData)
                            error('Fluctuation data is required for %s', featureName);
                        end
                        featureData = func(fluctuationData,'MaxEntropy',0);
                    case {'mirregularity', 'roughness'}
                        featureData = func(peaks);
                    case 'inharmonicity'
                        featureData = func(audio_spectrum_short_mono);
                    case {'pulse_clarity', 'mode_clarity'}
                        featureData = func(mono_audio, 'Frame', LongFrameLength, 's', Hop_Long_seconds, 's');
                    otherwise
                        featureData = func(audio_spectrum_short);
                end

                featureData = mirgetdata(featureData);
                if ndims(featureData)<3 && ismember(size(featureData, 2), [1,2])
                    featureData = reshape(featureData, 1, size(featureData,1), size(featureData,2));
                end
                
                rawFeatures{i} = permute(featureData, [2 1 3]);

            end
        end
        

    end
end
