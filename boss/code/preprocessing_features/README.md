### FEATURES EXTRACTION AND PREPROCESSING ####

## Essentia_extractor
python Essentia_extractor.py creates soundFeatDict.json, musicFeatDict.json(not working) and emoFeatDict.json
Essentia_extractor.py already arange the songs features in the order according to the *all_runscat.1D file (16 songs times 2 repetitions)

By reading the json with jsondecode(fileread('...'));
    N fields one for each feature,
        32 by 1 cells one for each (ordered) songs,
            timepoints by {1,2} channels

They include also the following fields {'onset_global', 'offset_global', 'time_end', 'time_init'} - 32 by 1 cells
THEY ARE NOT SAMPLED AT THE RIGHT SCALE

## Behavioural extractor
python behaviouralFeatureExtraction.py generates behavFeatDict.json and behavFeatDict_subj.mat
it already arrange the sliders data in the song order according to the *all_runscat.1D (16 songs times 2 rep)

behavFeatDict.json:
By reading the json with jsondecode(fileread('...'));
    N fields one for each slider,
        32 by 1 cells one for each (ordered) songs,
            timepoints by 1: timepoints are supposed ot have the same resolution of the acquisition mehtod (10Hz)

behavFeatDict_subj.mat
contains a matrix  nb_features by nb_participants by nb_songs by nb_max_timepoints (sampled at 10Hz)

behav_preprocessing.m will read behavFeatDict_subj.mat and generate behavFeatDict_avg.mat
it downsamples to match 1Hz by averaging (no overlapping) and it computes subject average, standard deviation and feature contrasts.
behavFeatDict_avg.mat contains 3 fields (avg, std, contrast)
    1 by 32 cells
        timepoints (in seconds) by 10 features (sliders) by 2 channels (identical)


## ISMIR and VGGISH extractor
### ISMIR

ismirExtractor = MirExtractor(audioPath, stimuliOrderFiles);
[mir_raw, mir_avg, mir_std, mir_std_long, durations] = ismirExtractor.extract_features();
[mir_pca_coeff{1}, mir_pca{1}, mir_pca_latent{1}, mir_pca_expl{1}] = ismirExtractor.get_pca_between_fields(mir_and_ess, data_global);

###VGGISH
vggishExtractor = VggishExtractor(audioPath, stimuliOrderFiles);
[vggish_raw, vggish_avg, vggish_std, vggish_std_long, durations] = vggishExtractor.extract_features();
featureNames = vggishExtractor.featureNames;
for k = 1:numel(vggishExtractor.featureNames)
    field = vggishExtractor.featureNames{k};
    [vgg_pca_coeff{k}, vgg_pca{k}, vgg_pca_latent{k}, vgg_pca_expl{k}] = vggishExtractor.get_pca_within_fields(vggish_avg.(field), data_global);
end






