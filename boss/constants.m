global pkg_folder data_path deriv_path atlas_path notebook_path ...
    stimuli_path res_path file_pattern cort_atlas_path cort_labels ...
    subcort_atlas_path subcort_labels feature_list meansubject_path ...
    commonmask_path HCP_atlas_path HCP_labels deriv_path_orig features_path;

pkg_folder = fileparts(mfilename('fullpath'));
data_path = fullfile(pkg_folder, 'data');
deriv_path = "/home/matteo.lionello/Code/CCABA_DVL/ccaba/data/derivatives/";
% deriv_path = fullfile(data_path, 'shorts');
atlas_path = fullfile(data_path, 'atlas');
meansubject_path = fullfile(data_path, 'AVG_EPI');
commonmask_path = fullfile(data_path, 'COMMON_MASK');
notebook_path = fullfile(pkg_folder, 'notebook');
stimuli_path = fullfile(data_path, 'stimuli');
% stimuli_path = fullfile(pkg_folder, '../../../CCABA_DVL/ccaba/data/stimuli');
% stimuli_path = fullfile(data_path, 'stimulifake');
features_path = fullfile(stimuli_path, 'features');
res_path = fullfile(data_path, 'results');
% res_path = fullfile("/home/matteo/data/", 'results');
file_pattern = '*_run-all_v_2mni_norm_preproc0_despiked.nii.gz';

cort_atlas_path = fullfile(atlas_path, "HarvardOxford-cort-maxprob-thr25-2mm_SPLIT_RESAMPLED.nii.gz");
cort_labels = fullfile(atlas_path, 'HarvardOxford-cortical_Split.txt');

subcort_atlas_path = fullfile(atlas_path, "HarvardOxford-sub-maxprob-thr50-2mmRESAMPLED_REDUCED_FINAL.nii.gz");
subcort_labels = fullfile(atlas_path, 'HarvardOxford-Subcortical.txt');

HCP_atlas_path = fullfile(atlas_path, "MNI_Glasser_HCP_v1.0_resampled.nii.gz");
HCP_labels = fullfile(atlas_path, 'labels_HCP.txt');

feature_list = {'vgg_pca', ...
    'emo_pca', ...
    'mir_pca', ...
    'full_pca', ...
    %'mir_avg', ...
    %'essentia_avg', ...
    };
