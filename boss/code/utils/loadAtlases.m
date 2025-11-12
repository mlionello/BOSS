function [cort_atlas, subcort_atlas] = loadAtlases(atlas_path)
    cort_atlas_path = fullfile(atlas_path, "Schaefer2018_1000Parcels_17Networks_order_FSLMNI152_2mmRESAMPLED.nii.gz");
    subcort_atlas_path = fullfile(atlas_path, "HarvardOxford-sub-maxprob-thr50-2mmRESAMPLED_REDUCED_FINAL.nii.gz");

    cort_atlas = niftiread(cort_atlas_path);
    subcort_atlas = niftiread(subcort_atlas_path);
end