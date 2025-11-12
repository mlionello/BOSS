function [masks, filesubid] = createMasks(roi1_id, roi2_id, ...
        atlas_path, deriv_path, cort_atlas_path, subcort_atlas_path, ...
        hcp_atlas_path, smooth_regions, smooth_cortex, ...
        roi1_label, roi2_label)

    atlas_dir = dir(fullfile(atlas_path));
    nifti_dir = dir(fullfile(deriv_path));

    cort_atlas = single(niftiread(cort_atlas_path));
    subcort_atlas = single(niftiread(subcort_atlas_path));
    hcp_atlas = single(niftiread(hcp_atlas_path));

    % Initialize cell arrays for masks
    roi1_mask = cell(2, 1);
    roi2_mask = cell(1);

    % Generate masks for left and right cases for source
    for side = 1:2
        roi1_mask{side} = extract_mask(roi1_label{side}, cort_atlas, subcort_atlas, hcp_atlas, roi1_id{side}, smooth_regions);
    end

    % Generate masks for left and right cases for target
    if ~matches(roi2_label{1}, "allbrain")
        roi2_mask{1} = extract_mask(roi2_label{1}, cort_atlas, subcort_atlas, hcp_atlas, roi2_id{1}, smooth_regions);
    else
        roi2_mask{1} = niftiread(fullfile(atlas_path, "template_maskRESAMPLED.nii.gz"));
    end

    corticalMask = cort_atlas;
    subCorticalMask = subcort_atlas;
    corticalMask(corticalMask ~= 0) = 1;
    subCorticalMask(subCorticalMask ~= 0) = 1;

    cortAndSubMaks = corticalMask+subCorticalMask;
    cortAndSubMaks(cortAndSubMaks>1) = 1;
    if smooth_cortex
        cortAndSubMaks = smooth_mask(cortAndSubMaks, smooth_cortex);
    end

    filesubid = string({nifti_dir(3:end).name});

    masks = struct('roi2_mask', roi2_mask, ...
        'corticalMask', corticalMask, ...
        'subCorticalMask', subCorticalMask, ...
        'cortAndSubMaks', cortAndSubMaks);
    masks.roi1_mask = roi1_mask;
end

function mask = extract_mask(pattern, cort_atlas, subcort_atlas, hcp_atlas, roi_id, smooth_regions)
    mask =  zeros(size(cort_atlas));
    for i = 1: length(pattern)
        split_pattern = strsplit(pattern{i}, ':');
        if strcmp(split_pattern{1}, 'cort')
            mask = mask | ismember(cort_atlas, single(roi_id(i)));
        elseif strcmp(split_pattern{1}, 'sub')
            mask = mask | ismember(subcort_atlas, single(roi_id(i)));
        elseif strcmp(split_pattern{1}, 'hcp')
            mask = mask | ismember(hcp_atlas, single(roi_id(i)));
        end
    end
    if smooth_regions
        mask = smooth_mask(mask, smooth_regions);
    end
end
