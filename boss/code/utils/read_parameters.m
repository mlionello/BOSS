function [corr_thr, channel, smooth_regions_padding, smooth_cortex_padding, ...
    cort_pattern, subcort_pattern, all_cortex_analysis, extra, feat_all, ...
    mask, targetRegionMask, sourceRegionMask, corticalMask, subCorticalMask, ...
    cortAndSubMaks] = read_parameters(output_folder)

    % Read JSON from file
    fid = fopen(fullfile(output_folder, 'params.json'), 'r'); 
    json_str_params = fread(fid, '*char').'; 
    fclose(fid);
    fid = fopen(fullfile(output_folder, 'masks.json'), 'r'); 
    json_str_masks = fread(fid, '*char').'; 
    fclose(fid);

    % Decode JSON to struct
    params = jsondecode(json_str_params);
    masks = jsondecode(json_str_masks);

    % Extract variables from the struct
    corr_thr = params.corr_thr;
    channel = params.channel;
    smooth_regions_padding = params.smooth_regions_padding;
    smooth_cortex_padding = params.smooth_cortex_padding;
    cort_pattern = params.cort_pattern;
    subcort_pattern = params.subcort_pattern;
    all_cortex_analysis = params.all_cortex_analysis;
    extra = params.extra;
    feat_all = params.feat_all;

    mask = masks.mask;
    targetRegionMask = masks.targetRegionMask;
    sourceRegionMask = masks.sourceRegionMask;
    corticalMask = masks.corticalMask;
    subCorticalMask = masks.subCorticalMask;
    cortAndSubMaks = masks.cortAndSubMaks;
end
