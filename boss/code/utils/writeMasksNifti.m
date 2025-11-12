function writeMasksNifti(output_folder, masks, info_nifti)
    mask_nifti_path = fullfile(output_folder, 'masks.nii.gz');
    if ~exist(mask_nifti_path, 'file')
        cat_mask = cat(4, masks.roi1_mask{:}, masks.roi2_mask);
        out_info = info_nifti;
        out_info.ImageSize = size(cat_mask);
        niftiwrite(cat_mask, mask_nifti_path, out_info, 'Compressed', true);
        fprintf("Saved masks as NIfTI to %s\n", mask_nifti_path);
    else
        fprintf("Masks NIfTI already exists at %s\n", mask_nifti_path);
    end
end
