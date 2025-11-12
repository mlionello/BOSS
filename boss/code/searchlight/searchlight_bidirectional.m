function results = searchlight_bidirectional(sourceMask1D, source2D, sourceTemplate, sourcePropInVx, targetPropInVx, targetMask1D, targetMask, target2D, targetTemplate,  analyser, mdata)
    availablecentres2 = sourceMask1D;
    init_availablereg2 = sum(availablecentres2, 'all');
    jj = 1;
    cr = "";
    for sourceCenter = 1: numel(sourceMask1D)
        if ~availablecentres2(sourceCenter)
            continue
        end
        source = extract_roi_from_template(source2D, sourceCenter, sourceTemplate);
        if size(source, 2) < int32(sourcePropInVx*sourceTemplate.maxnbvx)
            continue
        end
          
        results{jj} = searchlight(source, targetPropInVx, targetMask1D, targetMask, target2D, targetTemplate, analyser, mdata);
        results.sourceCenter{jj} = sourceCenter;
        availablecentres2 = update_availableCenters(availablecentres2, sourceCenter, sourceTemplate);
        jj = jj + 1;
    end
    fprintf(' Done\n')
end
