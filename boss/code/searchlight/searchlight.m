function [results, center_indices] = searchlight(sourceRoi, ...
    targetPropInVx, targetRegionMask, targetRegion2D, ...
    targetRoiTemplate, analyser, mdata)
    
    availablecentres = reshape(targetRegionMask, [], 1);
    init_availablereg = sum(availablecentres, 'all');
    cr = "";
    init_available_centers_list = find(targetRegionMask);
    results = cell(numel(init_available_centers_list), 1);
    center_indices = zeros(numel(init_available_centers_list), 1);
    tic;
    last_toc = toc;
    
    for cursor = 1: numel(init_available_centers_list)
        targetRoiCenter = init_available_centers_list(cursor);
        if ~availablecentres(targetRoiCenter) 
            continue
        end
        targetRoi = extract_roi_from_template(targetRegion2D, targetRoiCenter, targetRoiTemplate);
        if size(targetRoi, 2) < int32(targetPropInVx * targetRoiTemplate.maxnbvx)
            continue
        end

        port_undone = sum(availablecentres, 'all') / init_availablereg;
        port_done = 1 - port_undone;

        if toc > last_toc + 1
            msg = compose("\tprocessing targetRoi: %.1f/100 eta: %s elapsed: %s", ...
                100*port_done, string(duration(0, 0, toc * port_undone / port_done)), ...
                string(duration(0, 0, toc)));
            fprintf(cr + msg);
            cr = repmat('\b', 1, strlength(msg));
            last_toc = toc;
        end

        % compute your analysis
        if ~isempty(analyser)
            results{cursor} = analyser(sourceRoi, targetRoi, mdata); %TODO: support for svm!!
        end
        center_indices(cursor) = targetRoiCenter;

        % update the volumes of the candidate coordinates for the
        % center of targetRoi by considering the given padding
        availablecentres = update_availableCenters(availablecentres, ...
            targetRoiCenter, targetRoiTemplate); %this is why i cannot move to parfor

    end
    fprintf('\n')
    keep_indices = find(center_indices); % center_indices is initialized with zeros
    center_indices = center_indices(keep_indices);
    if ~isempty(analyser)
        results = results(keep_indices);
    end

end
