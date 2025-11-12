clear;
maxNumCompThreads(28)
% import constants defined in ../constants.m
addpath ../; constants;
% import music regressors from ./musicregrExtract.m
musicregrExtract;
% import support functions;
addpath ./utils/
% winkler correction for cca not relevant when working only on the first component
winkler=false;
fprAnalysis = 0;
corr_thr=0.5;

atlas_dir = dir(fullfile(atlas_path));
nifti_dir = dir(fullfile(deriv_path));

cort_atlas_path = fullfile(atlas_path, "Schaefer2018_1000Parcels_17Networks_order_FSLMNI152_2mmRESAMPLED.nii.gz");
subcort_atlas_path = fullfile(atlas_path, "HarvardOxford-sub-maxprob-thr50-2mmRESAMPLED_REDUCED_FINAL.nii.gz");

cort_atlas = niftiread(cort_atlas_path);
subcort_atlas = niftiread(subcort_atlas_path);

labels = arrayfun(@(x) replace(x, '_', "\_"), labels);
[maxvalue , rm_feats] = max(sum((corrcoef(musicRegr) - eye(size(musicRegr,2))) > corr_thr ,1));
fprintf("Removing columns: ");
while maxvalue > 0
    fprintf("%s; ", labels(rm_feats));
    musicRegr(:, rm_feats) = [];
    labels(rm_feats) = [];
    [maxvalue , rm_feats] = max(sum((corrcoef(musicRegr) - eye(size(musicRegr,2))) > corr_thr ,1));
end
fprintf("\n");


%% DATASET LOADING
cort_labels = fullfile(atlas_path, 'Schaefer2018_1000Parcels_17Networks_order.txt');
subcort_labels = fullfile(atlas_path, 'HarvardOxford-Subcortical.txt');
cort_pattern = "RH_SomMotB_Aud"; % RH_TempPar
subcort_pattern = "rightamyg";
[cort_roi_id, cort_roi_label] = extractROIIds(cort_labels, cort_pattern);
[subcort_roi_id, subcort_roi_label] = extractROIIds(subcort_labels, subcort_pattern);

% Display the result
disp('ROI IDs corresponding to lines with "' + cort_pattern + '":');
disp(cort_roi_id');
% Display the result
disp('ROI IDs corresponding to lines with "' + subcort_pattern + '":');
disp(subcort_roi_id');


%% THEMASK - sfumeggiante
targetRegionMask = ismember(cort_atlas, cort_roi_id);
sourceRegionMask = ismember(subcort_atlas, subcort_roi_id);
corticalMask = cort_atlas;
corticalMask(corticalMask ~= 0) = 1;
filesubid = string({nifti_dir(3:end).name});
fprRepRes = cell(1,fprAnalysis);
fprRepRDiff = cell(1,fprAnalysis);

for fprRep = 1: max(1, fprAnalysis)
    t_init = tic;
    fprintf("repetition %d out of %d\n", fprRep, max(1, fprAnalysis))
    tsum = []; ssum = [];
    subj = {};
    for sub_i = 1:length(filesubid)
        % loading participant nifti image
        nifti_image = niftiread(fullfile(nifti_dir(3).folder, filesubid(sub_i), 'func', filesubid(sub_i)+"_run-all_preproc_despkd.nii.gz"));
        info = niftiinfo(fullfile(nifti_dir(3).folder, filesubid(sub_i)));
        if ~fprAnalysis
            % targetRegionMask_indices = find(targetRegionMask);
            % target_region_new = nan(length(targetRegionMask_indices), size(nifti_image, 4));
            % for indx = 1: length(targetRegionMask_indices)
            %     target_region_new(indx, :) = nifti_image( ...
            %         targetRegionMask_indices(indx) : ...
            %         numel(targetRegionMask) : ...
            %         end);
            % end
            targetRegion = nifti_image .* targetRegionMask;
        else
            % if analysing fpr for this method
            targetRegion = zeros(size(nifti_image));
            corticalMaskIndices = find(corticalMask);
            [xinit, yinit, zinit] = ind2sub(size(corticalMask), corticalMaskIndices);

            permMaskInd = randperm(numel(corticalMaskIndices));
            permutedCorticalMaskIndices = corticalMaskIndices(permMaskInd);
            [x,y,z] = ind2sub(size(corticalMask), permutedCorticalMaskIndices);
            for ind = 1 : length(corticalMaskIndices)
                intargetarea = targetRegionMask(xinit(ind), yinit(ind), zinit(ind));
                ind2 = ind;
                while intargetarea && ~any(nifti_image(x(ind2), y(ind2), z(ind2), :), 'all')
                    ind2 = randi(length(corticalMaskIndices));
                end
                targetRegion(xinit(ind), yinit(ind), zinit(ind), :) = nifti_image(x(ind2), y(ind2), z(ind2), :) .* intargetarea;
            end

        end
        % sourceRegionMask_indices = find(sourceRegionMask);
        % source_region_new = nan(length(sourceRegionMask_indices), size(nifti_image, 4));
        % for indx = 1: length(sourceRegionMask_indices)
        %     source_region_new(indx, :) = nifti_image( ...
        %         sourceRegionMask_indices(indx) : ...
        %         numel(sourceRegionMask) : ...
        %         end);
        % end
        sourceRegion = nifti_image .* sourceRegionMask;
        tsum(sub_i) = sum(targetRegion~=0, 'all');
        ssum(sub_i) = sum(sourceRegion~=0, 'all');
        

        %% SEARCHLIGHT 1D TEMPLATE
        targetRadius = 1; targetPadding = 2;
        sourceRadius = 1; sourcePadding = 2;
        sourcePropInVx = 0.6;
        targetPropInVx = 0.8;

        targetRegion2D = reshape(targetRegion, [], size(targetRegion, 4));
        sourceRegion2D = reshape(sourceRegion, [], size(sourceRegion, 4));
        targetRegionMask1D = reshape(targetRegionMask, [], 1);
        sourceRegionMask1D = reshape(sourceRegionMask, [], 1);

        targetRoiTemplate = create_template(size(targetRegionMask), ...
            targetRadius, targetPadding);
        sourceRoiTemplate = create_template(size(sourceRegionMask), ...
            sourceRadius, sourcePadding);


        %% SEARCHLIGHT
        permSchema=[];
        permSchema(:, 1) = 1:size(musicRegr{1}, 1);
        for i = 1: 0
            permSchema(:, i+1) = randperm(size(musicRegr{1}, 1));
        end

        mdata = struct('musicRegr', musicRegr, 'labels', labels, ...
            'permSchema', permSchema, 'winkler', winkler);

        tic
        results = searchlight(sourceRegionMask1D, sourceRegion2D, ...
            sourceRoiTemplate, sourcePropInVx, targetPropInVx, ...
            targetRegionMask1D, targetRegionMask, targetRegion2D, ...
            targetRoiTemplate, perform_cca_on_searchlight, mdata);
        toc


        %% cell2StructMat
        results = resCell2Mat(cca_values, winkler);
        subj{sub_i}.rSourceVsTarget = results.rSourceVsTarget;
        subj{sub_i}.rSourceVsTarget = results.resCorrPerm;
        subj{sub_i}.rDiff = results.rSourceVsTarget - results.resCorrPerm;
    end
    %%
    % rDiff field non-empty structs
    rDiffCells = cellfun(@(x) reshape(x.rDiff, 1, size(x.rDiff,1), size(x.rDiff,2)), subj(~cellfun('isempty', subj)), 'UniformOutput', false);
    resCorrPerm = cellfun(@(x) reshape(x.rSourceVsTarget, 1, size(x.rSourceVsTarget,1), size(x.rSourceVsTarget,2)), subj(~cellfun('isempty', subj)), 'UniformOutput', false);
    allRDiff = cat(1, rDiffCells{:});
    allresCorrPerm = cat(1, resCorrPerm{:});
    allrSourceVsTarget = allRDiff + allresCorrPerm;

    allresCorrPermReshaped =  reshape(allresCorrPerm, size(allresCorrPerm, 1), []);
    allrSourceVsTargetReshaped =  reshape(allrSourceVsTarget, size(allrSourceVsTarget, 1), []);

    axisToTest = 1;
    numPermutations = min(10, 2^size(allRDiff, axisToTest));

    signFlipSchema = sign(randn(numPermutations, size(allRDiff, axisToTest)));
    % signFlipSchema = reshape(signFlipSchema, size(signFlipSchema, 1), 1, 1, size(signFlipSchema, 2));
    % signFlipSchema = repmat(signFlipSchema, [1, size(allRDiff, 1), size(allRDiff, 2), 1]);

    cr = "";
    nulldist = zeros(numPermutations, size(allrSourceVsTargetReshaped, 2));
    for perm = 1: numPermutations
        msg = sprintf("computing permutation %d out of %d", perm , numPermutations);
        fprintf(cr + msg);
        cr = repmat('\b', 1, strlength(msg));
        % nulldist(perm, :, :) = mean( signFlipSchema(perm, :)' .* allRDiff, 1);
        for pair = 1: length(allresCorrPermReshaped)
            x = cat(2, allrSourceVsTargetReshaped(:, pair), allresCorrPermReshaped(:, pair));
            permute = signFlipSchema(perm, :)>0;
            x(permute, :) = x(permute, [2, 1]);
            [~, ~, ~, stats] = ttest(x(:, 1), x(:, 2), 'Tail','right');
            nulldist(perm, pair) = stats.tstat;
        end
    end
    fprintf("...  done\n")

    % for pair = 1: length(allresCorrPermReshaped)
    %     x = cat(2, allrSourceVsTargetReshaped(:, pair), allresCorrPermReshaped(:, pair));
    %     [~, ~, ~, stats] = ttest(x(:, 1), x(:, 2));
    %     tscore(1, pair) = stats.tstat;
    % end



    fprintf("computing tiedrank")
    allRDiffreshape = reshape(allRDiff, size(allRDiff, 1), []);
    concat = cat(1, mean(allRDiff, 1), nulldist);
    maxstat = repmat(max(nulldist, [], [2, 3]), 1, size(nulldist, 2), size(nulldist, 3));
    pval = tiedrank(-concat)/size(concat, 1);
    pval_corr = tiedrank(-cat(1, mean(allRDiff, 1), maxstat))/size(maxstat, 1);
    fprintf("...  %d:%d done\n",  floor(toc(t_init)/60), rem(toc(t_init),60))
    fprRepRes{fprRep} = squeeze(pval(1, :, :));
    fprRepRDiff{fprRep} = allRDiff;
end


%% Anlaysis
for j = 1:size(cca_values.r, 1)
    for i = 1:size(cca_values.r, 2)
        r(j, i) = cca_values.r{j, i}(1);
        p(j, i) = cca_values.p{j, i}(1);
        roi1_coord(j, i, :) = cca_values.roi_center{j, i};
    end
    if ~winkler
        for i = 1:size(cca_values.r, 2)
            wilks(j,i) = cca_values.wilks{j, i}(1);
            pF(j,i) = cca_values.pF{j, i}(1);
            F(j,i) = cca_values.F{j, i}(1);
        end
    end
    coord(j, :) = cca_values.roi2_center{j};
end
% [sign_subj, sign_subi] = ind2sub(size(cca_values.r), find((p/size(cca_values.r, 2))<0.001));
% for n = 1: length(sign_subj)
%     r_significant(n) = r(sign_subj(n), sign_subi(n));
% end
%
% for j = 1:size(cca_values.r, 1)
%     sign_ind = find((p(j, :)/size(cca_values.r, 2))<0.001);
%     [arg, roi2_back(i)] = max(r(j, sign_ind));
% end

for i = 1:size(cca_values.r, 2)
    [r_winning(i), arg] = max(r(:, i));
    center_winning(i, :) = coord(arg, :); % amygdala coordiantes of winning roi
end

sourceSphere = 1;
figure;
subplot(311); plot(squeeze(results.rSourceVsTarget(sourceSphere, :)));
xlabel('spheres along target region')
ylabel(compose('r with %dth source sphere', sourceSphere))
title('residual correlations')
subplot(312); plot(squeeze(results.resCorrPerm(sourceSphere, :, 1, :)));
xlabel('spheres along target region')
ylabel(compose('r with %dth source sphere', sourceSphere))
title('full correlations')
subplot(313); plot((squeeze(results.rSourceVsTarget(sourceSphere, :))-squeeze(results.resCorrPerm(sourceSphere, :, 1, 1))));
xlabel('spheres along target region')
ylabel(compose('r with %dth source sphere', sourceSphere))
title('abs diff')

%% analysis on residual
resCorrPerm = permute(results.resCorrPerm, [3, 1, 2, 4]);

% pval calculation:
ranked = tiedrank(-resCorrPerm);
resCorrPVal = squeeze(ranked(1, :, :, :)/size(ranked, 1));
resSignificant = resCorrPVal<0.05;
MusFeatInd = sum(resSignificant, [1, 2]);
[MusFeatIndSorted, sortInd] = sort(MusFeatInd, 'descend');
for i = 1:length(sortInd)
    fprintf('%-12s --> %d signifcant datapoints; \n', ...
        labels(sortInd(i)), MusFeatIndSorted(i));
end
resCorrPermMostSignFeat = squeeze(resCorrPerm(1,:,:, sortInd(1)));
resCorrSignificant = squeeze(resCorrPerm(1, :, :, sortInd(1)));
resCorrSignificant(~resSignificant(:, :, sortInd(1))) = 0;
winnerMaxTarget = max(resCorrSignificant, [], 1);
niftiSignCorr = zeros([sourceRoiTemplate.size, 2]);
niftiSignCorr=drawnifti(niftiSignCorr,results.targetRoiCoord,targetRoiTemplate,max(results.rSourceVsTarget,[],1),1);
niftiSignCorr=drawnifti(niftiSignCorr,results.targetRoiCoord,targetRoiTemplate,winnerMaxTarget,2);
info_outniftiAllTargetRois = info;
info_outniftiAllTargetRois.ImageSize = size(niftiSignCorr);
niftiwrite(single(niftiSignCorr), 'niftiSignCorr_sub1', info_outniftiAllTargetRois, 'Compressed', true);

% corrMax calculation:
resCorrPermMax = max(resCorrPerm(2: end, :, :, :), [], [2, 3]);
resCorrPermMax = repmat(resCorrPermMax, [1, size(resCorrPerm, 2), size(resCorrPerm, 3), 1]);
ranked  = tiedrank(-cat(1, resCorrPerm(1,:,:,:), resCorrPermMax));
fw_significant = squeeze(ranked(1, :, :, :)/size(ranked,1)<0.05);
FW_MusFeatInd = sum(fw_significant, [1, 2]);
[FW_MusFeatIndSorted, FW_sortInd] = sort(FW_MusFeatInd, 'descend');
for i = 1:length(FW_sortInd)
    fprintf('%-12s --> %d signifcant datapoints; \n', ...
        labels(FW_sortInd(i)), FW_MusFeatIndSorted(i));
end

%% Analysis: Eigendecomposition over Laplacian
minDist = 404;
[mappedPoints, mappingMatrix] = analyse_connectivity(r, minDist);
resOnRMapped = rres'*mappingMatrix;
for perm_i = 1:size(resCorrperm, 3)
    resOnRMappedPERM{perm_i} = resCorrperm(:,:,perm_i)'*mappingMatrix;
    distRvsResPerm(perm_i) = norm(resOnRMappedPERM{perm_i}(:, 1) - mappedPoints(:, 1));
end

[mappedPointsRes, mappingMatrixRes] = analyse_connectivity(rres, minDist);
rOnResMapped = r'*mappingMatrixRes;
distInRSpace = norm(resOnRMapped(:, 1) - mappedPoints(:, 1));
pnormInRSpace = tiedrank(-[ distInRSpace, distRvsResPerm]);
pnormInRSpace = pnormInRSpace(1)/(numel(distRvsResPerm)+1);
anovaResultInRSpace = anova1(cat(2, resOnRMapped(:, 1), mappedPoints(:, 1)), [], 'off');
distInResSpace = norm(resOnRMapped(:, 2) - mappedPoints(:, 2));
pnormInResSpace = tiedrank(-[distInResSpace, distRvsResPerm]);
pnormInResSpace = pnormInResSpace(1)/(numel(distRvsResPerm)+1);
anovaResultInResSpace = anova1(cat(2, rOnResMapped(:, 1), mappedPointsRes(:, 1)), [], 'off');
figure;
subplot(211)
scatter(mappedPoints(:, 1), mappedPoints(:, 2), 'o'); hold on;
scatter(resOnRMapped(:, 1), resOnRMapped(:, 2), 'xr');
xlabel('Mapped Dimension 1');
ylabel('Mapped Dimension 2');
title(compose('R mapped space: anova: %.3f; mean norm: %.2f (p=%.2f)', anovaResultInRSpace, distInRSpace, pnormInRSpace))
subplot(212)
scatter(mappedPointsRes(:, 1), mappedPointsRes(:, 2), 'o'); hold on;
scatter(rOnResMapped(:, 1), rOnResMapped(:, 2), 'xr');
xlabel('Mapped Dimension 1');
ylabel('Mapped Dimension 2');
title(compose('Residual mapped space: anova: %.3f; mean norm: %.2f (p=%.2f)', anovaResultInResSpace, distInResSpace, pnormInResSpace))

ConnectivityMatrix = abs(mappingMatrix * mappingMatrix');
ConnectivityMatrixRes = abs(mappingMatrixRes * mappingMatrixRes');

subplot(1, 2, 1);
imagesc(log(ConnectivityMatrix));
colorbar;
title('Connectivity Strength - r');
xlabel('Voxel Index');
ylabel('Voxel Index');
subplot(1, 2, 2);
imagesc(log(ConnectivityMatrixRes));
colorbar;
title('Connectivity Strength - r\_res');
xlabel('Voxel Index');
ylabel('Voxel Index');

info_outniftiAllTargetRois = info;

outniftiMappedPointsSource = zeros([sourceRoiTemplate.size, 2]);
outniftiMappedPointsSource = drawnifti(outniftiMappedPointsSource, targetRoiCoord, sourceRoiTemplate, mappedPoints(:, 1), 1);
outniftiMappedPointsSource = drawnifti(outniftiMappedPointsSource, targetRoiCoord, sourceRoiTemplate, mappedPoints(:, 2), 2);

outniftiMappedPointsSourceRes = zeros([sourceRoiTemplate.size, 2]);
outniftiMappedPointsSourceRes = drawnifti(outniftiMappedPointsSourceRes, targetRoiCoord, sourceRoiTemplate, mappedPointsRes(:, 1), 1);
outniftiMappedPointsSourceRes = drawnifti(outniftiMappedPointsSourceRes, targetRoiCoord, sourceRoiTemplate, mappedPointsRes(:, 2), 2);

maxValue = max(max(outniftiMappedPointsSource, [], [1, 2, 3]), max(outniftiMappedPointsSourceRes, [], [1, 2, 3]));
minValue = min(min(outniftiMappedPointsSource, [], [1, 2, 3]), min(outniftiMappedPointsSourceRes, [], [1, 2, 3]));

outniftiMappedPointsSource = (outniftiMappedPointsSource - minValue) ./ (maxValue - minValue);
outniftiMappedPointsSourceRes = (outniftiMappedPointsSourceRes - minValue) ./ (maxValue - minValue);

bckgrnd = find(~targetRegionMask);
for i = 1:length(bckgrnd)
    outniftiMappedPointsSource(bckgrnd(i):numel(targetRegionMask):end)=0;
    outniftiMappedPointsSourceRes(bckgrnd(i):numel(targetRegionMask):end)=0;
end
info_outniftiAllTargetRois.ImageSize = size(outniftiMappedPointsSourceRes);
niftiwrite(single(outniftiMappedPointsSource), 'outniftiMappedPointsTarget_sub1', info_outniftiAllTargetRois, 'Compressed',true)
niftiwrite(single(outniftiMappedPointsSourceRes), 'outniftiMappedPointsTargetRes_sub1', info_outniftiAllTargetRois, 'Compressed',true)

%%
significance = (p/numel(cca_values.r))<0.05;

% all components - come together right now over me
outniftiAllTargetRois = zeros([sourceRoiTemplate.size, size(r,1)]);
for sourceRoi = 1:size(r, 1)
    outniftiAllTargetRois = drawnifti(outniftiAllTargetRois, coord(sourceRoi), sourceRoiTemplate, 1, sourceRoi);
    outniftiAllTargetRois = drawnifti(outniftiAllTargetRois, targetRoiCoord, targetRoiTemplate, r(sourceRoi, :), sourceRoi);
end
info_outniftiAllTargetRois = info;
info_outniftiAllTargetRois.ImageSize = size(outniftiAllTargetRois);
niftiwrite(single(outniftiAllTargetRois), 'outniftiAllTargetRois_sub1', info_outniftiAllTargetRois, 'Compressed',true)

outniftiAllTargetRoisRES = zeros([sourceRoiTemplate.size, size(rres,1)]);
for sourceRoi = 1:size(rres, 1)
    outniftiAllTargetRoisRES = drawnifti(outniftiAllTargetRoisRES, coord(sourceRoi), sourceRoiTemplate, 1, sourceRoi);
    outniftiAllTargetRoisRES = drawnifti(outniftiAllTargetRoisRES, targetRoiCoord, targetRoiTemplate, rres(sourceRoi, :), sourceRoi);
end
niftiwrite(single(outniftiAllTargetRoisRES), 'outniftiAllTargetRoisRES_sub1', info_outniftiAllTargetRois, 'Compressed',true)

outniftiMusSrcTarget = zeros([sourceRoiTemplate.size, 2]);
outniftiMusSrcTarget = drawnifti(outniftiMusSrcTarget, coord, sourceRoiTemplate, rMusSource, 1);
outniftiMusSrcTarget = drawnifti(outniftiMusSrcTarget, targetRoiCoord, targetRoiTemplate, rMusTarget, 2);

info_outniftiAllTargetRois.ImageSize = size(outniftiMusSrcTarget);
niftiwrite(single(outniftiMusSrcTarget), 'outniftiMusSrcTarget_sub1', info_outniftiAllTargetRois, 'Compressed',true)

% the winner takes it all - no time for losers 'cause we are the champions
[r_winning , arg_r] = max(r, [], 1);
sourceRoiWinning = coord(arg_r);
uniqueWinningRoi = unique(sourceRoiWinning);

outnifti = zeros([sourceRoiTemplate.size, numel(uniqueWinningRoi)]);
for winnignRoi_i = 1: numel(uniqueWinningRoi)
    sourceWinningRoiIndices = uniqueWinningRoi(winnignRoi_i);
    winnignRoi_i_indx = find(sourceRoiWinning == sourceWinningRoiIndices);

    numVoxSecodnary(winnignRoi_i) = numel(winnignRoi_i_indx);

    outnifti = drawnifti(outnifti, sourceWinningRoiIndices, sourceRoiTemplate, 1, winnignRoi_i);
    %outnifti = drawnifti(outnifti, targetRoiCoord(winnignRoi_i_indx), targetRoiTemplate, r_winning(winnignRoi_i_indx), winnignRoi_i);
    outnifti = drawnifti(outnifti, targetRoiCoord(winnignRoi_i_indx), targetRoiTemplate, r_winning(winnignRoi_i_indx), winnignRoi_i);
end
info_outnifti = info;
info_outnifti.ImageSize = size(outnifti);
niftiwrite(single(outnifti), 'outnifti_sub1', info_outnifti, 'Compressed', true)


significant = find(all(significance, 1));
significantR = r_winning(significant)';
[CWS_source(:, 1), CWS_source(:, 2), CWS_source(:, 3)] = ind2sub(targetRoiTemplate.size, sourceRoiWinning(significant));
[CS_target(:, 1), CS_target(:, 2), CS_target(:, 3)] = ind2sub(sourceRoiTemplate.size, targetRoiCoord(significant));

figure
sgtitle('only significant rois')
for dim_n = 1: 3
    fprintf('\ndim %d: ', dim_n)
    subplot(3, 1, dim_n)

    [sortedCols, sortign_ind] = sortrows([significantR, CWS_source(:, dim_n), CS_target(:, dim_n)], [2,3]);
    significantRSorted(dim_n, :) = sortedCols(:, 1);
    CWS_source_sorted(dim_n, :) = sortedCols(:, 2);
    CS_target_sorted(dim_n, :) = sortedCols(:, 3);

    % signSorted = significant(sortign_ind);
    % CWS_source_sorted = sourceRoiWinning(signSorted, dim_n);
    % [CS_target_sorted(:, 1), CS_target_sorted(:, 2), CS_target_sorted(:, 3)] = ind2sub(sourceRoiTemplate.size, targetRoiCoord(signSorted));
    % significantRSorted(dim_n, :) = r_winning(signSorted);

    plot(squeeze(significantRSorted(dim_n, :)));
    hold on;

    roiChange = find(diff(CWS_source_sorted(dim_n, :)));

    xline(roiChange);
    xticks(1:6:length(CS_target_sorted(dim_n,:)))
    xticklabels(CS_target_sorted(dim_n, 1:6:length(CS_target_sorted(dim_n,:))))

    roiInd_range =  [1, roiChange, length(CWS_source_sorted(dim_n, :))];
    clear rRange rRangeStat
    for roiInd = 1:length(roiInd_range)-1
        rRange{roiInd} = significantRSorted(roiInd_range(roiInd):roiInd_range(roiInd+1));
        rRangeStat(roiInd).avg = mean(rRange{roiInd});
        rRangeStat(roiInd).stdev = std(rRange{roiInd});
        fprintf('%.2f%s%.2f; ', rRangeStat(roiInd).avg, char(177), rRangeStat(roiInd).stdev);
    end
    fprintf('\n')

    p_values = zeros(length(rRange), length(rRange));
    for i = 1:length(rRange)
        for j = 1:length(rRange)
            if i ~= j
                [~, p_values(i, j)] = kstest2(rRange{i}, rRange{j});
            end
        end
    end
    disp(p_values)
end

figure
sgtitle('all rois, x for non significant')
for dim_n = 1: 3
    subplot(3, 1, dim_n)
    [sortedWC, sortign_ind] = sort(CWS_source(:, dim_n));
    nonSignificant = find(any(~significance(:, sortign_ind), 1));
    rSorted(dim_n, :) = r_winning(sortign_ind);
    plot(squeeze(rSorted(dim_n, :)));
    hold on;
    roiChange = find(diff(sortedWC(significant)))+1;
    roiInd_range = [1, roiChange', length(significantRSorted)];
    xline(roiChange);
    for i = 1: length(nonSignificant)
        scatter(nonSignificant(i), squeeze(rSorted(dim_n, nonSignificant(i))), 'xr')
    end
end


%% SUPPORT FUNCTIONS
function cort_roi_id = extractROIIds(filePath, pattern)
cort_roi_id = [];

fid = fopen(filePath, 'r');
if fid == -1
    error('Error opening the file.');
end

while ~feof(fid)
    currentLine = fgetl(fid);
    roi_id = textscan(currentLine, '%f');

    lastString = split(currentLine, ' ');
    if contains(lastString{end}, pattern)
        cort_roi_id = [cort_roi_id; roi_id{1}(1)];
    end
end
fclose(fid);
end

function values_out = get_roi(region, center, radius)
% Extracts a spherical region of interest (ROI) from a 3D binary region
% centered at the specified coordinates with the given radius.

% Ensure the input coordinates are valid
if ~isnumeric(center) || numel(center) ~= 3
    error('Invalid center coordinates. Must be a 3-element numeric array.');
end

% Ensure the input region is a 4D array
if ndims(region) ~= 4
    error('Invalid region. Must be a 3D binary array.');
end

% Get the size of the region
[cols, rows, slices, t_points] = size(region);

% Ensure the center coordinates are within the valid range
if any(center < 1) || any(center > [rows, cols, slices])
    error('Center coordinates are outside the valid range.');
end

% Create a grid of coordinates
[X, Y, Z] = ndgrid(1:cols, 1:rows, 1:slices);

% Compute the distance from each point to the center
distances = sqrt((X - center(1)).^2 + (Y - center(2)).^2 + (Z - center(3)).^2);

% Create a binary mask for the spherical ROI
roi_mask = distances <= radius;

% Extract the values within the ROI
nonzeros_indices = find(roi_mask);
numel_mask = numel(roi_mask);
values_out = zeros(size(region, 4), numel(nonzeros_indices));
for n = 1:numel(nonzeros_indices)
    values_out(:, n) = squeeze(region(nonzeros_indices(n):numel_mask:end));
end
end

function updatedCenters = update_availableCenters(availableCenters, currentCenter, padding)
% Get the size of the availableCenters matrix
[rows, cols, slices] = size(availableCenters);

% Create a grid of coordinates
[X, Y, Z] = ndgrid(1:rows, 1:cols, 1:slices);

% Compute the distance from each point to the center
distances = sqrt((X - currentCenter(1)).^2 + (Y - currentCenter(2)).^2 + (Z - currentCenter(3)).^2);

% Define a spherical mask for the ROI
roi_mask = distances < padding;

% Update the availability status of neighboring centers
updatedCenters = availableCenters & ~roi_mask;
end
