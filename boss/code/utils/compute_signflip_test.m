function [pval, maxstat] = compute_signflip_test(resCorr, raw_corr, numPermutations)

    signFlipSchema = sign(randn(numPermutations, size(resCorr, 1)));
    signFlipSchema = signFlipSchema';

    nulldist = cell(1, numPermutations);
    ttest_data = zeros(1, size(resCorr, 2));

    ttest_data(1, :) = my_ttest(raw_corr - resCorr);
    cr = "";
    maxstat = zeros(numPermutations, 1);
    for perm = 1: numPermutations
        msg = compose("%d out of %d", perm, numPermutations);
        fprintf(cr + msg);
        cr = repmat('\b', 1, strlength(msg));
        nulldist = my_ttest((raw_corr - resCorr).*signFlipSchema(:,perm));
        maxstat(perm) = max(nulldist);
    end
    fprintf(cr); cr="";

    %maxstat = cellfun(@(x) max(x), nulldist);

    pval.fivecent = ttest_data > prctile(maxstat, 95);
    pval.cent = ttest_data > prctile(maxstat, 99);
    pval.mill = ttest_data > prctile(maxstat, 99.9);
    pval.fivecent_twotails = ttest_data > prctile(maxstat, 97.5);
    pval.cent_twotails = ttest_data > prctile(maxstat, 99.5);
    pval.mill_twotails = ttest_data > prctile(maxstat, 99.95);
%     pval = nan(1, numel(ttest_data));
%     for i = 1:numel(ttest_data)
%         msg = compose("pval: %d out of %d", i, numel(ttest_data));
%         fprintf(cr + msg);
%         cr = repmat('\b', 1, strlength(msg));
%         pva_res = tiedrank(-cat(1, ttest_data(i), maxstat))/(size(maxstat, 1)+1);
%         pval(1, i) = pva_res(1);
%     end

end


function tstat = my_ttest(x)
    tstat = mean(x) ./ (std(x) ./ sqrt(size(x,1)));
end
