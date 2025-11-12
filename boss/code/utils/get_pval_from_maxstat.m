function pval = get_pval_from_maxstat(corr)
    if size(corr, 1)>1
        
        nulldist = max(corr(2:end, :), [], 2);
        pval = tiedrank( -cat(1, corr(1, :), repmat(nulldist, 1, size(corr,2))))/size(corr,1);
        pval = pval(1,:);
    else
        pval = nan;
    end
end