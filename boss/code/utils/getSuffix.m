function suffix = getSuffix(in_suff)
    suffix = "";
    if ~isempty(in_suff)
        suffix = "_" + in_suff;
    end
end
