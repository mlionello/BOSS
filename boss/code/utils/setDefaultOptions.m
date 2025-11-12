function opts = setDefaultOptions(opts)
    % Set default options if not provided
    if ~isfield(opts, 'channel'), opts.channel = []; end
    if ~isfield(opts, 'nb_perm_step1'), opts.nb_perm_step1 = 0; end
    if ~isfield(opts, 'smooth_regions_padding'), opts.smooth_regions_padding = 0; end
    if ~isfield(opts, 'smooth_cortex_padding'), opts.smooth_cortex_padding = 0; end
    if ~isfield(opts, 'numPermutations_individual'), opts.numPermutations_individual = 0; end
    if ~isfield(opts, 'nb_perm'), opts.nb_perm = 0; end
    if ~isfield(opts, 'nb_signflip'), opts.nb_signflip = 10000; end
    if ~isfield(opts, 'nb_permAVG'), opts.nb_permAVG = 0; end
    if ~isfield(opts, 'extra'), opts.extra = 0; end
    if ~isfield(opts, 'winkler'), opts.winkler = 0; end
    if ~isfield(opts, 'checkpoint'), opts.checkpoint = ""; end
    if ~isfield(opts, 'suffix'), opts.suffix = ""; end
    if ~isfield(opts, 'cpus'), opts.cpus = 6; end
    if ~isfield(opts, 'target_model'), opts.target_model = "CanonCorr"; end
    if ~isfield(opts, 'radius'), opts.radius = [1 2 3 4]; end
    if ~isfield(opts, 'padding'), opts.padding = [1 1 2 3]; end
    if ~isfield(opts, 'pVxIn'), opts.pVxIn = 1; end
end
