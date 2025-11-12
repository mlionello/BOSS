function template = create_template(volumeShape, radius, padding)
    [X, Y, Z] = ndgrid(1:volumeShape(1), 1:volumeShape(2), 1:volumeShape(3));
    tmpl_coord = [floor(max(X, [], 'all')/2), floor(max(Y, [], 'all')/2), floor(max(Z, [], 'all')/2)];
    % Compute the distance from each point to the center
    distances = sqrt((X - tmpl_coord(1)).^2 + (Y - tmpl_coord(2)).^2 + (Z - tmpl_coord(3)).^2);
    % Create a binary mask for the spherical ROI
    template_3D = distances <= radius;
    % get the set list of indices from the template for a roi centered at 0
    template.indices = find(template_3D)-sub2ind(size(template_3D), ...
        tmpl_coord(1), ...
        tmpl_coord(2), ...
        tmpl_coord(3));
    template.size = size(template_3D);
    template.radius = radius;
    template.padding = padding;
    template.maxnbvx = sum(template_3D, 'all');
end
