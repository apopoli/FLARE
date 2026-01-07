function BCval_p = set_BCs_on_nodes(msh, BC)
    % Number of mesh nodes
    np = size(msh.POS, 1);

    % Boundary edges and tags
    edges = msh.LINES(:, 1:2);   % node indices
    tags  = msh.LINES(:, 3);     % physical tags

    % Initialize BC table for ALL nodes
    % [node_id, BC_type, BC_value]
    % BC_type: 1 = Neumann, 2 = Dirichlet
    BCval_p = zeros(np, 3);
    BCval_p(:,1) = (1:np).';     % global node numbering

    % Boundary nodes only
    boundary_nodes = unique(edges(:));

    % Loop on boundary nodes
    for k = 1:length(boundary_nodes)
        node = boundary_nodes(k);

        % All boundary edges touching this node
        rows = (edges(:,1) == node) | (edges(:,2) == node);
        node_tags = tags(rows);

        % ---- Dirichlet has priority ----
        idir = ismember(node_tags, BC.D.tag);
        if any(idir)
            tag = node_tags(find(idir,1));
            BCval_p(node,2) = 2;  % Dirichlet
            BCval_p(node,3) = BC.D.val(BC.D.tag == tag);
            continue
        end

        % ---- Neumann ----
        ineu = ismember(node_tags, BC.N.tag);
        if any(ineu)
            tag = node_tags(find(ineu,1));
            BCval_p(node,2) = 1;  % Neumann
            BCval_p(node,3) = BC.N.val(BC.N.tag == tag);
        else
            % default: homogeneous Neumann
            BCval_p(node,2) = 1;
            BCval_p(node,3) = 0;
        end
    end
end
