function lines = compute_fieldlines(msh, Ex_e, Ey_e, start_pts, opts)
% COMPUTE_FIELDLINES
% Exact forward-only field-line tracing for piecewise-constant vector fields
% on triangular meshes.
%
% INPUTS
%   msh       : mesh struct with
%               - msh.POS       [np x 2 or np x 3]
%               - msh.TRIANGLES [nt x >=3]
%   Ex_e,Ey_e : element-wise constant field components [nt x 1]
%   start_pts : starting points [Ns x 2]
%   opts      : optional struct
%       .smax      : max geometric arc length (default: inf)
%       .Ecrit     : min |E| below which tracing stops (default: 1e-14)
%       .eps       : geometric tolerance (default: auto)
%       .maxSteps  : max crossed elements per line (default: 10*nt)
%
% OUTPUT
%   lines{k} with fields:
%       .x, .y     : polyline points
%       .s         : cumulative geometric arc length
%       .elems     : crossed elements
%       .Ex, .Ey   : physical field in crossed elements
%       .E         : physical |E| in crossed elements

    if nargin < 5
        opts = struct();
    end
    if ~isfield(opts,'smax'),     opts.smax = inf; end
    if ~isfield(opts,'Ecrit'),    opts.Ecrit = 1e-14; end
    if ~isfield(opts,'maxSteps'), opts.maxSteps = []; end

    p = msh.POS(:,1:2);
    T = msh.TRIANGLES(:,1:3);
    nt = size(T,1);

    % Tolleranza geometrica automatica
    if ~isfield(opts,'eps') || isempty(opts.eps)
        bbox = max(p,[],1) - min(p,[],1);
        hdom = max(bbox);
        if hdom <= 0
            hdom = 1;
        end
        opts.eps = 1e-12 * hdom;
    end

    if isempty(opts.maxSteps)
        opts.maxSteps = 10 * nt;
    end

    % ------------------------------------------------------------
    % Enforce CCW orientation
    % ------------------------------------------------------------
    for t = 1:nt
        P = p(T(t,:),:);
        A = [P(2,:) - P(1,:); P(3,:) - P(1,:)];
        if det(A) < 0
            T(t,[2 3]) = T(t,[3 2]);
        end
    end

    % ------------------------------------------------------------
    % Precompute geometry & adjacency
    % ------------------------------------------------------------
    geom = precompute_geometry(p,T);

    % ------------------------------------------------------------
    % Physical field and field direction
    % ------------------------------------------------------------
    Eraw = [Ex_e(:), Ey_e(:)];
    Emag = hypot(Eraw(:,1), Eraw(:,2));

    Edir = zeros(size(Eraw));
    nz = Emag > opts.Ecrit;
    Edir(nz,:) = Eraw(nz,:) ./ Emag(nz);

    % ------------------------------------------------------------
    % Locate starting elements
    % ------------------------------------------------------------
    TR = triangulation(T,p);
    elem0 = pointLocation(TR, start_pts(:,1), start_pts(:,2));

    valid = ~isnan(elem0) & (elem0 > 0);
    start_pts = start_pts(valid,:);
    elem0 = elem0(valid);

    Ns = size(start_pts,1);
    lines = cell(Ns,1);

    % ------------------------------------------------------------
    % Trace each field line
    % ------------------------------------------------------------
    for k = 1:Ns
        lines{k} = trace_single_line( ...
            start_pts(k,:).', elem0(k), geom, Eraw, Edir, Emag, opts);
    end

    lines = lines(~cellfun(@isempty,lines));
end


% =====================================================================
% ======================= CORE TRACER =================================
% =====================================================================
function line = trace_single_line(x0, elem0, geom, Eraw, Edir, Emag, opts)

    x = x0(:);
    elem = elem0;
    s = 0;

    X = x.';
    S = s;
    ELEM = elem;

    nSteps = 0;

    while true
        nSteps = nSteps + 1;
        if nSteps > opts.maxSteps
            warning('compute_fieldlines:maxSteps', ...
                    'Maximum number of crossed elements reached.');
            break
        end

        if elem <= 0 || Emag(elem) < opts.Ecrit || s >= opts.smax
            break
        end

        v = Edir(elem,:).';   % unit direction vector

        tmin = inf;
        hit  = 0;
        qhit = [];

        % Search the first exiting edge
        for k = 1:3
            a = geom.a{elem,k};
            b = geom.b{elem,k};
            n = geom.n{elem,k};
            e = b - a;

            denom = dot(n, v);

            % Must be exiting through this edge
            if denom <= opts.eps
                continue
            end

            % Time to intersection with the supporting line
            t = dot(n, a - x) / denom;
            if t <= opts.eps
                continue
            end

            % Intersection point
            q = x + t * v;

            % Check if q lies on the segment
            ee = dot(e,e);
            if ee <= 0
                continue
            end

            proj = dot(q - a, e) / ee;
            if proj < -opts.eps || proj > 1 + opts.eps
                continue
            end

            if t < tmin
                tmin = t;
                hit = k;
                qhit = q;
            end
        end

        if hit == 0 || isempty(qhit)
            break
        end

        % Truncate if smax is exceeded inside current element
        if s + tmin > opts.smax
            tleft = opts.smax - s;
            if tleft > 0
                x = x + tleft * v;
                s = opts.smax;
                X(end+1,:) = x.';
                S(end+1,1) = s;
                ELEM(end+1,1) = elem;
            end
            break
        end

        % Move to exit point
        x = qhit;
        s = s + tmin;

        % Store point still associated with current element
        X(end+1,:) = x.';
        S(end+1,1) = s;
        ELEM(end+1,1) = elem;

        % Move to neighboring element
        elem_next = geom.adj(elem, hit);
        if elem_next <= 0
            break
        end

        % Small push inside the next element to avoid sticking on the edge
        x = x + opts.eps * v;
        elem = elem_next;
    end

    line = struct( ...
        'x',    X(:,1), ...
        'y',    X(:,2), ...
        's',    S, ...
        'elems', ELEM, ...
        'Ex',   Eraw(ELEM,1), ...
        'Ey',   Eraw(ELEM,2), ...
        'E',    Emag(ELEM) );
end


% =====================================================================
% ===================== GEOMETRY SETUP ================================
% =====================================================================
function geom = precompute_geometry(p,T)

    nt = size(T,1);

    geom.adj = zeros(nt,3);
    geom.a   = cell(nt,3);
    geom.b   = cell(nt,3);
    geom.n   = cell(nt,3);

    edge_map = containers.Map('KeyType','char','ValueType','any');

    % ------------------------------------------------------------
    % Adjacency
    % Edge k of triangle t is from local node k to local node mod(k,3)+1
    % ------------------------------------------------------------
    for t = 1:nt
        for k = 1:3
            i1 = T(t,k);
            i2 = T(t,mod(k,3)+1);

            key = sprintf('%d_%d', min(i1,i2), max(i1,i2));

            if ~isKey(edge_map, key)
                edge_map(key) = [t, k];
            else
                val = edge_map(key);
                t2 = val(1);
                k2 = val(2);

                geom.adj(t, k) = t2;
                geom.adj(t2, k2) = t;
            end
        end
    end

    % ------------------------------------------------------------
    % Edge geometry and outward normals
    % ------------------------------------------------------------
    for t = 1:nt
        P = p(T(t,:), :);
        c = mean(P,1);

        for k = 1:3
            a = P(k,:);
            b = P(mod(k,3)+1,:);
            e = b - a;

            ne = norm(e);
            if ne <= 0
                error('Degenerate edge found in triangle %d.', t);
            end

            % Candidate normal
            n = [e(2), -e(1)];
            n = n / norm(n);

            % Make it outward
            mid = (a + b) / 2;
            if dot(n, c - mid) > 0
                n = -n;
            end

            geom.a{t,k} = a(:);
            geom.b{t,k} = b(:);
            geom.n{t,k} = n(:);
        end
    end
end