function lines = compute_fieldlines(msh, Ex_e, Ey_e, start_pts, opts)
% COMPUTE_FIELDLINES
% Forward-only, exact field-line tracing for element-wise constant E
% One single self-contained file

% ------------------------------------------------------------
% Defaults
% ------------------------------------------------------------
if nargin < 5, opts = struct(); end
if ~isfield(opts,'normalize'), opts.normalize = true; end
if ~isfield(opts,'smax'),      opts.smax = inf; end
if ~isfield(opts,'Ecrit'),     opts.Ecrit = 1e-14; end
if ~isfield(opts,'eps'),       opts.eps = 1e-14; end

p = msh.POS(:,1:2);
T = msh.TRIANGLES(:,1:3);
nt = size(T,1);

% ------------------------------------------------------------
% Enforce CCW orientation (MANDATORY)
% ------------------------------------------------------------
for t = 1:nt
    P = p(T(t,:),:);
    if det([P(2,:)-P(1,:); P(3,:)-P(1,:)]) < 0
        T(t,[2 3]) = T(t,[3 2]);
    end
end

% ------------------------------------------------------------
% Precompute geometry & adjacency
% ------------------------------------------------------------
geom = precompute_geometry(p,T);

% ------------------------------------------------------------
% Element-wise field
% ------------------------------------------------------------
E = [Ex_e(:), Ey_e(:)];
Emag = hypot(E(:,1),E(:,2));

if opts.normalize
    nz = Emag > opts.Ecrit;
    E(nz,:) = E(nz,:) ./ Emag(nz);
end

% ------------------------------------------------------------
% Locate starting elements
% ------------------------------------------------------------
TR = triangulation(T,p);
elem0 = pointLocation(TR,start_pts(:,1),start_pts(:,2));

valid = elem0 > 0 & ~isnan(elem0);
start_pts = start_pts(valid,:);
elem0 = elem0(valid);

Ns = size(start_pts,1);
lines = cell(Ns,1);

% ============================================================
% Trace each field line (forward only)
% ============================================================
for k = 1:Ns
    lines{k} = trace_single_line( ...
        start_pts(k,:).', elem0(k), geom, E, Emag, opts );
end

lines = lines(~cellfun(@isempty,lines));
end

% =====================================================================
% ======================= CORE TRACER =================================
% =====================================================================
function line = trace_single_line(x, elem, geom, E, Emag, opts)

x = x(:);
s = 0;

X = x.';
S = s;
ELEM = elem;

while true
    if elem <= 0 || Emag(elem) < opts.Ecrit || s >= opts.smax
        break
    end

    v = E(elem,:).';

    tmin = inf;
    hit  = 0;
    qhit = [];

    for k = 1:3
        a = geom.a{elem,k};
        b = geom.b{elem,k};
        n = geom.n{elem,k};
        e = b - a;

        denom = dot(n,v);
        if denom <= opts.eps
            continue   % not exiting through this edge
        end

        t = dot(n,a-x)/denom;
        if t <= opts.eps, continue; end

        q = x + t*v;
        proj = dot(q-a,e)/dot(e,e);

        if proj < -opts.eps || proj > 1+opts.eps
            continue
        end

        if t < tmin
            tmin = t;
            hit  = k;
            qhit = q;
        end
    end

    if hit == 0
        break
    end

    x = qhit;
    s = s + tmin;
    
    % store CURRENT element (still valid)
    X(end+1,:) = x.';
    S(end+1,1) = s;
    ELEM(end+1,1) = elem;
    
    % then move to neighbor
    elem = geom.adj(elem,hit);
end

line = struct( ...
    'x',X(:,1), ...
    'y',X(:,2), ...
    's',S, ...
    'elems',ELEM, ...
    'Ex',E(ELEM,1), ...
    'Ey',E(ELEM,2), ...
    'E', hypot(E(ELEM,1),E(ELEM,2)) );
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

% --- adjacency ---
for t = 1:nt
    for k = 1:3
        i1 = T(t,k);
        i2 = T(t,mod(k,3)+1);
        key = sprintf('%d_%d',min(i1,i2),max(i1,i2));

        if ~isKey(edge_map,key)
            edge_map(key) = [t k];
        else
            val = edge_map(key);
            t2 = val(1); k2 = val(2);
            geom.adj(t,k)   = t2;
            geom.adj(t2,k2)= t;
        end
    end
end

% --- edge geometry ---
for t = 1:nt
    P = p(T(t,:),:);
    c = mean(P,1);

    for k = 1:3
        a = P(k,:);
        b = P(mod(k,3)+1,:);
        e = b - a;

        n = [e(2), -e(1)];
        n = n / norm(n);

        mid = (a+b)/2;
        if dot(n,c-mid) > 0
            n = -n;
        end

        geom.a{t,k} = a(:);
        geom.b{t,k} = b(:);
        geom.n{t,k} = n(:);
    end
end
end
