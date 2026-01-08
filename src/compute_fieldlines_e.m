function lines = compute_fieldlines_e(msh, Ex_e, Ey_e, start_pts, opts)
%COMPUTE_FIELDLINES_E  Compute field lines using element-wise constant E-field
%
% lines = compute_fieldlines_e(msh, Ex_e, Ey_e, start_pts, opts)
%
% INPUT:
%   msh.POS(:,1:2)   node coordinates [np x 2]
%   msh.TRIANGLES    element connectivity [nt x 3]
%   Ex_e, Ey_e       element-wise E-field components [nt x 1]
%   start_pts        [Ns x 2] starting points
%   opts (optional)  struct with fields:
%       .ds          max integration step (default: char. size/100)
%       .smax        max arc length per direction (default: 1)
%       .normalize   true/false (default: true) → integrate direction only
%       .direction   'forward', 'backward', or 'both' (default: 'forward')
%       .Ecrit       critical field magnitude to stop (default: 1e-12)
%
% OUTPUT:
%   lines{k}.x,y,s   coordinates and arc-length parameter
%   lines{k}.Ex,Ey,E field values along the line (element-wise constant)
%   lines{k}.elems   element indices along the line (for debugging)

% ---------------- defaults ----------------
if nargin < 5, opts = struct(); end
if ~isfield(opts,'normalize'), opts.normalize = true; end
if ~isfield(opts,'smax'),      opts.smax = 1; end
if ~isfield(opts,'direction'), opts.direction = 'forward'; end
if ~isfield(opts,'Ecrit'),     opts.Ecrit = 1e-12; end

% Characteristic mesh size (min edge length)
p = msh.POS(:,1:2);
T = msh.TRIANGLES;
edge_lengths = [
    sqrt(sum((p(T(:,1),:) - p(T(:,2),:)).^2,2));
    sqrt(sum((p(T(:,2),:) - p(T(:,3),:)).^2,2));
    sqrt(sum((p(T(:,3),:) - p(T(:,1),:)).^2,2))];
h = min(edge_lengths);
if ~isfield(opts,'ds'), opts.ds = h/100; end

% Triangulation object for point location
TR = triangulation(T(:,1:3), p);

% ---------------- ODE setup ----------------
% Field evaluator (element-wise constant)
F_E = @(x,y) eval_E_at_point(x, y, TR, Ex_e, Ey_e);

% ODE function handles both normalization and direction sign
odefun = @(s,x,direction_sign) field_ode(x, F_E, opts.normalize, direction_sign);

% Events: 1) exit mesh, 2) critical point (|E| small)
eventfun = @(s,x) events_fieldline(x, F_E, TR, opts.Ecrit);

odeopts = odeset(...
    'MaxStep', opts.ds, ...
    'Events', eventfun, ...
    'RelTol', 1e-6, ...
    'AbsTol', 1e-9);

% ---------------- integrate all field lines ----------------
Ns = size(start_pts,1);
lines = cell(Ns,1);

% Directions to integrate
directions = {'forward'};
if strcmp(opts.direction, 'both')
    directions = {'forward', 'backward'};
elseif strcmp(opts.direction, 'backward')
    directions = {'backward'};
end

for k = 1:Ns
    x0 = start_pts(k,:).';
    
    % Check if start point is inside mesh
    loc = pointLocation(TR, x0(1), x0(2));
    if loc < 0 || isnan(loc) || loc == 0
        warning('Start point (%.3f, %.3f) outside mesh. Skipping.', x0(1), x0(2));
        continue;
    end
    
    % Integrate for each direction
    all_segments = struct('x',{}, 'y',{}, 's',{});
    for d = 1:numel(directions)
        dir_str = directions{d};
        dir_sign = strcmp(dir_str, 'backward') * -2 + 1; % +1 forward, -1 backward
        
        % Integrate in this direction
        [s_seg, X_seg, te, xe, ie] = ode45(@(s,x) odefun(s,x,dir_sign), [0 opts.smax], x0, odeopts);
        
        % Compute arc-length parameter (starting from 0 at x0)
        x_path = X_seg(:,1);
        y_path = X_seg(:,2);
        ds = hypot(diff(x_path), diff(y_path));
        s_path = [0; cumsum(ds)];
        
        % Store segment
        seg = struct();
        seg.x = x_path;
        seg.y = y_path;
        seg.s = s_path;
        all_segments(d) = seg;
    end
    
    % Combine segments if bidirectional
    if numel(directions) == 2
        % Flip backward segment (so s increases from x0)
        bw = all_segments(2);
        bw.x = flipud(bw.x(2:end)); % skip duplicate x0
        bw.y = flipud(bw.y(2:end));
        bw.s = flipud(bw.s(2:end));
        bw.s = max(bw.s) - bw.s; % make s increasing from x0
        
        % Concatenate: backward (reversed) + forward
        x_comb = [bw.x; all_segments(1).x];
        y_comb = [bw.y; all_segments(1).y];
        s_comb = [bw.s; all_segments(1).s + bw.s(end)];
        
        % Evaluate field along combined path
        [Ex_path, Ey_path] = F_E(x_comb, y_comb);
        E_path = hypot(Ex_path, Ey_path);
        elems_path = pointLocation(TR, x_comb, y_comb);
        
        lines{k}.x = x_comb;
        lines{k}.y = y_comb;
        lines{k}.s = s_comb;
        lines{k}.Ex = Ex_path;
        lines{k}.Ey = Ey_path;
        lines{k}.E = E_path;
        lines{k}.elems = elems_path;
    else
        seg = all_segments(1);
        [Ex_path, Ey_path] = F_E(seg.x, seg.y);
        E_path = hypot(Ex_path, Ey_path);
        elems_path = pointLocation(TR, seg.x, seg.y);
        
        lines{k}.x = seg.x;
        lines{k}.y = seg.y;
        lines{k}.s = seg.s;
        lines{k}.Ex = Ex_path;
        lines{k}.Ey = Ey_path;
        lines{k}.E = E_path;
        lines{k}.elems = elems_path;
    end
end

% Remove empty lines
lines = lines(~cellfun(@isempty, lines));
end

% ==========================================================
function [Ex, Ey] = eval_E_at_point(x, y, TR, Ex_e, Ey_e)
% Evaluate element-wise constant field at points (x,y)
if isscalar(x)
    x = x(:); y = y(:);
end

% Find element for each point (-1 if outside)
loc = pointLocation(TR, x, y);

% Initialize outputs
Ex = NaN(size(x));
Ey = NaN(size(x));

% Valid points (inside mesh)
valid = (loc > 0) & ~isnan(loc);
Ex(valid) = Ex_e(loc(valid));
Ey(valid) = Ey_e(loc(valid));
end

% ==========================================================
function dxds = field_ode(x, F_E, normalize, direction_sign)
% ODE RHS: dx/ds = direction_sign * (E / |E|) if normalize, else direction_sign * E
x = x(:);
[Ex, Ey] = F_E(x(1), x(2));

% Handle outside/invalid points
if isnan(Ex) || isnan(Ey)
    dxds = [0; 0];
    return;
end

E_vec = [Ex; Ey];
E_mag = norm(E_vec);

% Critical point: field too small → stop moving
if E_mag < eps
    dxds = [0; 0];
    return;
end

% Normalize if requested
if normalize
    E_vec = E_vec / E_mag;
end

% Apply direction sign and return
dxds = direction_sign * E_vec;
end

% ==========================================================
function [value, isterminal, direction] = events_fieldline(x, F_E, TR, Ecrit)
% Events: 
%   (1) Exit mesh (point not in any element)
%   (2) Critical point (|E| < Ecrit)
x = x(:);

% Event 1: exit mesh
loc = pointLocation(TR, x(1), x(2));
exit_mesh = (loc <= 0) | isnan(loc);

% Event 2: critical point
[Ex, Ey] = F_E(x(1), x(2));
E_mag = hypot(Ex, Ey);
crit_point = (E_mag < Ecrit);

% Combine events
value = [exit_mesh; crit_point];
isterminal = [1; 1]; % Stop integration for both
direction = [0; 0];  % Trigger from any direction
end