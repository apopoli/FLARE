function lines = compute_fieldlines(msh, Ex, Ey, start_pts, opts)
%COMPUTE_FIELDLINES  Compute electric field lines on a 2D FEM mesh
%
% lines = compute_fieldlines(msh, Ex, Ey, start_pts, opts)
%
% INPUT:
%   msh.POS(:,1:2)   node coordinates
%   Ex, Ey           nodal electric field components
%   start_pts        [Ns x 2] starting points
%   opts (optional)  struct with fields:
%       .ds          max integration step (default: characteristic size/100)
%       .smax        max arc length (default: 1)
%       .normalize   true/false (default: true)
%
% OUTPUT:
%   lines{k}.x       x coordinates
%   lines{k}.y       y coordinates
%   lines{k}.s       curvilinear coordinate
%   lines{k}.Ex,Ey,E field values along the line

% ---------------- defaults ----------------
if nargin < 5, opts = struct(); end
if ~isfield(opts,'normalize'), opts.normalize = true; end
if ~isfield(opts,'smax'),      opts.smax = 1; end

% characteristic mesh size
p = msh.POS(:,1:2);
h = mean(sqrt(sum((p(msh.TRIANGLES(:,1),:) - ...
                   p(msh.TRIANGLES(:,2),:)).^2,2)));
if ~isfield(opts,'ds'), opts.ds = h/100; end

% ---------------- FEM interpolants ----------------
Fx = scatteredInterpolant(p(:,1), p(:,2), Ex, 'linear', 'none');
Fy = scatteredInterpolant(p(:,1), p(:,2), Ey, 'linear', 'none');

% ---------------- ODE definition ----------------
odefun = @(s,x) field_ode(x, Fx, Fy, opts.normalize);

odeopts = odeset( ...
    'MaxStep', opts.ds, ...
    'Events', @(s,x) stop_outside(x, Fx) );

% ---------------- integrate all field lines ----------------
Ns = size(start_pts,1);
lines = cell(Ns,1);

for k = 1:Ns
    x0 = start_pts(k,:).';
    if isnan(Fx(x0(1),x0(2)))
        continue
    end

    [s,X] = ode45(odefun, [0 opts.smax], x0, odeopts);

    x = X(:,1);
    y = X(:,2);

    Exl = Fx(x,y);
    Eyl = Fy(x,y);
    El  = hypot(Exl,Eyl);

    % arc-length
    ds = hypot(diff(x),diff(y));
    sline = [0; cumsum(ds)];

    lines{k}.x  = x;
    lines{k}.y  = y;
    lines{k}.s  = sline;
    lines{k}.Ex = Exl;
    lines{k}.Ey = Eyl;
    lines{k}.E  = El;
end

end

% ==========================================================
function dxds = field_ode(x, Fx, Fy, normalize)
Ex = Fx(x(1),x(2));
Ey = Fy(x(1),x(2));

if isnan(Ex) || isnan(Ey)
    dxds = [0;0];
    return
end

if normalize
    n = hypot(Ex,Ey);
    if n == 0
        dxds = [0;0];
    else
        dxds = [Ex;Ey] / n;
    end
else
    dxds = [Ex;Ey];
end
end

% ==========================================================
function [value,isterminal,direction] = stop_outside(x, Fx)
value = ~isnan(Fx(x(1),x(2)));  % 0 when outside
isterminal = 1;
direction = 0;
end
