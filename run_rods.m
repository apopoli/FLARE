clear variables
close all

% --- Parametri gas / breakdown ---
P     = 1 * 101325;       % Pa  (esempio: 0.1 bar)
Tgas  = 300;                % K
gamma = 0.01;               % coeff. emissione secondaria
kB    = 1.380649e-23;       % J/K
Ngas  = P/(kB*Tgas);        % densità numerica [m^-3]

% --- Dati BOLSIG+ per la miscela in questione ---
bolsigFile = "11_out_He.dat";   % 04_out_Air

% MESH
mesh_rods_half_rsmall;

% facing rods
% before executing: gmsh .\mesh\mesh_rods_half.geo
utils_FEM;
ndom = num_regions(msh);
% BC (Dirichlet)
BC.D.tag = [11,12];
BC.D.val = [1,0]*1;
BC.N.tag = 13;
BC.N.val = 0;
% materials
[opts.materials] = set_materials('mesh_unit_circle',ndom); % borrow from unit_circle example
% PROBLEM KIND
opts.ProblemKind = 'Electrostatic'; % [Electrostatic][Magnetostatic][QMagnetostaticSin][MagTimeDependent]
opts.source = 0;
% DIAGNOSTICS
opts.flag.print_measured_time = 0; 

% solution
[out] = fesolve(msh,BC,opts);

x = msh.POS(:,1); y = msh.POS(:,2); % get mesh coordinates

%% field lines_e
BCval_p = out.BCval_p;
idx = find(BCval_p(:,2)==2 & BCval_p(:,3)==1);
start_pts = msh.POS(idx(:),1:2)+1E-10; % start_pts = msh.POS(idx(1:10:end),1:2);

lines = compute_fieldlines(msh,out.field_e.Ex,out.field_e.Ey,start_pts);

figure;
% Plot come patch 2D (z = 0)
p = patch('Faces', msh.TRIANGLES(:,1:3), ...
          'Vertices', [msh.POS(:,1:2), zeros(size(msh.POS,1),1)], ... % z=0
          'FaceVertexCData', out.field.phi, ...
          'FaceColor', 'interp', ...
          'EdgeColor', 'w', ...
          'CDataMapping', 'scaled');
colorbar;
view(2); axis equal tight;
xlabel('x (m)'); ylabel('y (m)');

hold on;
for k = 1:length(lines)
    if isempty(lines{k}), continue; end
    plot(lines{k}.x, lines{k}.y, 'k-', 'LineWidth', 1); % bianco per contrasto
end
hold off;

figure
trisurf(msh.TRIANGLES(:,1:3),msh.POS(:,1),msh.POS(:,2),out.field.phi,'edgecolor','none'); view(2); colorbar; axis equal;

figure
trisurf(msh.TRIANGLES(:,1:3),msh.POS(:,1),msh.POS(:,2),sqrt(out.field.Ex.^2+out.field.Ey.^2),'edgecolor','none'); view(2); colorbar; axis equal;
xlabel('x (m)'); ylabel('y (m)'); title('|E| (V/m)')
field_unif_theory = 1/(0.006); % V/d

% field lines_e
BCval_p = out.BCval_p;
idx = find(BCval_p(:,2)==2 & BCval_p(:,3)==1);
start_pts = msh.POS(idx(:),1:2)+1E-10; % o quello che usavi tu [2:13]

lines = compute_fieldlines(msh, out.field_e.Ex, out.field_e.Ey, start_pts, []);

% ==============================================================
%   Townsend non-uniforme: Vb per ogni linea di campo
% ==============================================================
B = read_Bolsig(bolsigFile);

EN_grid   = B.Transp.E_N_Td;  % [Td]
alphaN_g  = B.Transp.A18_Townsend_ioniz_coef_alphaN_m2;  % [m^2]
if isfield(B.Transp,'A19_Townsend_attach_coef_etaN_m2')
    etaN_g = B.Transp.A19_Townsend_attach_coef_etaN_m2;  % [m^2]
else
    etaN_g = zeros(size(alphaN_g));
end

% interpolanti αN(E/N), ηN(E/N). Fuori range → 0
alphaN_fun = @(EN_Td) interp1(EN_grid, alphaN_g, EN_Td, 'pchip', 0);
etaN_fun   = @(EN_Td) interp1(EN_grid, etaN_g, EN_Td, 'pchip', 0);

logFactor = log(1 + 1/gamma);        % lato destro del criterio Townsend
nLines    = numel(lines);

Vb_line   = nan(nLines,1);           % breakdown voltage per linea [V]
I_at_Vb   = nan(nLines,1);           % integrale alla soluzione
has_root  = false(nLines,1);         % flag: breakdown trovato

% Limiti di ricerca per Vb (puoi regolarli)
Vmin = 0;           % sempre F(0) < 0 (integrale = 0)
Vmax_global = 20e3;  % [V] tensione massima che sei disposto a considerare

for k = 1:nLines
    % Coordinate lungo la linea
    X = lines{k}.x(:);
    Y = lines{k}.y(:);
    E0_line = lines{k}.E(:);   % |E| per 1 V in ogni punto

    if numel(X) < 2 || all(E0_line == 0)
        continue;
    end

    % Lunghezza geometrica dei segmenti (questa è davvero ds in metri)
    dX = diff(X);
    dY = diff(Y);
    ds = hypot(dX, dY);        % [m]

    % Valore di E per il segmento (uso il valore nel punto iniziale)
    E0_seg = E0_line(1:end-1); % [V/m] per 1 V
    % Funzione F(V) = ∫ α_eff(x; V) dx - logFactor
    F = @(V) nonuniform_integral_minus_log( ...
                V, ds, E0_seg, Ngas, alphaN_fun, etaN_fun, logFactor);

    % Cerchiamo un intervallo [0, Vhigh] con cambiamento di segno
    Vhigh = 1e3;  % partiamo da 1 kV
    F0    = F(0); % dovrebbe essere -logFactor < 0

    if F0 > 0
        % Caso patologico (non dovrebbe succedere)
        warning('Linea %d: F(0) > 0, qualcosa non torna.', k);
        continue;
    end

    Fhigh = F(Vhigh);

    while Fhigh <= 0 && Vhigh < Vmax_global
        Vhigh = 2*Vhigh;
        Fhigh = F(Vhigh);
    end

    if Fhigh <= 0
        % Nessun breakdown fino a Vmax_global
        Vb_line(k) = NaN;
        has_root(k) = false;
        I_at_Vb(k) = F(Vmax_global) + logFactor; % integrale a Vmax (per info)
        continue;
    end

    % Ora F(0) < 0 e F(Vhigh) > 0 → intervallo valido per fzero
    try
        Vb = fzero(F, [0, Vhigh]);
        Vb_line(k)  = Vb;
        has_root(k) = true;
        I_at_Vb(k)  = F(Vb) + logFactor;   % dovrebbe essere ≈ logFactor
    catch ME
        warning('Linea %d: fzero fallito (%s).', k, ME.message);
        Vb_line(k)  = NaN;
        has_root(k) = false;
    end
end

% Breakdown globale del sistema
if any(has_root)
    idx_valid = find(has_root);                 % indici delle linee valide
    [Vb_system,i_min] = min(Vb_line(idx_valid)); 
    idx_crit = idx_valid(i_min);                % indice reale della linea
else
    Vb_system = NaN;
    idx_crit = NaN;
end

fprintf('\n=== NON-UNIFORM BREAKDOWN ===\n');
for k = 1:nLines
    if has_root(k)
        fprintf('Line %3d: Vb = %8.5f kV\n', k, Vb_line(k)*1e-3);
    else
        fprintf('Line %3d: Vb > %8.5f kV (nessun breakdown trovato)\n', ...
                k, Vmax_global*1e-3);
    end
end
if ~isnan(Vb_system)
    fprintf('\nGlobal Breakdown ~ %.4f kV (min Vb_line)\n', Vb_system*1e-3);
else
    fprintf('\nNo line is in breakdown until %.1f kV\n', Vmax_global*1e-3);
end

function F = nonuniform_integral_minus_log(V, ds, E0_seg, Ngas, alphaN_fun, etaN_fun, logFactor)
% NONUNIFORM_INTEGRAL_MINUS_LOG
%   Calcola F(V) = ∫ α_eff(x; V) dx - log(1 + 1/gamma)
%
%   ds      : lunghezze dei segmenti lungo la linea [m]
%   E0_seg  : |E| sui segmenti per 1 V di differenza [V/m]
%   Ngas    : densità numerica del gas [m^-3]
%   alphaN_fun, etaN_fun: interpolanti αN(E/N), ηN(E/N) -> [m^2]
%   logFactor: log(1 + 1/gamma) (costante)

    % Campo reale lungo i segmenti per la tensione V
    E_seg = E0_seg * V;          % [V/m]

    % E/N e E/N in Td
    EN    = E_seg ./ Ngas;       % [V·m^2]
    EN_Td = EN / 1E-21;          % [Td]

    % αN(E/N), ηN(E/N)
    alphaN = alphaN_fun(EN_Td);  % [m^2]
    etaN   = etaN_fun(EN_Td);    % [m^2]

    % α_eff(x) = (αN - ηN) * N  → [m^2]
    alpha_eff = (alphaN - etaN) * Ngas;

    % integrale di Townsend lungo la linea (somma a tratti)
    I = sum(alpha_eff .* ds);    % adimensionale

    % F(V) = I - log(1+1/gamma)
    F = I - logFactor;
end

% ==============================================================
%   Plot linee a Vb_system: rosso se breakdown, nero altrimenti
% ==============================================================

figure;
patch('Faces', msh.TRIANGLES(:,1:3), ...
      'Vertices', [msh.POS(:,1:2), zeros(size(msh.POS,1),1)], ...
      'FaceVertexCData', sqrt(out.field.Ex.^2+out.field.Ey.^2), ...
      'FaceColor', 'interp', ...
      'EdgeColor', 'w', ...
      'CDataMapping', 'scaled');
colorbar;
view(2); axis equal tight;
xlabel('x (m)'); ylabel('y (m)');
f = gcf; % colormap(f, ap.map.red_white_blue);
hold on;

for k = 1:nLines
    if isempty(lines{k}), continue; end

    if k == idx_crit
        plot(lines{k}.x, lines{k}.y, 'r-', 'LineWidth', 2);
    else
        plot(lines{k}.x, lines{k}.y, 'k-', 'LineWidth', 1);
    end
end

hold off;
