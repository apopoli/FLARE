function [out] = fesolve(msh,opts)
%FESOLVE Summary of this function goes here
%   Detailed explanation goes here

%% MESH
p0 = msh.POS; p = p0(:,1:2);
t0 = msh.TRIANGLES; t = t0(:,1:3);
ireg = t0(:,4);
ed = msh.LINES;
edgeBound = ed(:,1:2); % lati (di elemento) sul contorno
tagbe = ed(:,3);

%% BCs
bcflag_e = zeros(length(edgeBound),1);
bcval_e = zeros(length(edgeBound),1);

for i=1:length(opts.tag_boundary)
    ibe = find(tagbe==opts.tag_boundary(i));  % trova tutti i lati sul contorno con tag (physical group assegnato in gmsh)
    if ~any(ibe)
        error('Controllare lati su contorno')
    end
    % 1: Neumann 2: Dirichlet
    bcflag_e(ibe) = 2 * ones(length(ibe),1);  % assegna il flag per la condizione al controno 
    bcval_e(ibe) = 0 * ones(length(ibe),1);  % assegna il valore per la condizione al controno
end

iregbe = ones(size(tagbe));

%% SOLVER CALL
ProblemKind = opts.ProblemKind; % [Electrostatic][Magnetostatic][QMagnetostaticSin][MagTimeDependent]

switch ProblemKind
    case {'Electrostatic','Magnetostatic','QMagnetostaticSin','QMagnetostaticSin_LAPL'}
        [out] = FEM2D00(opts, p, t, edgeBound, bcflag_e, bcval_e, ireg, iregbe, opts.materials, opts.source);
    case {'MagTimeDependent'}
        switch opts.flag.decomp
            case (0)
                switch opts.flag.dt_auto
                    case (0)
                        [out] = FEM2D00_t(ProblemKind, opts.time_array, opts.sv, p, t, edgeBound, bcflag_e, bcval_e, ireg, iregbe, opts.materials, opts.source);
                    case (1)
                        [out] = FEM2D00_dt_auto(opts, p, t, edgeBound, bcflag_e, bcval_e, ireg, iregbe, opts.materials, opts.source);
                end
            case (1)
                [out] = FEM2D00_t_decomp(opts, p, t, edgeBound, bcflag_e, bcval_e, ireg, iregbe, opts.materials, opts.source);
        end
end

end

