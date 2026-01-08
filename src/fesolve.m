function [out] = fesolve(msh,BC,opts)
%FESOLVE Summary of this function goes here
%   Detailed explanation goes here

% Extract mesh info
tagbe = msh.LINES(:,3);

% Map each edge tag to flag/val
tagsAll = [BC.D.tag, BC.N.tag];
flagsAll = [2*ones(size(BC.D.tag)), 1*ones(size(BC.N.tag))];
valsAll = [BC.D.val, BC.N.val];
[~, idx] = ismember(tagbe, tagsAll);
if any(idx == 0)
    missing = unique(tagbe(idx == 0));
    error('Tags %s on boundary not defined in BC.', mat2str(missing));
end
BCflag_e = flagsAll(idx);
BCval_e  = valsAll(idx);
BC.BCflag_e = BCflag_e;
BC.BCval_e = BCval_e;

%% SOLVER CALL
ProblemKind = opts.ProblemKind; % [Electrostatic][Magnetostatic][QMagnetostaticSin][MagTimeDependent]

switch ProblemKind
    case {'Electrostatic','Magnetostatic','QMagnetostaticSin','QMagnetostaticSin_LAPL'}
        [out] = FEM2D00(msh, BC, opts);
    case {'MagTimeDependent'}
        p = msh.POS(:,1:2);
        t = msh.TRIANGLES(:,1:3);
        ireg = msh.TRIANGLES(:,4);
        edgeBound = msh.LINES(:,1:2); % lati (di elemento) sul contorno
        tagbe = msh.LINES(:,3);
        iregbe = ones(size(tagbe));
        switch opts.flag.decomp
            case (0)
                switch opts.flag.dt_auto
                    case (0)
                        [out] = FEM2D00_t(ProblemKind, opts.time_array, opts.sv, p, t, edgeBound, BCflag_e, BCval_e, ireg, iregbe, opts.materials, opts.source);
                    case (1)
                        [out] = FEM2D00_dt_auto(opts, p, t, edgeBound, BCflag_e, BCval_e, ireg, iregbe, opts.materials, opts.source);
                end
            case (1)
                [out] = FEM2D00_t_decomp(opts, p, t, edgeBound, BCflag_e, BCval_e, ireg, iregbe, opts.materials, opts.source);
        end
end

end

