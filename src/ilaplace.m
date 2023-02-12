function [I_FEM] = ilaplace(s,handle_sorgente)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

utils_FEM;
% MESH
mesh_corridor; ndom = num_regions(msh);
% BC (Dirichlet)
opts.tag_boundary = 71:74;
% materials
opts.materials = strings(ndom,1);
opts.materials(1:3) = 'aluminium'; % phase_1, 2 e 3
opts.materials(4) = 'air'; % OGW_1
opts.materials(5) = 'air'; % Aria
opts.materials(6) = 'soil'; % Terreno
opts.materials(7) = 'aMagSteel'; % Pipe
opts.materials(8) = 'air'; % Pipe centro
opts.materials(9:12) = 'soil'; % Mit 1, 2, 3, 4
% PROBLEM KIND
opts.ProblemKind = 'QMagnetostaticSin_LAPL'; % [Electrostatic][Magnetostatic][QMagnetostaticSin][MagTimeDependent]
opts.s_lapl = s;
opts.source = handle_sorgente;
% DIAGNOSTICA
opts.flag.print_measured_time = 0;

opts.source = opts.source(opts.s_lapl);

% solution
[out] = fesolve(msh,opts);

I_FEM = out.scal.I(7); % I_p

% fprintf('s = %g + %gi \n',real(s),imag(s));
% fprintf('I_ph1 %g \n',abs(out.scal.I(1)))
% fprintf('I_terr %g \n',abs(out.scal.I(6)))
% fprintf('I_p %g \n\n',abs(I_FEM))

end