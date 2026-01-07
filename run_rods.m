%% unit circle
% before executing: gmsh .\mesh\mesh_unit_circle.geo
utils_FEM;
% MESH
mesh_unit_circle; ndom = num_regions(msh);
% BC (Dirichlet)
opts.tag_boundary = 1; % domain boundary on edges marked with tag=1 in mesh files
% materials
[opts.materials] = set_materials('mesh_unit_circle',ndom);
% PROBLEM KIND
opts.ProblemKind = 'Electrostatic'; % [Electrostatic][Magnetostatic][QMagnetostaticSin][MagTimeDependent]
opts.source = 1;
% DIAGNOSTICS
opts.flag.print_measured_time = 0; 

% solution
[out] = fesolve(msh,opts);

out.field.phi = out.field.phi*8.8541878128E-12; % scaling
x = msh.POS(:,1); y = msh.POS(:,2); % get mesh coordinates
% Plot result
figure
trisurf(msh.TRIANGLES(:,1:3),x,y,out.field.phi,out.field.phi,edgecolor='none')
xlabel('x (m)'), ylabel('y (m)'); zlabel('solution - \phi'), axis tight;
ax = gca; ax.FontSize = 12;
f = gcf; colormap(f,ap.map.red_white_blue);

% Error with respect to analytical solution
figure
sol_e = (1 - x.^2 - y.^2)/4; % Whiteley p.52
trisurf(msh.TRIANGLES(:,1:3),x,y,(out.field.phi-sol_e),(out.field.phi-sol_e),edgecolor='none')
xlabel('x (m)'), ylabel('y (m)'); zlabel('error'); axis tight;
f = gcf; colormap(f,ap.map.red_white_blue);