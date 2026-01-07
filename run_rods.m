clear variables
close all

%% facing rods
% before executing: gmsh .\mesh\mesh_rods_half.geo
utils_FEM;
% MESH
mesh_rods_half; ndom = num_regions(msh);
% BC (Dirichlet)
BC.D.tag = [11,12];
BC.D.val = [1,0]*250;
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

% Plot potential
figure
trisurf(msh.TRIANGLES(:,1:3),x,y,out.field.phi,out.field.phi,edgecolor='none')
xlabel('x (m)'), ylabel('y (m)'); zlabel('solution - \phi'), axis tight; axis equal; view(2); colorbar;
ax = gca; ax.FontSize = 12;
f = gcf; colormap(f,ap.map.red_white_blue);

figure
quiver(msh.POS(:,1),msh.POS(:,2),out.field.Ex,out.field.Ey,'r')
