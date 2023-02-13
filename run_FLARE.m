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

%% Example: multiple_regions time-harmonic 50 Hz
% before executing: gmsh .\mesh\mesh_ex_regions.geo
utils_FEM;
% MESH
mesh_ex_regions; ndom = num_regions(msh);
% BC (Dirichlet)
opts.tag_boundary = [5,6,7,8];
% materials
[opts.materials] = set_materials('mesh_ex_regions',ndom);
% PROBLEM KIND
opts.ProblemKind = 'QMagnetostaticSin'; % [Electrostatic][Magnetostatic][QMagnetostaticSin][MagTimeDependent]
opts.source = [0, 1e6, -1e6, 0]; % j0 = 1A/mm^2
opts.freq = 50;
% DIAGNOSTICS
opts.flag.print_measured_time = 0; 

% solution
[out] = fesolve(msh,opts);

x = msh.POS(:,1); y = msh.POS(:,2); z = msh.POS(:,3); % get mesh coordinates
% Plot results
figure
trisurf(msh.TRIANGLES(:,1:3),x,y,z,abs(out.field.Az)); view(2);
xlabel('x (m)'), ylabel('y (m)'); zlabel('Az'), colorbar;
f = gcf; colormap(f,ap.map.red_white_blue); 

disp(table((1:4)',opts.materials,abs(out.scal.I),angle(out.scal.I),'VariableNames',{'Region','Material','I (A)','phase(I) (rad)'}));

%% Corridor, time-harmonic 50 Hz
% before executing: gmsh .\mesh\mesh_corridor.geo
utils_FEM;
% MESH
mesh_corridor; ndom = num_regions(msh);
% BC (Dirichlet)
opts.tag_boundary = 71:74;
% materials
[opts.materials] = set_materials('mesh_corridor',ndom);

% PROBLEM KIND
opts.ProblemKind = 'QMagnetostaticSin'; % [Electrostatic][Magnetostatic][QMagnetostaticSin][MagTimeDependent]
opts.source = [1*3.5d7/8.6165, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; % J0 = 3.5d7/8.6160 to get I_phase1 = 250 A
opts.freq = 50;
% DIAGNOSTICS
opts.flag.print_measured_time = 0; 

% solution
[out] = fesolve(msh,opts);

disp(table((1:7)',opts.materials(1:7),abs(out.scal.I(1:7)),'VariableNames',{'Region','Material','I (A)'}));

% Plot Az and field lines of B
close all
x = msh.POS(:,1); y = msh.POS(:,2); z = msh.POS(:,3);

t_plot = 0.2+10E-3; w = 2*pi*opts.freq;
vec = out.field.Az; 
Azt = abs(vec).*sin(w*t_plot+angle(vec));
trisurf(msh.TRIANGLES(:,1:3),x,y,z,Azt,edgecolor='none');
xplot = [-5 6]; yplot = [10 18];
xlim(xplot); ylim(yplot);
view(2); colorbar;
xlabel('x (m)'); ylabel('y (m)');

np = 1E3;
xq = linspace(xplot(1),xplot(2),np);
yq = linspace(yplot(1),yplot(2),np);

vec = out.field.Bx; Bxt = abs(vec).*sin(w*t_plot+angle(vec));
vec = out.field.By; Byt = abs(vec).*sin(w*t_plot+angle(vec));
Bx = tri2grid(msh.POS(:,1:2)',msh.TRIANGLES(:,1:3)',Bxt,xq,yq);
By = tri2grid(msh.POS(:,1:2)',msh.TRIANGLES(:,1:3)',Byt,xq,yq);

l = streamslice(xq,yq,Bx,By,1); % change 1 to 2 to double the streamline density
set(l,'LineWidth',0.5);
set(l,'Color','k');
f = gcf; colormap(f,ap.map.red_white_blue); % shading interp;

%% Corridor, sinusoidal, transient (constant dt)
% before executing: gmsh .\mesh\mesh_corridor.geo
utils_FEM;
% MESH
mesh_corridor; ndom = num_regions(msh);
% BC (Dirichlet)
opts.tag_boundary = 71:74;
% materials
[opts.materials] = set_materials('mesh_corridor',ndom);
% PROBLEM KIND
opts.ProblemKind = 'MagTimeDependent'; % [Electrostatic][Magnetostatic][QMagnetostaticSin][MagTimeDependent]
time_fin = 0.25;
opts.time_array = [linspace(0,time_fin,1E3)]; % [0, logspace(-8,-3,20), linspace(2E-3,time_fin,100)];
opts.source = @(time) [1*3.5d7*sin(2*pi*50*time), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
% SAVE
opts.sv.nit_skip = 1; % variables saved each <nit_skip> time iterations
% DECOMPOSITION
opts.flag.decomp = 1; % 0: constant dt; 1: variable dt

% solution
[out] = fesolve(msh,opts);

% post processing
ip = out.scal.I(7,:);
figure
plot(out.sv.time,out.scal.I(1,:)./8.6160,'k-',...
    out.sv.time,out.scal.I(2,:)./8.6160,'r-.',...
    out.sv.time,out.scal.I(3,:)./8.6160,'c--',...
    out.sv.time,2*ip./8.6160,'b:','LineWidth',1.2); 
grid on; lgd=legend('i_{Phase 1}(t)','I_{Phase 2}(t)','I_{Phase 3}(t)','2 \cdot I_{pipe}(t)','Location','northeast');
lgd.NumColumns = 2;
xlabel('t (s)');
ylabel('Current (A)');

% amplitude of waveforms
I_scaled = out.scal.I(:,end-npp:end)/8.6062;

npp = 50;
A_iph1 = amplitude(I_scaled(1,:));
A_iph2 = amplitude(I_scaled(2,:));
A_iph3 = amplitude(I_scaled(3,:));
A_p = amplitude(I_scaled(7,:));

disp(table([A_iph1;A_iph2;A_iph3;A_p],'RowNames',{'Line 1','Line 2','Line 3','Pipe'},'VariableNames',{'Current (A)'}));

%% Inverse Laplace
% before executing: gmsh .\mesh\mesh_corridor.geo
utils_FEM
addpath("796_matlab\")
help lapinv_796

% HANDLE FUNCTION of source term, for FEM(s)
omega = 2*pi*50;
fun_sorg = @(s) laplace_sine(s,omega);

% HANDLE FUNCTION of FEM solver, for lapinv796
fun_FEM = @(s) ilaplace(s,fun_sorg);

tarray = linspace(1E-6,0.25,50); % % linspace(0,4*pi,250);
tol = 1E-6; ssbar = 0; nmax = 550; sigma0 = 0.; % nmax = 550; tol = 10^-5;
[tarray,fzinv,~,~] = lapinv_796(fun_FEM,tarray,tol,ssbar,nmax,sigma0);

% Post processingt
ip = fzinv/8.6160; % scaling
A_ip_LapInv = (max(ip)-min(ip))/2;

plot(tarray,ip,'-x','LineWidth',1.2); 
grid on; legend('I_{pipe}(t) - Inverse Laplace');
xlabel('t (s)');
ylabel('Current (A)');
