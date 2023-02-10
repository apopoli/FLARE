%% unit circle
utils_FEM;
% MESH
unit_circle; ndom = num_regions(msh);
% BC (Dirichlet)
opts.tag_boundary = 1;
% MATERIALS
[opts.materials] = set_materials('unit_circle',ndom);
% PROBLEM KIND
opts.ProblemKind = 'Electrostatic'; % [Electrostatic][Magnetostatic][QMagnetostaticSin][MagTimeDependent]
opts.source = 1;
% DIAGNOSTICA
opts.flag.print_measured_time = 0; 

% solution
tic
[out] = fesolve(msh,opts);
toc
out.field.phi = out.field.phi*8.8541878128E-12;

x = msh.POS(:,1); y = msh.POS(:,2); z = msh.POS(:,3);

figure
trisurf(msh.TRIANGLES(:,1:3),x,y,z,out.field.phi,edgecolor='none')
view(2); axis equal; colorbar; title('FLARE');

sol_e = (1 - x.^2 - y.^2)/4; % Whiteley p.52 https://it.mathworks.com/help/pde/ug/solve-poissons-equation-on-a-unit-disk.html
figure
trisurf(msh.TRIANGLES(:,1:3),x,y,z,sol_e,edgecolor='none')
view(2); axis equal; colorbar; title('EXACT');

%% POISSON 2D
close all
figure
trisurf(msh.TRIANGLES(:,1:3),x,y,z,out.field.phi,edgecolor='none')
view(2); axis equal; colorbar; xlabel('x (m)'), ylabel('y (m)'); axis tight;
ax = gca; ax.FontSize = 12;

% interpolazione griglia non strutturata -> strutturata
% np = 1E3;
% xq = linspace(min(x),max(x),np);
% yq = linspace(min(y),max(y),np);
% uxy = tri2grid(msh.POS(:,1:2)',msh.TRIANGLES(:,1:3)',out.field.phi,xq,yq);
% [Ex,Ey] = gradient(-uxy,length(xq)/np);
% sc = 10; % choose any value suitable for required arrows
% hold on
% quiver(xq(1:sc:end,1:sc:end), yq(1:sc:end,1:sc:end), Ex(1:sc:end,1:sc:end), Ey(1:sc:end,1:sc:end),'k');
% hold off

np = 1E3;
xq = linspace(min(x),max(x),np);
yq = linspace(min(y),max(y),np);
uxy = tri2grid(msh.POS(:,1:2)',msh.TRIANGLES(:,1:3)',out.field.phi,xq,yq);
[Ex,Ey] = gradient(-uxy,length(xq)/np);

radius = 1E-1;
theta = 0:30:360;
startX = radius * cosd(theta);
startY = radius * sind(theta);
XY = stream2(xq,yq,Ex,Ey,startX,startY);
streamline(XY); % streamslice(xq,yq,Ex,Ey)

%% POISSON 3D
close all
figure
trisurf(msh.TRIANGLES(:,1:3),x,y,out.field.phi,out.field.phi,edgecolor='none')
alpha 0.8
xlabel('x (m)'), ylabel('y (m)'); zlabel('solution - \phi'), axis tight;
ax = gca; ax.FontSize = 12;
f = gcf; colormap(f,ap.map.red_white_blue); exportgraphics(f,strcat(figspath,'poisson_solution.png'),'Resolution',300);
% shading interp; % colormap(f,winter(20)); %view(2);
%% ERRORE 2D
figure
trisurf(msh.TRIANGLES(:,1:3),x,y,z,(out.field.phi-sol_e)) % edgecolor='none'
view(2); axis equal; colorbar; xlabel('x (m)'), ylabel('y (m)'); axis tight;
ax = gca; ax.FontSize = 12;
f = gcf; exportgraphics(f,strcat(figspath,'poisson_error.png'),'Resolution',300);
%% ERRORE 3D
figure
trisurf(msh.TRIANGLES(:,1:3),x,y,(out.field.phi-sol_e),(out.field.phi-sol_e),edgecolor='none')
xlabel('x (m)'), ylabel('y (m)'); zlabel('error'); axis tight;
ax = gca; ax.FontSize = 12;
f = gcf; colormap(f,ap.map.red_white_blue); exportgraphics(f,strcat(figspath,'poisson_error.png'),'Resolution',300);

%% VALIDATION OF TIME-DOMAIN SOLVER
% OneNote -> VALIDAZIONE CODICI FEM
% manca possibilità di definire RHS in base a coordinate oltre che per
% regione

%% ESEMPIO LABClass
utils_FEM;
% MESH
LABClass; ndom = num_regions(msh);
% BC (Dirichlet)
opts.tag_boundary = [5, 6, 7, 8];
% MATERIALS
[opts.materials] = set_materials('labClass',ndom);
% PROBLEM KIND
opts.ProblemKind = 'QMagnetostaticSin'; % [Electrostatic][Magnetostatic][QMagnetostaticSin][MagTimeDependent]
opts.source = [0, 1e6, -1e6, 0];
opts.freq = 50;
% DIAGNOSTICA
opts.flag.print_measured_time = 0; 

% solution
[out] = fesolve(msh,opts);

disp(table((1:4)',abs(out.scal.I),'VariableNames',{'Conductor','I (A)'}));

%% VERIFICA FEM TEMPODIPENDENTE
utils_FEM;
% MESH
LABClass; ndom = num_regions(msh);
% BC (Dirichlet)
opts.tag_boundary = [5, 6, 7, 8];
% MATERIALS
[opts.materials] = set_materials('labClass',ndom);

% PROBLEM KIND
opts.ProblemKind = 'MagTimeDependent'; % [Electrostatic][Magnetostatic][QMagnetostaticSin][MagTimeDependent]
opts.source = @(time) [0, 1E6*sin(2*pi*50*time), -1E6*sin(2*pi*50*time), 0];
% TIME
time_fin = 0.1;
opts.time_array = [linspace(0,time_fin,400)]; % [0, logspace(-8,-3,20), linspace(2E-3,time_fin,100)];
% SAVE
opts.sv.nit_skip = 2; % variables saved each <nit_skip> time iterations
% DECOMPOSITION
opts.flag.decomp = 1; % 0: constant dt; 1: variable dt

% solution
[out] = fesolve(msh,opts); % 15.78 s, speedup ~200x

plot(out.sv.time,out.scal.I(2,:)); grid on; legend('I_{pipe}','Location','best');

%% PIPE SINUSOIDALE stazionario
utils_FEM;
% MESH
mesh_LAPLACE; ndom = num_regions(msh);
% BC (Dirichlet)
opts.tag_boundary = 70;
% MATERIALS
[opts.materials] = set_materials('Laplace',ndom);

% PROBLEM KIND
opts.ProblemKind = 'QMagnetostaticSin'; % [Electrostatic][Magnetostatic][QMagnetostaticSin][MagTimeDependent]
opts.source = [1*3.5d7/8.6165, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; % /8.6160 per avere I_phase1 = 250 A
opts.freq = 50;
% DIAGNOSTICA
opts.flag.print_measured_time = 0; 

% solution
[out] = fesolve(msh,opts);

disp(table((1:7)',opts.materials(1:7),abs(out.scal.I(1:7)),'VariableNames',{'Conductor','Material','I (A)'}));

% F90:
% f=50 Hz
%  CORRENTI NEI CONDUTTORI DI LINEA [MODULO E FASE]
%            1   2154.10690307392      A  -84.0146246459637      gradi
%            2   421.539540225472      A   98.0396871398149      gradi
%            3   607.742500436835      A   99.2918065158481      gradi
%  CORRENTE NEL TERRENO [MODULO E FASE]
%    80.2298462568369      A   174.111118772464      gradi
%  CORRENTE NELLA PIPELINE [MODULO E FASE]
%    638.235153598192      A   93.2347168411173      gradi

% MATLAB:
%     Conductor     Material      I (A) 
%     _________    ___________    ______
% 
%         1        "aluminium"    2154.1
%         2        "aluminium"    421.53
%         3        "aluminium"    607.69
%         4        "air"               0
%         5        "air"               0
%         6        "soil"         80.152
%         7        "aMagSteel"    638.51
%---------------------------------------------------------------------

% F90:
% f=50*1E5 Hz
%  CORRENTI NEI CONDUTTORI DI LINEA [MODULO E FASE]
%            1  2.364737239225055E-002 A  -89.7143274847057      gradi
%            2  3.346534909834859E-003 A   88.4607081052082      gradi
%            3  5.487951290631285E-003 A   89.1908612629313      gradi
%  CORRENTE NEL TERRENO [MODULO E FASE]
%   1.329687899293348E-002 A   94.9826763778011      gradi
%  CORRENTE NELLA PIPELINE [MODULO E FASE]
%   1.482018204292748E-003 A   55.1496600105430      gradi

% MATLAB:
%     Conductor     Material        I (A)  
%     _________    ___________    _________
% 
%         1        "aluminium"     0.023651
%         2        "aluminium"    0.0033471
%         3        "aluminium"    0.0054889
%         4        "air"                  0
%         5        "air"                  0
%         6        "soil"          0.013686
%         7        "aMagSteel"     0.001482
%---------------------------------------------------------------------

% F90
% f=50*1E8 Hz
%  CORRENTI NEI CONDUTTORI DI LINEA [MODULO E FASE]
%            1  2.380672608461276E-005 A  -89.9873039853960      gradi
%            2  3.227674715979254E-006 A   89.9302949354160      gradi
%            3  5.386608377051516E-006 A   89.9641477027882      gradi
%  CORRENTE NEL TERRENO [MODULO E FASE]
%   1.486473653394252E-005 A   90.0512782386711      gradi
%  CORRENTE NELLA PIPELINE [MODULO E FASE]
%   8.774831218624900E-011 A  -61.4231784532681      gradi

% MATLAB:
%     Conductor     Material        I (A)   
%     _________    ___________    __________
% 
%         1        "aluminium"    2.3823e-05
%         2        "aluminium"    3.2298e-06
%         3        "aluminium"    5.3902e-06
%         4        "air"                   0
%         5        "air"                   0
%         6        "soil"         2.2523e-05
%         7        "aMagSteel"    8.7779e-11
%---------------------------------------------------------------------

% F90:
% f=50*1E10 Hz
%  CORRENTI NEI CONDUTTORI DI LINEA [MODULO E FASE]
%            1  2.381066012605793E-007 A  -89.9983689332428      gradi
%            2  3.224827163752478E-008 A   89.9972338746277      gradi
%            3  5.384201063941230E-008 A   89.9986185014418      gradi
%  CORRENTE NEL TERRENO [MODULO E FASE]
%   1.487399624756786E-007 A   90.0020273566644      gradi
%  CORRENTE NELLA PIPELINE [MODULO E FASE]
%   5.151965053070691E-025 A  -79.5386307977180      gradi

% MATLAB:
%     Conductor     Material        I (A)   
%     _________    ___________    __________
% 
%         1        "aluminium"    2.3827e-07
%         2        "aluminium"     3.227e-08
%         3        "aluminium"    5.3878e-08
%         4        "air"                   0
%         5        "air"                   0
%         6        "soil"          2.875e-07
%         7        "aMagSteel"     5.152e-25

%% Plot densità Az e linee di campo per B
close all
x = msh.POS(:,1); y = msh.POS(:,2); z = msh.POS(:,3);

t_plot = 0; w = 2*pi*opts.freq; % 0.2+10E-3
vec = out.field.Az; 
Azt = abs(vec).*sin(w*t_plot+angle(vec));
trisurf(msh.TRIANGLES(:,1:3),x,y,z,Azt,edgecolor='none');
xplot = [-5 6]; yplot = [10 18]; % xplot = [10 100]; yplot = [-100 100];
xlim(xplot); ylim(yplot);
view(2); colorbar;
xlabel('x (m)'); ylabel('y (m)');
ax = gca; ax.FontSize = 12;

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
exportgraphics(f,strcat(figspath,'qstaz_50Hz_field_phase.png'),'Resolution',300);

figure
I = out.scal.I(1);
Iph1t = @(t) abs(I)*sin(w*t + angle(I));
I = out.scal.I(7);
Ip = @(t) abs(I)*sin(w*t + angle(I));
time = linspace(0.2,0.25,100);
plot(time,Iph1t(time(:)),time,Ip(time),'LineWidth',1)
xlim([0.2 0.25]); ylim([-300 300]); grid on; legend('i_{phase 1}(t)','i_p(t)')

%% PIPE SINUSOIDALE tempodipendente
utils_FEM;
% MESH
mesh_LAPLACE; ndom = num_regions(msh);
% BC (Dirichlet)
opts.tag_boundary = 70;
% MATERIALS
[opts.materials] = set_materials('Laplace',ndom);
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
tic
[out] = fesolve(msh,opts);
toc

% SCALARS
%%
load TimeMarch_sinusoidale.mat
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
ax = gca; ax.FontSize = 12;
f = gcf; exportgraphics(f,strcat(figspath,'sinusoidal_t.png'),'Resolution',300);

%% CALCOLO AMPIEZZE
npp = 50;
iph1 = out.scal.I(1,end-npp:end)./8.6062;
A_iph1 = (max(iph1)-min(iph1))/2
iph2 = out.scal.I(2,end-npp:end)./8.6062;
A_iph2 = (max(iph2)-min(iph2))/2
iph3 = out.scal.I(3,end-npp:end)./8.6062;
A_iph3 = (max(iph3)-min(iph3))/2
ip = out.scal.I(7,end-npp:end)./8.6062;
A_ip = (max(ip)-min(ip))/2

% A_iph1 =
%   250.0006
% A_iph2 =
%    48.9731
% A_iph3 =
%    70.6019
% A_ip =
%    74.0930

%% PIPE lightning
utils_FEM;
load Piantini_digitized.csv
% applico fattori di scala descritti nel paper
V_piantini(:,1) = Piantini_digitized(:,1)*50*1e-9; % scalo time di 50 e converto da ns a s
V_piantini(:,2) = Piantini_digitized(:,2)*18000 ; % scalo TENSIONE di 18000
[pp] = spline(V_piantini(:,1),V_piantini(:,2));
V_spline = ppval(pp, V_piantini(:,1));
plot(V_piantini(:,1),V_piantini(:,2),'o-',V_piantini(:,1),V_spline,'--');
legend('digitized','spline')

% MESH
mesh_LAPLACE; ndom = num_regions(msh);
% BC (Dirichlet)
opts.tag_boundary = 70;
% MATERIALS
[opts.materials] = set_materials('Laplace',ndom);
% PROBLEM KIND
opts.ProblemKind = 'MagTimeDependent'; % [Electrostatic][Magnetostatic][QMagnetostaticSin][MagTimeDependent]
% TIME
opts.time_array = [linspace(0,1E-3,1E5)]; % [0, logspace(-8,-3,20), linspace(2E-3,time_fin,100)];
opts.source = @(time) lightning(time,pp);
% SAVE
opts.sv.nit_skip = 10; % variables saved each <nit_skip> time iterations
% DECOMPOSITION
opts.flag.decomp = 1; % 0: constant dt; 1: variable dt
% DT
opts.flag.dt_auto = 0;

% solution
[out] = fesolve(msh,opts);

semilogx(out.sv.time,out.scal.I(7,:)); grid on; legend('I_{pipe}','Location','best'); xlabel('t (s)'); ylabel('I_p (A)');

out.sv.phi = []; out.field = [];
save('TimeMarch_sinusoidale.mat','out')

%% PIPE lightning dt automatico
utils_FEM;
load Piantini_digitized.csv
% applico fattori di scala descritti nel paper
V_piantini(:,1) = Piantini_digitized(:,1)*50*1e-9; % scalo time di 50 e converto da ns a s
V_piantini(:,2) = Piantini_digitized(:,2)*18000 ; % scalo TENSIONE di 18000
[pp] = spline(V_piantini(:,1),V_piantini(:,2));
V_spline = ppval(pp, V_piantini(:,1));
plot(V_piantini(:,1),V_piantini(:,2),'o-',V_piantini(:,1),V_spline,'--');
legend('digitized','spline')

% MESH
mesh_LAPLACE_v3; ndom = num_regions(msh);
% BC (Dirichlet)
opts.tag_boundary = (71:74);
% MATERIALS
[opts.materials] = set_materials('Laplace',ndom);
% PROBLEM KIND
opts.ProblemKind = 'MagTimeDependent'; % [Electrostatic][Magnetostatic][QMagnetostaticSin][MagTimeDependent]
% TIME
opts.source = @(time) lightning(time,pp);
% SAVE
opts.sv.nit_skip = 1; % variables saved each <nit_skip> time iterations
% DECOMPOSITION
opts.flag.decomp = 0; % 0: constant dt; 1: variable dt
% DT
opts.flag.dt_auto = 1;

% solution
[out] = fesolve(msh,opts);

semilogx(out.sv.time,out.scal.I(7,:)); grid on; legend('I_{pipe}','Location','best'); xlabel('t (s)'); ylabel('I_p (A)');

out.field = [];
out.sv.phi = [];
save('fulm_FEM_LT_v3','out')

%% PIPE lightning dt automatico Vlightning
utils_FEM;
load Vlightning.mat

% MESH
mesh_LAPLACE_v3; ndom = num_regions(msh);
% BC (Dirichlet)
opts.tag_boundary = (71:74);
% MATERIALS
[opts.materials] = set_materials('Laplace',ndom);
% PROBLEM KIND
opts.ProblemKind = 'MagTimeDependent'; % [Electrostatic][Magnetostatic][QMagnetostaticSin][MagTimeDependent]
% TIME
time = 0;
opts.source = @(time) lightning_interp(time,Vlightning);
% SAVE
opts.sv.nit_skip = 1; % variables saved each <nit_skip> time iterations
% DECOMPOSITION
opts.flag.decomp = 0; % 0: constant dt; 1: variable dt
% DT
opts.flag.dt_auto = 1;

opts.t_end = 0.1;
opts.absvar = 0.5;
opts.relvar = 0.1;
opts.dt_first = 2;

% solution
[out] = fesolve(msh,opts);

semilogx(out.sv.time,out.scal.I(7,:)); grid on; legend('I_{pipe}','Location','best'); xlabel('t (s)'); ylabel('I_p (A)');

out.field = [];
out.sv.phi = [];
save('fulm_FEM_LT_v3_interp','out')

%% PIPE step dt automatico
utils_FEM;

% MESH
mesh_LAPLACE_v3; ndom = num_regions(msh);
% BC (Dirichlet)
opts.tag_boundary = (71:74);
% MATERIALS
[opts.materials] = set_materials('Laplace',ndom);
% PROBLEM KIND
opts.ProblemKind = 'MagTimeDependent'; % [Electrostatic][Magnetostatic][QMagnetostaticSin][MagTimeDependent]
% TIME
opts.source = @(time) step(time);
% SAVE
opts.sv.nit_skip = 1; % variables saved each <nit_skip> time iterations
% DECOMPOSITION
opts.flag.decomp = 0; % 0: constant dt; 1: variable dt
% DT
opts.flag.dt_auto = 1;

opts.t_end = 0.1;
opts.absvar = 0.5;
opts.relvar = 0.1;
opts.dt_first = 2;

% solution
[out] = fesolve(msh,opts);

semilogx(out.sv.time,out.scal.I(7,:)); grid on; legend('I_{pipe}','Location','best'); xlabel('t (s)'); ylabel('I_p (A)');

out.field = [];
out.sv.phi = [];
save('step_FEM_LT','out')

%% Confronto risposte al gradino (t) e (L)
addpath("..\F90_mesh_LAPL\Qstaz_Lapl_Inv_F_ramp")
load results_lapl.csv

load step_FEM_LT.mat;

semilogx(out.sv.time,out.scal.I(7,:),results_lapl(:,1),results_lapl(:,2),'o')
legend('time-dependent','inverse laplace','Location','best'); grid on; xlabel('t (s)'); ylabel('I_p (A)'); xlim([1E-7 inf]);

%% PIPE sinus dt automatico
utils_FEM;

% MESH
mesh_LAPLACE_v3; ndom = num_regions(msh);
% BC (Dirichlet)
opts.tag_boundary = (71:74);
% MATERIALS
[opts.materials] = set_materials('Laplace',ndom);
% PROBLEM KIND
opts.ProblemKind = 'MagTimeDependent'; % [Electrostatic][Magnetostatic][QMagnetostaticSin][MagTimeDependent]
% TIME
w = 2*pi*50;
opts.source = @(time) sine(time,w);
% SAVE
opts.sv.nit_skip = 1; % variables saved each <nit_skip> time iterations
% DECOMPOSITION
opts.flag.decomp = 0; % 0: constant dt; 1: variable dt
% DT
opts.flag.dt_auto = 1;

opts.t_end = 0.1;
opts.absvar = 0.5;
opts.relvar = 0.1;
opts.dt_first = 1E-9;

% solution
[out] = fesolve(msh,opts);

semilogx(out.sv.time,out.scal.I(7,:)); grid on; legend('I_{pipe}','Location','best'); xlabel('t (s)'); ylabel('I_p (A)');

out.field = [];
out.sv.phi = [];
save('sine_FEM_LT','out')

%% PIPE LAPLACE (validazione con LAPLACE F90)
utils_FEM;
% MESH
mesh_LAPLACE; ndom = num_regions(msh);
% BC (Dirichlet)
opts.tag_boundary = 70;
% MATERIALS
[opts.materials] = set_materials('Laplace',ndom);
% PROBLEM KIND
opts.ProblemKind = 'QMagnetostaticSin_LAPL'; % [Electrostatic][Magnetostatic][QMagnetostaticSin][MagTimeDependent]
opts.s_lapl = complex(12202590.5306814,1532484.22126331); 
omega = 2*pi*50;
opts.source = @(s_lapl) laplace_sine(s_lapl,omega); % @(s) laplace_impulse(); % @(s_lapl) laplace_sine(s_lapl,omega);
% DIAGNOSTICA
opts.flag.print_measured_time = 0; 

% solution
opts.source = opts.source(opts.s_lapl);
[out] = fesolve(msh,opts);

ii_t = [1,6,7];
disp(table((ii_t)',opts.materials(ii_t),abs(out.scal.I(ii_t)),'VariableNames',{'Conductor','Material','I (A)'}));

% C:\Users\arturo.popoli2\Documents\GitHub\Papers\2022 Applied Sciences Laplace\dev\Qstaz_Lapl_Inv_F
% s = complex(101.688254422345,0)
% I_ph1     14.7490102783044
% I_terr    0.233358207903236
% I_p       4.73789499993623

% MATLAB
%     Conductor     Material       I (A) 
%     _________    ___________    _______
% 
%         1        "aluminium"     14.824
%         6        "soil"         0.20253
%         7        "aMagSteel"     4.7229
%---------------------------------------------------------------------

% C:\Users\arturo.popoli2\Documents\GitHub\Papers\2022 Applied Sciences Laplace\dev\Qstaz_Lapl_Inv_F
% s = complex(12202590.5306814,0.000000000000000E+000)
%  I_ph1  1.278498799849864E-013
%  I_terr  6.685293046938621E-014
%  I_p  1.009080090543553E-014

% MATLAB
%     Conductor     Material        I (A)   
%     _________    ___________    __________
% 
%         1        "aluminium"    1.2787e-13
%         6        "soil"         6.8502e-14
%         7        "aMagSteel"    1.0091e-14
%---------------------------------------------------------------------

% C:\Users\arturo.popoli2\Documents\GitHub\Papers\2022 Applied Sciences Laplace\dev\Qstaz_Lapl_Inv_F
% s = complex(12202590.5306814,1532484.22126331)
%  I_ph1  1.248887409979681E-013
%  I_terr  6.534584710814666E-014
%  I_p  9.837143451562327E-015

% MATLAB
%     Conductor     Material        I (A)   
%     _________    ___________    __________
% 
%         1        "aluminium"    1.2491e-13
%         6        "soil"         6.6961e-14
%         7        "aMagSteel"    9.8371e-15
%---------------------------------------------------------------------

%% CONFRONTO MATRICI
addpath("..\F90_mesh_LAPL\Qstaz_Lapl_Inv_F_impulse\diagn")
[out_f90]=loadFEM(1);

utils_FEM;
% MESH
mesh_LAPLACE_v3; ndom = num_regions(msh);
% BC (Dirichlet)
opts.tag_boundary = 71:74;
% MATERIALS
[opts.materials] = set_materials('Laplace',ndom);
% PROBLEM KIND
opts.ProblemKind = 'QMagnetostaticSin_LAPL'; % [Electrostatic][Magnetostatic][QMagnetostaticSin][MagTimeDependent]
% DIAGNOSTICA
opts.flag.print_measured_time = 0; 
opts.flag.diagnostics_linear_system = 1;

opts.s_lapl = out_f90.s;
opts.source = @(s_lapl) laplace_impulse(s_lapl); opts.source = opts.source(opts.s_lapl);
[out] = fesolve(msh,opts);

ii_t = [1,6,7];
disp(table((ii_t)',opts.materials(ii_t),abs(out.scal.I(ii_t)),abs(out_f90.I),'VariableNames',{'Conductor','Material','I (A)','I_FEM (A)'}));

% compar_mat(out_f90.sol,out.field.Az);

f = @() fesolve(msh,opts);
fprintf('FINESSE_time = %f s for %i nodes\n',timeit(f),size(msh.POS,1));

%% PIPE LAPLACE VERIFICA IMPULSO
utils_FEM;
% MESH
mesh_LAPLACE; ndom = num_regions(msh);
% BC (Dirichlet)
opts.tag_boundary = 70;
% MATERIALS
[opts.materials] = set_materials('Laplace',ndom);
% PROBLEM KIND
opts.ProblemKind = 'QMagnetostaticSin_LAPL'; % [Electrostatic][Magnetostatic][QMagnetostaticSin][MagTimeDependent]
opts.s_lapl = complex(1220259053068.14,0.000000000000000D+000);
opts.source = @(s_lapl) laplace_impulse(s_lapl);
% DIAGNOSTICA
opts.flag.print_measured_time = 0;

% solution
opts.source = opts.source(opts.s_lapl);
tic
[out] = fesolve(msh,opts);
toc
ii_t = [1,6,7];
disp(table((ii_t)',opts.materials(ii_t),abs(out.scal.I(ii_t)),'VariableNames',{'Conductor','Material','I (A)'}));

%  s_lapl (1220259053068.14,0.000000000000000E+000)
%  I_ph1  6.130101159959831E-007
%  I_terr  3.829196602735363E-007 @@@
%  I_p  2.138318947087393E-019

%     Conductor     Material        I (A)   
%     _________    ___________    __________
% 
%         1        "aluminium"    6.1341e-07
%         6        "soil"         7.4904e-07 @@@
%         7        "aMagSteel"    2.1384e-19

%% Pipe D'Amore impulso
load("..\Convolution\pipe_impulso_mur_1.csv");
t_imp_f90 = pipe_impulso_mur_1(:,1);
I_p_imp_f90 = pipe_impulso_mur_1(:,2);

addpath("..\Convolution\796_matlab\")
help lapinv_796

% HANDLE FUNCTION source PER FEM(s)
fun_sorg = @(s) laplace_impulse(s);
% NB se avessi bisogno di s (fornito da D'Amore) ed omega (costante) come input per il termine di source:
% omega = 50;
% fun_sorg = @(s) sorgente_sine(s,omega);

% HANDLE FUNCTION FEM SOLVER PER lapinv796
fun_FEM = @(s) ilaplace(s,fun_sorg);

tarray = t_imp_f90(1:2)'; % linspace(0,4*pi,250);
tol = 10^-5; ssbar = 0; nmax = 50; sigma0 = 0.; % nmax = 550;
[tarray,fzinv,~,~] = lapinv_796(fun_FEM,tarray,tol,ssbar,nmax,sigma0);

%% Pipe D'Amore sinusoidale
addpath("..\Convolution\796_matlab\")
help lapinv_796

% HANDLE FUNCTION source PER FEM(s)
omega = 2*pi*50;
fun_sorg = @(s) laplace_sine(s,omega);

% HANDLE FUNCTION FEM SOLVER PER lapinv796
fun_FEM = @(s) ilaplace(s,fun_sorg);

tarray = linspace(1E-6,0.25,50); % % linspace(0,4*pi,250);
tol = 1E-6; ssbar = 0; nmax = 550; sigma0 = 0.; % nmax = 550; tol = 10^-5;
[tarray,fzinv,~,~] = lapinv_796(fun_FEM,tarray,tol,ssbar,nmax,sigma0);

% save('InvLap_sinusoidale.mat','tarray','fzinv','omega','tol','nmax','ssbar','sigma0','fun_FEM','fun_sorg')

%% Elab
close all

load InvLap_sinusoidale.mat
ip = fzinv/8.6160; % 8.6062; 
A_ip_LapInv = (max(ip)-min(ip))/2;

load TimeMarch_sinusoidale.mat
plot(out.sv.time,out.scal.I(7,:)/8.6160,'b:',tarray,ip,'x','LineWidth',1.2); 
grid on; lgd=legend('I_{pipe}(t) - Time Marching','I_{pipe}(t) - Inverse Laplace','Location','southeast');
xlabel('t (s)');
ylabel('Current (A)');
ax = gca; ax.FontSize = 12;
f = gcf; exportgraphics(f,strcat(figspath,'TimeMarch_vs_InvLap.png'),'Resolution',300);
