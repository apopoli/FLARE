function [out] = FEM2D00(msh, BC, opts)
%FEM2D Finite element solution of div(prop(x,y) grad(phi)) = source(x,y)
% material(id_phys_reg) = MaterialKind assigned to the region id_phys_reg

p = msh.POS(:,1:2);
t = msh.TRIANGLES(:,1:3);
ireg = msh.TRIANGLES(:,4);
edgeBound = msh.LINES(:,1:2); % lati (di elemento) sul contorno
ndom = length(unique(ireg)); % number of regions within the domain

% BC for each node
BCval_p = set_BCs_on_nodes(msh, BC);

eps0 = 8.854187817e-12;
mu0 = 4*pi*1.e-7;
np = size(p,1); % number of mesh points
nt = size(t,1); % number of triangles

prop_e = zeros(nt,1);
sigma = zeros(nt,1);

% Length of Neumann edges
i_EdgeNeumann = find(BC.BCflag_e == 1);
dp = p(edgeBound(i_EdgeNeumann,1),:) - p(edgeBound(i_EdgeNeumann,2),:);
L = (dp(:,1).^2 + dp(:,2).^2).^0.5;

% Indexes for BCs
ii_BCflag_p_2 = find(BCval_p(:,2) == 2); % Dirichlet nodes indices (corresponds to ipc)
ind_lin_p_2 = sub2ind([np np],ii_BCflag_p_2,ii_BCflag_p_2); % linear indexes of Dirichlet nodes, must be forced to 1 in K

% cell array with regions
el_ireg = cell(ndom,1); % element indexes for each region
tic
for i = 1:ndom
    PROPel = MatLib(opts.materials(i)); % region properties
    el_ireg{i} = find(ireg==i); % elements in region

    switch opts.ProblemKind
    case {"Electrostatic"}
        prop_e(el_ireg{i}) = PROPel(1); % dielectric constant, from MatLib
        j_omega = 0;
    case {"Magnetostatic"}
        prop_e(el_ireg{i}) = 1/PROPel(2); % inverse of permeability, from MatLib
        j_omega = 0;
    case {"QMagnetostaticSin"}
        prop_e(el_ireg{i}) = 1/PROPel(2);
        sigma(el_ireg{i}) = PROPel(3);
        j_omega = 1i*2*pi*opts.freq;
    case {"QMagnetostaticSin_LAPL"}
        prop_e(el_ireg{i}) = 1/PROPel(2);
        sigma(el_ireg{i}) = PROPel(3);
        j_omega = opts.s_lapl;
    end
end
tm.assembl_cell_array_reg = toc;

tic
[K, Area,gradN] = assembling_steady_state(opts.ProblemKind,p,t,prop_e(:),sigma(:),j_omega);
tm.assembl_matrici_globali = toc;

% aux vectors for sparse assembly
ii_sp = repelem(1:nt,3); % 1 1 1 2 2 2 3 3 3 ...
jj_sp = reshape(t', [], 1); % rows of "tri" (vertexes of triangles)
tic
K_rhs = sparse(jj_sp,ii_sp,1,np,nt);
tm.assembly_matrice_aux_RHS = toc;

% assembly RHS
sorg = opts.source().';
src = sorg(ireg(:));  % element source term, based on element region
switch opts.ProblemKind
    case "Electrostatic"
        prop_el = 1/eps0;
    case "Magnetostatic"
        prop_el = mu0;
    case {"QMagnetostaticSin","QMagnetostaticSin_LAPL"}
        prop_el = mu0;
end
prop_el = prop_el * 1/3*src.*Area;
RHSg =  K_rhs * prop_el;

% Neumann BCs
for k = 1:length(i_EdgeNeumann)
    e  = i_EdgeNeumann(k);
    n1 = edgeBound(e,1);
    n2 = edgeBound(e,2);

    g  = BC.BCval_e(e);
    Le = L(k);

    RHSg(n1) = RHSg(n1) + g * Le / 2;
    RHSg(n2) = RHSg(n2) + g * Le / 2;
end

% Dirichlet BCs
K(ii_BCflag_p_2,:) = 0;
K(ind_lin_p_2) = 1;
RHSg(ii_BCflag_p_2) = BCval_p(ii_BCflag_p_2,3);

% solve
phi = K\RHSg;

% POST PROCESSING
% Mat_grdn_x e Mat_grdn_y using SPARSE
gradN_x = squeeze(gradN(:,1,:));
gradN_y = squeeze(gradN(:,2,:));

v_sp_x = reshape(gradN_x', [], 1);
v_sp_y = reshape(gradN_y', [], 1);

Mat_grdn_x = sparse(ii_sp,jj_sp,v_sp_x,nt,np);
Mat_grdn_y = sparse(ii_sp,jj_sp,v_sp_y,nt,np);

IntSource = zeros(ndom,1);
for i = 1:ndom
    IntSource(i) = dot(src(el_ireg{i}),Area(el_ireg{i})); % integral source J0
    for j = 1:3
        dA_dt = j_omega*phi(t(el_ireg{i}));
        IntSource(i) = IntSource(i) - dot(sigma(el_ireg{i}),dA_dt.*Area(el_ireg{i})/3); % integral source dA/dt sum((sigma(el_ireg{i}).*dA_dt).* Area(el_ireg{i})/3); 
    end
end

gradphi_e_x = Mat_grdn_x * phi;
gradphi_e_y = Mat_grdn_y * phi;

% repeat element quantities three times, for accumarray 
src_rep = repelem(src,3);
gradphi_e_x_rep = repelem(gradphi_e_x,3);
gradphi_e_y_rep = repelem(gradphi_e_y,3);
prop_e_rep = repelem(prop_e,3);
sigma_rep = repelem(sigma,3);

icount = accumarray(jj_sp,1); % # of elements within the support domain of each node
icounts = accumarray(jj_sp,double(sigma_rep~=0)); % # the ones with non-null sigma
ips = find(icounts ~= 0); % indexes of nodes belonging only to elements with sigma = 0

gradphi_x = accumarray(jj_sp,gradphi_e_x_rep)./icount;
gradphi_y = accumarray(jj_sp,gradphi_e_y_rep)./icount;

gradphi_prop_x = accumarray(jj_sp,gradphi_e_x_rep.*prop_e_rep)./icount;
gradphi_prop_y = accumarray(jj_sp,gradphi_e_y_rep.*prop_e_rep)./icount;

sigma_p = accumarray(jj_sp,sigma_rep);
Jz_p = accumarray(jj_sp,src_rep);
if (opts.ProblemKind == "QMagnetostaticSin" || opts.ProblemKind == "QMagnetostaticSin_LAPL")
    dA_dt_p = j_omega*phi;
    Jz_p = Jz_p - sigma_p.*dA_dt_p;
end
Jz_p(ips) = Jz_p(ips)./icounts(ips);

switch opts.ProblemKind
    case 'Electrostatic'
        field.phi = phi; % electric potential
        field.Ex = -gradphi_x; % electric field Ex
        field.Ey = -gradphi_y; % electric field Ey
        field.Dx = -gradphi_prop_x*eps0; % electric induction Dx
        field.Dy = -gradphi_prop_y*eps0; % electric induction Dy      
    case 'Magnetostatic'
        field.Bx = gradphi_y; % magnetic flux density Bx
        field.By = -gradphi_x; % magnetic flux density By
        field.Hx = gradphi_prop_y/mu0; % magnetic field Hx
        field.Hy = -gradphi_prop_x/mu0; % magnetic field Hy
        field.Jz = Jz_p; % current density Jz
    case {'QMagnetostaticSin','QMagnetostaticSin_LAPL'}
        field.Az = phi; % magnetic vector potential Az
        field.Bx = gradphi_y; % magnetic flux density Bx       
        field.By = -gradphi_x; % magnetic flux density By
        field.Hx = gradphi_prop_y/mu0; % magnetic field Hx
        field.Hy = -gradphi_prop_x/mu0; % magnetic field Hy
        field.Jz = Jz_p; % current density Jz
        scal.I = IntSource; % current in each region
end

% output
out.BCval_p = BCval_p;

out.field = field;
if (exist('scal','var'))
    out.scal = scal;
end

if opts.flag.print_measured_time == 1
    fprintf('assembly of cell array for regions  = %g s \n',tm.assembl_cell_array_reg);
    fprintf('global matrices assembly - vectorized = %g s \n',tm.assembl_matrici_globali);
    fprintf('assembly of aux matrix for RHS = %g s \n\n',tm.assembly_matrice_aux_RHS);
end

