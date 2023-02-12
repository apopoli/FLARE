function [out] = FEM2D00_t(ProblemKind, time_array, p, t, edgeBound, bcflag_e, bcval_e, ireg, iregbe, materials, source)
%FEM2D Finite element solution of div(prop(x,y) grad(phi)) = source(x,y)
% material(id_phys_reg) = MaterialKind assigned to the region id_phys_reg

ndom = length(unique(ireg)); % number of regions within the domain

eps0 = 8.854187817e-12;
mu0 = 4*pi*1.e-7;
np = size(p,1); % number of mesh points
nt = size(t,1); % number of triangles
nbedges = size(edgeBound,1); % number (element) edges on the boundary

prop_e = zeros(nt,1);
sigma = zeros(nt,1);
src = zeros(nt,1);

BCflag_p = zeros(np,1); % init. flags for BC type on nodes
BCval_p = zeros(np,1); % values for BC on nodes

EdgeNeumann = find(bcflag_e == 2);
dp = p(edgeBound(EdgeNeumann,1),:) - p(edgeBound(EdgeNeumann,2),:);
L = (dp(:,1).^2 + dp(:,2).^2).^0.5;

% homogeneous Dirichlet BCs, regardless of input
ipc = unique(edgeBound);

BCflag_p(ipc) = 2;
BCval_p(ipc) = 0.0;
RHSg = zeros(np,1); % init RHS array

% cell array with regions
el_ireg = cell(ndom,1); % element indexes for each region
tic
for i = 1:ndom
    PROPel = MatLib(materials(i)); % region properties
    el_ireg{i} = find(ireg==i); % elements in region
    prop_e(el_ireg{i}) = 1/PROPel(2);
    sigma(el_ireg{i}) = PROPel(3);
end
fprintf('assembly of cell array for regions  = %g s \n',toc);

tic
[K, S, Area,gradN] = assembling(p,t,prop_e(:),sigma(:));
fprintf('global matrices assembly - vectorized = %g s \n',toc);

% Vectorized version
% aux vectors for sparse assembly
ii_sp = repelem(1:nt,3); % 1 1 1 2 2 2 3 3 3 ...
jj_sp = reshape(t', [], 1); % unrolled rows of t -> vertexes of all triangles 
tic
K_rhs = sparse(jj_sp,ii_sp,1,np,nt);
fprintf('assembly of aux matrix for RHS = %g s \n',toc);

% Indexes for BCs
ii_BCflag_p_2 = find(BCflag_p == 2); % Dirichlet nodes indices (corresponds to ipc)
ind_lin_p_2 = sub2ind([np np],ii_BCflag_p_2,ii_BCflag_p_2); % linear indexes of Dirichlet nodes, must be forced to 1 in Kg

% SETTINGS TIME LOOP
disp('time  loop')
mode_assembly_rhs = 2; % 1=for loop; 2=matXvec
mode_apply_BC = 2; % 1=for loop; 2=vectorized

dt = diff(time_array);
phi = zeros(np,1); % TODO: phi(t=0) given by initial conditions

% SETTINGS SAVE
sv.nit_skip = 1; % variables saved each <nit_skip> time iterations
sv.nit_save = length(time_array); % numero salvataggi
sv.phi = zeros(np,sv.nit_save); % saved potential
sv.phi(:,1) = phi; % initial value of phi 
sv.time = zeros(sv.nit_save,1);
sv.ii_skip = 0; % initialization

% init time-intervals for profiling
tm.assembly_for = 0; tm.assembly_matXvec = 0; tm.BC = 0; tm.system_solve = 0;

sv.phi(:,1) = phi; % CN yields phi_{k+1}, so phi_{1} is the starting solution

for i_it=1:length(time_array)-1
    
    fprintf('it %i / %i\n',i_it,length(time_array)-1);

    % matrices for CN
    M1 = 0.5*K + 1/dt(i_it)*S;
    M2 = -0.5*K + 1/dt(i_it)*S;

    sorg = source((time_array(i_it)+time_array(i_it+1))/2); % evaluate source terms at time k+1/2
    switch mode_assembly_rhs
        case (1)
            tic
            RHSg = zeros(np,1);
            for iel=1:nt
                src(iel) = sorg(ireg(iel)); % assign the source to given element depending on the region 
                PV = [ones(3,1), p(t(iel,:),:)];
                AreaDoppia = det(PV);
                Area = 0.5 * AreaDoppia;
                RHSel = mu0 * src(iel) * Area / 3 * ones(3,1);
                % assembly di RHSel in RHSg
                RHSg(t(iel,:)) = RHSg(t(iel,:)) + RHSel;
            end
            tm.assembly_for = tm.assembly_for + toc;
        case (2)
            tic
            src = sorg(ireg(:))';
            prop_el = 1/3*mu0*src.*Area;
            RHSg =  K_rhs * prop_el;
            tm.assembly_matXvec = tm.assembly_matXvec + toc;
    end
    
    tic
    RHSg = RHSg + M2*phi;

    Kg = M1;

    switch mode_apply_BC % BCs
        case 1 % for loop
        for ip = 1:np
            switch BCflag_p(ip)
                case(1) % Neumann BC
                    % TODO
                case(2) % Dirichlet 
                    Kg(ip,:) = 0.0;
                    Kg(ip, ip) = 1.;
                    RHSg(ip) = BCval_p(ip);
            end
        end
        case 2 % vectorized
            Kg(ii_BCflag_p_2,:) = 0;
            Kg(ind_lin_p_2) = 1;
            RHSg(ii_BCflag_p_2) = BCval_p(ii_BCflag_p_2);
    end
    tm.BC = tm.BC + toc;

    tic
    phi = Kg\RHSg;
    tm.system_solve = tm.system_solve + toc;

    if mod(i_it,sv.nit_skip)==0
        sv.ii_skip = sv.ii_skip + 1;
        sv.time(sv.ii_skip+1) = time_array(i_it+1);
        sv.phi(:,sv.ii_skip+1) = phi; % +1 since column 1 is filled with initial condition, CN yields variables at next time-step
    end
end

disp('Starting post-processing...')

%% POST PROCESSING
phi = sv.phi;

% Mat_grdn_x e Mat_grdn_y using SPARSE
gradN_x = squeeze(gradN(:,1,:));
gradN_y = squeeze(gradN(:,2,:));

tic
v_sp_x = reshape(gradN_x', [], 1);
v_sp_y = reshape(gradN_y', [], 1);
fprintf('collecting indexes for sparse = %g s \n',toc);

tic
Mat_grdn_x = sparse(ii_sp,jj_sp,v_sp_x,nt,np);
Mat_grdn_y = sparse(ii_sp,jj_sp,v_sp_y,nt,np);
fprintf('assembly sparse = %g s \n',toc);

% cell array with regions
el_ireg = cell(ndom,1); % element indexes for each region
tic
for i = 1:ndom
    el_ireg{i} = find(ireg==i);
end
fprintf('assembly of cell array for regions  = %g s \n',toc);

% init time-intervals for profiling
tm.grad = 0; tm.switch = 0; tm.int = 0; tm.node = 0; tm.save = 0;

% vector fields 
field1 = 'Az';  value1 = zeros(np,sv.nit_save);
field2 = 'Bx';  value2 = zeros(np,sv.nit_save);
field3 = 'By';  value3 = zeros(np,sv.nit_save);
field4 = 'Hx';  value4 = zeros(np,sv.nit_save);
field5 = 'Hy';  value5 = zeros(np,sv.nit_save);
field = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5);

% scalar fields 
field1 = 'I';  value1 = zeros(ndom,sv.nit_save);
scal = struct(field1,value1);

IntSource_AP = zeros(ndom,1);
sorg = source(sv.time(1)); % evaluate source terms at present time
src = sorg(ireg(:))'; % assign the source to given element depending on the region 
for i = 1:ndom
    IntSource_AP(i) = sum(src(el_ireg{i}).*Area(el_ireg{i})); % integral source J0
    % assumo dA/dt = 0 all'istante iniziale (sarebbe da sommare...)
end
scal.I(:,1) = IntSource_AP; % current in each region

for i_t = 1:sv.nit_save % time  loop
    tic
    switch i_t
        case (1) % first instant
            delta_phi = -3*phi(:,i_t)+4*phi(:,i_t+1)-1*phi(:,i_t+2);
            delta_t = sv.time(2)-sv.time(1) + sv.time(3)-sv.time(2);
        case (sv.nit_save) % last instant
            delta_phi = 3*phi(:,i_t)-4*phi(:,i_t-1)+1*phi(:,i_t-2);
            delta_t = sv.time(i_t)-sv.time(i_t-1) + sv.time(i_t-1)-sv.time(i_t-2);
        otherwise % central instants
            delta_phi = phi(:,i_t+1)-phi(:,i_t-1);
            delta_t = sv.time(i_t+1)-sv.time(i_t) + sv.time(i_t)-sv.time(i_t-1);
    end
    tm.switch = tm.switch + toc;

    sorg = source(sv.time(i_t)); % evaluate source terms at present time
    src = sorg(ireg(:))'; % assign the source to given element depending on the region 

    % integral of element J0
    IntSource_AP = zeros(ndom,1);
    tic
    for i = 1:ndom
        IntSource_AP(i) = sum(src(el_ireg{i}).*Area(el_ireg{i})); % integral source J0
        for j = 1:3
            dA_dt = - delta_phi(t(el_ireg{i},j))/delta_t;
            IntSource_AP(i) = IntSource_AP(i) + sum((sigma(el_ireg{i}).*dA_dt).* Area(el_ireg{i})/3); % integral source dA/dt
        end
    end
    tm.int = tm.int + toc;

    tic
    gradphi_e_x = Mat_grdn_x * phi(:,i_t);
    gradphi_e_y = Mat_grdn_y * phi(:,i_t);
    tm.grad = tm.grad + toc;

    % repeat element quantities three times, to use accumarray for nodal values
    tic
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
    dA_dt_p = delta_phi/delta_t;
    Jz_p = Jz_p - sigma_p.*dA_dt_p;
    
    Jz_p(ips) = Jz_p(ips)./icounts(ips);
    tm.node = tm.node + toc;
    
    tic
    switch ProblemKind
        case 'MagTimeDependent'
            field.Az(:,i_t) = phi(:,i_t); % magnetic vector potential Az
            field.Bx(:,i_t) = gradphi_y; % magnetic flux density Bx       
            field.By(:,i_t) = -gradphi_x; % magnetic flux density By
            field.Hx(:,i_t) = gradphi_prop_y/mu0; % magnetic field Hx
            field.Hy(:,i_t) = -gradphi_prop_x/mu0; % magnetic field Hy
            field.Jz(:,i_t) = Jz_p; % current density Jz
            scal.I(:,i_t) = IntSource_AP; % current in each region
    end
    tm.save = tm.save + toc;
end

out.tm = tm;
out.field = field;
out.scal = scal;
out.sv = sv;

end

