function [K,area,gradL] = assembling_steady_state(ProblemKind,node,elem,prop,sigma,j_omega)
%Computes the global matrices K and S
%   [K,S,area] = assembling(node,elem,prop,sigma)

mu0 = 4*pi*1.e-7;

% Sadiku, Numerical Techniquesin Electromagnetics, 3ed Edition 
% Stiffness matrix [K]: p.384
% Mass matrix [S]: p.426

%%
N = size(node,1); NT = size(elem,1);
ii = zeros(9*NT,1); jj = zeros(9*NT,1); sK = zeros(9*NT,1); sS = zeros(9*NT,1);

x1 = node(elem(:,1),1);
x2 = node(elem(:,2),1);
x3 = node(elem(:,3),1);
y1 = node(elem(:,1),2);
y2 = node(elem(:,2),2);
y3 = node(elem(:,3),2);

l(:,:,1) = [x3-x2 y2-y3];
l(:,:,2) = [x1-x3 y3-y1];
l(:,:,3) = [x2-x1 y1-y2];

area = 0.5*((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1)); % anti-clockwise node numbering

switch ProblemKind
    case {"QMagnetostaticSin","QMagnetostaticSin_LAPL"}
        switch_S = 1;
    otherwise
        switch_S = 0;
end

index = 0;
for i = 1:3
    for j = 1:3
        % indexes for sparse
        ii(index+1:index+NT) = elem(:,i);
        jj(index+1:index+NT) = elem(:,j);
        % values of [K] for sparse
        sK(index+1:index+NT) = prop.*dot(l(:,:,i),l(:,:,j),2)./(4*area);

        % Add contribution given by mass matrix [S]
        if switch_S==1
            if i==j
                coeff = 1/6;
            else
                coeff = 1/12;
            end
            sS(index+1:index+NT) =  sigma .* area * mu0 * coeff; % sigma
            sK(index+1:index+NT) = sK(index+1:index+NT) + j_omega*sS(index+1:index+NT); % Kel = Kel + 1i*omega*Sel;
        end

        index = index + NT;
    end
end
K = sparse(ii,jj,sK,N,N); % global matrix assembly

% shape function gradient
gradL = fliplr(l./(2*area)); % swap columns and change sign
% gradL(:,2,:) = gradL(:,2,:)*(-1); % change sign to column 2

end