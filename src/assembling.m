function [K,S,area,gradL] = assembling(node,elem,prop,sigma)
%Computes the global matrices K and S
%   [K,S,area] = assembling(node,elem,prop,sigma)

mu0 = 4*pi*1.e-7;

% Marini - Metodi Numerici 
% Sadiku, Numerical Techniquesin Electromagnetics, 2ed
% Edition p.393 (pdf)
% Sadiku, Numerical Techniquesin Electromagnetics, 3ed
% Edition p.422

%%
N = size(node,1); NT = size(elem,1);
ii = zeros(9*NT,1); jj = zeros(9*NT,1); sK = zeros(9*NT,1); sS = zeros(9*NT,1); 
l(:,:,1) = node(elem(:,3),:)-node(elem(:,2),:); % l_1 = [x_3 - x_2, y_3 - y_2]
l(:,:,2) = node(elem(:,1),:)-node(elem(:,3),:); % l_2 = [x_1 - x_3, y_1 - y_3]
l(:,:,3) = node(elem(:,2),:)-node(elem(:,1),:); % l_3 = [x_2 - x_1, y_2 - y_1]
area = 0.5*abs(-l(:,1,3).*l(:,2,2)+l(:,2,3).*l(:,1,2)); % 1/2*(x_2-x_1)*(y_3-y_1)-(x_3-x_1)*(y_2-y_1)

index = 0;
for i = 1:3
    for j = 1:3
        ii(index+1:index+NT) = elem(:,i);
        jj(index+1:index+NT) = elem(:,j);
        
        sK(index+1:index+NT) = prop.*dot(l(:,:,i),l(:,:,j),2)./(4*area);
        
        % Sadiku, Numerical Techniquesin Electromagnetics, 3ed Edition p.426
        if i==j
            coeff = 1/6;
        else
            coeff = 1/12;
        end
        sS(index+1:index+NT) = sigma .* area * mu0 * coeff;

        index = index + NT;
    end
end
K = sparse(ii,jj,sK,N,N); % global K matrix assembly
S = sparse(ii,jj,sS,N,N); % assemblo matrice S globale

% gradiente funzioni di forma
gradL = fliplr(-l./(2*area)); % swap columns
gradL(:,2,:) = gradL(:,2,:)*(-1); % change sign to column 2

end