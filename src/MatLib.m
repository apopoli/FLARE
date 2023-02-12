function [Mprop] = MatLib(MatKind)
% Material library
% prop(1) = material relative permittivity [Adim]
% prop(2) = material relative permeability [Adim]
% prop(3) = material electric conductivity [S/m]
switch MatKind
    case "air"
        Mprop = [1, 1, 0];
    case "iron"
        Mprop = [1, 2000, 1.03e+007];
    case "copper"
        Mprop = [1, 0.999991, 5.8e+007];
    case "aluminium"
        Mprop = [1, 1.000022, 3.5d7];
    case "porcelain"
        Mprop = [5.7, 1, 1.e-15]; 
    case "Teflon" 
        Mprop = [2.5, 1, 0]; 
    case "MagSteel" 
        Mprop = [1, 2000, 2.e+006]; 
    case "aMagSteel" 
        Mprop = [1, 1, 5.5E6]; % for pipe, from https://www.mdpi.com/2076-3417/12/2/872
    case "soil" 
        Mprop = [1, 1, 1E-2];   
end

