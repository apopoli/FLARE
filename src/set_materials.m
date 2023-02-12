function [material] = set_materials(flag,ndom)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

material = strings(ndom,1);

switch flag
    case ('mesh_unit_circle')
        material(1) = 'air';
    case ('mesh_ex_regions')
        material(1) = 'air';
        material(2) = 'copper';
        material(3) = 'copper';
        material(4) = 'MagSteel';
    case ('mesh_corridor')
        material(1:3) = 'aluminium'; % phase_1, 2 e 3
        material(4) = 'air'; % OGW_1
        material(5) = 'air'; % Aria
        material(6) = 'soil'; % Terreno
        material(7) = 'aMagSteel'; % Pipe aMagSteel
        material(8) = 'air'; % Pipe centro
        material(9:12) = 'air'; % Mit 1, 2, 3, 4
end