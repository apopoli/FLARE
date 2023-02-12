% utility functions
num_regions = @(msh) length(unique(msh.TRIANGLES(:,4))); % number of regions within the domain
compar_mat = @(sparse_mat_1,sparse_mat_2) fprintf('Max difference is %g\n',max(max(abs(full(sparse_mat_1)-full(sparse_mat_2)))));
amplitude = @(x) (max(x)-min(x))/2;

addpath("src\");
addpath("sources\");
addpath("mesh\");

load("utils\colormaps.mat") % https://it.mathworks.com/matlabcentral/fileexchange/69470-custom-colormap