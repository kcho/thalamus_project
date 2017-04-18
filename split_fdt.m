% https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=fsl;d3f65c4e.1405
% Load Matrix2
x = load('fdt_matrix2.dot');
M = full(spconvert(x));

% Load coordinate information to save results
addpath([getenv('FSLDIR') '/etc/matlab']);
[mask,~,scales] = read_avw('fdt_paths');
mask = 0*mask;

coord = load('tract_space_coords_for_fdt_matrix2')+1;
ind  = sub2ind(size(mask), coord(:,1), coord(:,2), coord(:,3));

for i = 1:size(M,1);
    mask3D = repmat(mask,1);
    mask3D(ind) = M(i,:);
    save_avw(mask3D,strcat(int2str(i), '.nii.gz'),'f',scales);
end
