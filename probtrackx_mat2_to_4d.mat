% Load Matrix2
x = load('fdt_matrix2.dot');
M = full(spconvert(x));

% Load coordinate information to save results
addpath([getenv('FSLDIR') '/etc/matlab']);
[mask,~,scales] = read_avw('fdt_paths');
mask = 0*mask;

coord = load('tract_space_coords_for_fdt_matrix2')+1;
ind  = sub2ind(size(mask), coord(:,1), coord(:,2), coord(:,3));

% resample points
[new_x, new_y, new_z] = ...
    ndgrid(linspace(1, size(mask,1), size(mask,1)/3), ...
           linspace(1, size(mask,1), size(mask,1)/3), ...
           linspace(1, size(mask,1), size(mask,1)/3));

resample_mask = interp3(mask, new_x, new_y, new_z);

mask4D = repmat(resample_mask, [1 1 1 size(M,1)]);

for i = 1:size(M,1);
    mask(ind) = M(i,:);
    resample_M = interp3(mask, new_x, new_y, new_z);
    mask4D(:,:,:,i) = resample_M;
end
save_avw(mask4D, 'merged.nii.gz', 'f', [3 3 3 1]);
