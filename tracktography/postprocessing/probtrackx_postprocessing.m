% Load Matrix2
x = load('fdt_matrix2.dot');
M = full(spconvert(x));

% Load coordinate information to save results
addpath([getenv('FSLDIR') '/etc/matlab']);
[mask,~,scales] = read_avw('fdt_paths');

coord = load('tract_space_coords_for_fdt_matrix2')+1;
ind  = sub2ind(size(mask), coord(:,1), coord(:,2), coord(:,3));

% k-means clustering
% https://github.com/durai23/mri_zips/blob/master/probtrack/probtrackx_seg_kmeans.m

% Correlation coefficient 
cc=corrcoef(M');

% Number of clusters
k=8;
idx=kmeans(cc,k);
[~,~,j] = unique(idx);
mask(ind) = j;
save_avw(mask,'kmeans_segmentation','i',scales);
!fslcpgeom fdt_paths segmentation

%save 4d images
mask = 0*mask;
mask4D = repmat(mask, [1 1 1 size(M,1)]);
for i = 1:size(M,1);
    mask(ind) = M(i,:);
    mask4D(:,:,:,i) = mask;
    %disp(i)
end
save_avw(mask4D, 'merged_no_resample.nii.gz', 'f', scales);